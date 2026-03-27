from distance import *
from crn import *
from sbspc import *

import sys
import math

import queue
import random

from stormpy import SparseMatrixBuilder, StateLabeling, SparseModelComponents
import stormpy

class SolverSettings:
	DESIRED_NUMBER_STATES = 2
	ABSORBING_INDEX = 0
	PRINT_FREQUENCY = 100000
	COMPUTE_UPPER_BOUND = False
	CYCLE_REPEAT = 1

all_states = []
state_ids = {}

class Entry:
	def __init__(self, col : int, val : float):
		self.col : int = col
		self.val : float = val

	def __eq__(self, other):
		return other.col == self.col

	def __gt__(self, other):
		return self.col > other.col

	def __le__(self, other):
		return not self > other

	def __lt__(self, other):
		return self.col < other.col

	def __ge__(self, other):
		return not self < other

	def __str__(self):
		return f"{self.col}:{self.val}"

class RandomAccessSparseMatrixBuilder:
	'''
	A wrapper class for random entry into storm's sparse matrix builder
	'''

	def __init__(self):
		# Self loop for the absorbing state
		self.from_list = [[Entry(0, 1.0)]]
		self.exit_rates = [1.0]

	def add_next_value(self, row : int, col : int, val : float):
		while len(self.from_list) <= row:
			self.from_list.append([])
		if row == col:
			if len(self.from_list[row]) == 1:
				print(f"[Warning]: overwriting self-loop rate on state {row}")
				self.from_list[row][0].val = val
				return
			elif len(self.from_list[row]) > 1:
				raise Exception("Cannot add self-loop on CTMC with exit edges")
		self.from_list[row].append(Entry(col, val))

	def remove_self_loop(self, row: int):
		'''
	Removes a self-loop if it exists. If not, does nothing
		'''
		if len(self.from_list) <= row:
			return
		self.from_list[row] = [e for e in self.from_list[row] if e.col != row]
		# self.exit_rates[row] = None

	def clear_row(self, row: int):
		if len(self.from_list) <= row:
			return
		self.from_list[row] = []

	def remove_absorbing_edge(self, row: int):
		if len(self.from_list) <= row:
			return
		self.from_list[row] = [e for e in self.from_list[row] if e.col != 0]

	def subtract_from_abs(self, row: int, val_to_subtract: float):
		'''
	Subtracts from the rate going to the absorbing state for a particular state.
		'''
		abs_edges = [i for i, e in enumerate(self.from_list[row]) if e.col == 0]
		assert len(abs_edges) == 1
		self.from_list[row][abs_edges[0]].val -= val_to_subtract
		assert self.from_list[row][abs_edges[0]].val >= 0

	def has_entry(self, row: int, col: int) -> bool:
		if len(self.from_list) <= row:
			return False
		for entry in self.from_list[row]:
			if entry.col == col:
				return True
		return False

	def create_exit_rates(self) -> list:
		return [self.row_sum(i) for i in range(len(self.from_list))]

	def to_smb(self):
		'''
		Creates a stormpy.SparseMatrixBuilder
		'''
		matrix_builder = SparseMatrixBuilder()
		for row in range(len(self.from_list)):
			self.from_list[row].sort()
			if len(self.from_list[row]) == 0:
				# with open("model.tra", 'a') as f:
				# 	f.write(f"{row} {row} 1.0\n")
				matrix_builder.add_next_value(row, row, 1.0)
			for entry in self.from_list[row]:
				col = entry.col
				val = entry.val
				if row == col:
					if len(self.from_list[row]) != 1:
						raise Exception(f"State {row} should only have one edge: a self loop. Got" + \
							f"{len(self.from_list[row])} edges!\n{','.join([str(e) for e in self.from_list[row]])}")
					matrix_builder.add_next_value(row, col, 1.0)
					# with open("model.tra", 'a') as f:
					# 	f.write(f"{row} {col} 1.0\n")
					break
				# with open("model.tra", 'a') as f:
				# 	f.write(f"{row} {col} {val}\n")
				matrix_builder.add_next_value(row, col, val)
		return matrix_builder

	def deadlocks(self) -> list:
		def is_deadlock(row: list, state_idx: int):
			if len(row) > 1 or len(row) == 0:
				return False
			else:
				return row[0].col == state_idx
		return [i for i, e in enumerate(self.from_list) if is_deadlock(e, i)]

	def build(self):
		'''
		Creates a stormpy SparseMatrix
		'''
		matrix_builder = self.to_smb()
		return matrix_builder.build()

	def size(self):
		return len(self.from_list)

	def add_exit_rate(self, idx : int, rate : float):
		while len(self.exit_rates) <= idx:
			self.exit_rates.append(None)
		self.exit_rates[idx] = rate

	def row_sum(self, row: int) -> float:
		assert row < len(self.from_list)
		sm = 0.0
		for e in self.from_list[row]:
			sm += e.val
		return sm

	def assert_all_entries_correct(self):
		# assert(len(self.exit_rates) == len(self.from_list))
		for i in range(len(self.from_list)):
			# print(f"{i}: {self.exit_rates[i]}, {[str(entry) for entry in self.from_list[i]]}")
			if len(self.from_list[i]) == 0:
				# self.exit_rates[i] = None
				# assert (self.exit_rates[i] is None)
				continue
			max_entry = max(self.from_list[i])
			max_rate = max_entry.val
			# if self.exit_rates[i] is not None and not self.exit_rates[i] >= max_rate and not np.isclose(self.exit_rates[i] - max_rate, 0):
			# print(f"Error: {self.exit_rates[i]} < {max_rate} (state index {i})")
			# assert (self.exit_rates[i] is None or (self.exit_rates[i] >=
		        # max_rate or math.isclose(max_rate, self.exit_rates[i])))
			# TODO: remove this extra check
			# rs = self.row_sum(i)
			# if not np.isclose(self.exit_rates[i], rs):
			# 	self.exit_rates[i] = rs
			# 	# raise Exception(f"State {i} has differing exit rates and row sums! {
			# 	#                self.exit_rates[i]} vs {rs}")

def stotup(state) -> tuple:
	return tuple(map(lambda i: int(i), state))

def min_probability_subsp(crn, dep, number=1, print_when_done=False, write_when_done=False, time_bound=None, expand_all_states=False, single_order=False, cnc=False, cycle_count=3):
	global all_states
	global state_ids
	State.initialize_static_vars(crn, dep, single_order=single_order, cnc=cnc, cycle_count=cycle_count)
	state_ids = {}
	all_states = []
	# Add the absorbing state
	matrixBuilder = RandomAccessSparseMatrixBuilder()
	last_index = SolverSettings.ABSORBING_INDEX + 1
	all_states.append(None)
	# Other stuff
	boundary = crn.boundary
	sat_states = []
	# Min queue
	pq = queue.PriorityQueue()
	curr_state = None
	# Create and enqueue
	init_state = crn.init_state
	init_state = State(init_state, last_index)
	all_states.append(init_state)
	pq.put((init_state))
	last_index += 1
	# The number of explored and satisfying states
	num_satstates = 0
	num_explored = 0
	while (not pq.empty()) and num_satstates < number:
		num_explored += 1
		if num_explored % SolverSettings.PRINT_FREQUENCY == 0:
			print(f"Explored {num_explored} states. Have {num_satstates} satisfying")
		curr_state_data = pq.get()
		# print(f"Exploring state with index {curr_state_data.idx}")
		curr_state = curr_state_data.vec
		if satisfies(curr_state, boundary):
			# print(f"Found satisfying state {tuple(curr_state)}")
			num_satstates += 1
			sat_states.append(curr_state_data.idx)
			# We will create a self-loop later, so declare the total exit rate as 1.0
			# matrixBuilder.add_exit_rate(curr_state_data.idx, 1.0)
			matrixBuilder.add_next_value(curr_state_data.idx, curr_state_data.idx, 1.0)
			curr_state_data.perimeter = False
			continue

		# Total expanded rate: the rate of transitions we EXPANDED in the graph
		# Total full rate: the total rate of all POSSIBLE enabled transitions from this state.
		successors, total_expanded_rate = curr_state_data.successors(all_successors=expand_all_states, incl_rand=True)
		total_full_rate = curr_state_data.get_total_outgoing_rate()
		# print(total_full_rate, total_expanded_rate)
		assert (total_full_rate + 1e-5 >= total_expanded_rate)
		if len(successors) == 0:
			print("No successors")
			# Introduce a self-loop
			matrixBuilder.add_next_value(curr_state_data.idx, curr_state_data.idx, 1.0)
			continue
		# If this is true there are some transitions we didn't expand that we must lead
		# to the absorbing state. We do this since we only take reactions in that subspace
		if total_full_rate > total_expanded_rate:
			if not cnc:
				matrixBuilder.add_next_value(curr_state_data.idx, 0, total_full_rate - total_expanded_rate)
		# matrixBuilder.add_exit_rate(curr_state_data.idx, total_full_rate)
		for s, rate in successors:
			next_state = s.vec
			# If the state is new, we explore it
			next_state_tuple = stotup(next_state)
			if next_state_tuple not in state_ids:
				next_index = len(all_states)
				all_states.append(s)
				assert (last_index == next_index)
				# Assign new index
				state_ids[next_state_tuple] = last_index
				s.idx = last_index
				last_index += 1
				next_state_tuple = stotup(next_state)
				# Only explore new states
				pq.put(s)
			# If this state already exists, use the state data we already have
			else:
				# update reachability for re-explored states
				s_old = s
				s = all_states[state_ids[next_state_tuple]]
				s.reach += s_old.reach
			assert (s.idx is not None)
			# Place the transition in the matrix
			if not cnc:
				matrixBuilder.add_next_value(curr_state_data.idx, s.idx, rate)
	if print_when_done:
		print(f"Explored {len(matrixBuilder.from_list)} states" + \
			f"(expanded {num_explored}). Found {num_satstates} satisfying states.\n")
	if num_satstates == 0:
		print(f"Could not find any satisfying states!")
		return
	if cnc:
		assert last_index == len(all_states)
		# Make an intelligent guess if necessary on the times to repeat cycles
		if SolverSettings.CYCLE_REPEAT is None:
			# If the seed state space is less than 10k states, perform cycle repetition
			# at a rate of 5. If less than 20k, perform it at a rate of 3. Otherwise, it
			# is too large of a state space for this heuristic to work well
			if num_explored <= 4000:
				if len(State.cycles) < 3:
					SolverSettings.CYCLE_REPEAT = 10
				else:
					SolverSettings.CYCLE_REPEAT = 7
			elif num_explored <= 10000:
				SolverSettings.CYCLE_REPEAT = 5
			elif num_explored <= 20000:
				SolverSettings.CYCLE_REPEAT = 3
			else:
				SolverSettings.CYCLE_REPEAT = 1 # No repetition
			print(f"[INFO] Cycle repetition not specified. Choosing based on model. Repetition factor: {SolverSettings.CYCLE_REPEAT}")
		apply_cycles(matrixBuilder, crn, last_index, sat_states)
	sanity_check()
	finalize_and_check(matrixBuilder, sat_states, time_bound, crn)

def apply_cycles(matrixBuilder, crn, next_available_idx, sat_indecies):
	if len(State.cycles) == 0:
		print("Cannot apply cycles! No cycles exist!")
		return
	print(f"Applying {len(State.cycles)} cycles...", end="", flush=True)
	start_time = time.time()
	global all_states
	global state_ids
	assert next_available_idx == len(all_states)

	# First we will apply all of the cycles, creating internal connections, and then create connections
	# between states that have been newly created if there is a one-step transition between them.
	repeat_cycles = [c for c in State.cycles if c.is_commute_cycle or len(c.ordered_reactions) < 5]
	last_start_idx = 1
	for i in range(SolverSettings.CYCLE_REPEAT):
		if i != 0:
			print("repeating...", end="", flush=True)
			# Only apply the commute cycles on the next iterations to reduce internal complexity
			next_available_idx = apply_specific_cycles(repeat_cycles, crn, sat_indecies, start_idx=last_start_idx)
		else:
			last_start_idx = next_available_idx
			next_available_idx = apply_specific_cycles(State.cycles, crn, sat_indecies)

	# Re-construct edges
	for state in all_states[1::]:
		assert next_available_idx == len(all_states)
		# state.perimeter = False
		state_id = state.idx
		state.perimeter = False
		assert state_id != 0
		# We will re-build the list of transitions here
		matrixBuilder.clear_row(state_id)
		total_full_rate = state.get_total_outgoing_rate()
		successors, total_exit_rate = state.successors(True, all_successors=True)
		# states not expanded will go to the absorbing state
		rate_to_abs = total_full_rate - total_exit_rate
		for stup, rate in successors:
			if stup in state_ids:
				next_idx = state_ids[stup]
				assert next_idx != 0
				matrixBuilder.add_next_value(state.idx, next_idx, rate)
			# If the successor is a satisfying state we should add it anyway
			elif satisfies(stup, crn.boundary):
				# Add a new state
				next_state = State(np.matrix(stup).T, next_available_idx, need_compute_order=False)
				next_state.perimeter = False
				state_ids[stup] = next_state.idx
				next_available_idx += 1
				all_states.append(next_state)
				matrixBuilder.add_next_value(state.idx, next_state.idx, rate)
				matrixBuilder.add_next_value(next_state.idx, next_state.idx, 1.0)
				sat_indecies.append(next_state.idx)
			else:
				rate_to_abs += rate
		if rate_to_abs > 0.0:
			matrixBuilder.add_next_value(state.idx, 0, rate_to_abs)

	print(f"...finished after {time.time() - start_time} seconds.")

def apply_specific_cycles(cycles, crn, sat_indecies, start_idx = 1) -> int:
	'''
Applies a specific list of cycles to the state graph and returns the new total state count. The start index 
is the first state index to start applying cycles to. Generally this is the initial state (index 1).
	'''
	global all_states
	global state_ids
	next_available_idx = len(all_states)
	for state in all_states[start_idx::]:
		# Do not apply cycles to satisfying states
		if satisfies(state.vec, crn.boundary):
			continue
		# Do not apply cycles to perimeter states (we do not want to expand them anyway)
		if state.perimeter and SolverSettings.CYCLE_REPEAT == 1:
			continue
		# We need to iterate over the states first, then the cycles.
		for cycle in cycles:
			for directed_cycle in cycle.perms():
				# We can apply the cycle forward and backward. We apply forward first.
				cur_state = state
				cur_state_idx = state.idx
				for t in directed_cycle:
					assert next_available_idx == len(all_states)
					# if cur_state.perimeter:
					# 	break
					if not t.enabled(cur_state.vecm):
						break
					# Get the next state and see if it's new or not
					next_state_vec = cur_state.vec + t.vector
					satisfying = satisfies(next_state_vec, crn.boundary)
					if np.any(next_state_vec < 0.0):
						break

					nsvt = stotup(next_state_vec)
					# print(f"Next state tuple: {nsvt}")
					if nsvt in state_ids:
						# State is not new
						next_idx = state_ids[nsvt]
						assert cur_state_idx != next_idx
						cur_state = all_states[next_idx]
						cur_state_idx = next_idx
						if satisfying or cur_state.perimeter:
							break
					else:
						# State is new, so we have to add it. We only have to actually add the edge if the previous state was already
						# in the graph, since otherwise finalize_and_check() will take care of that
						# matrixBuilder.add_next_value(cur_state_idx, next_available_idx, rate)
						next_state = State(next_state_vec, next_available_idx, need_compute_order=False)
						state_ids[nsvt] = next_state.idx

						next_available_idx += 1
						all_states.append(next_state)
						cur_state = next_state
						cur_state_idx = next_state.idx
						if satisfying:
							# NOTE: since it is a perimeter state, its index will be added to satisfying_state_idxs in finalize_and_check()
							# satisfying_state_idxs.append(next_available_idx)
							# We do not need to continue down this cycle
							sat_indecies.append(next_state.idx)
							break
	return next_available_idx


# This can become a lemma when we eventually use Nagini to verify this
def sanity_check():
	print("Performing sanity check...", end="", flush=True)
	global all_states
	global state_ids
	# Check our indecies
	idx = 0
	for state in all_states:
		assert (state is None or state.idx == idx)
		idx += 1
	# TODO: remove this intensive check when we've found the bug
	# for _, sid in state_ids.items():
	# 	assert sid < len(all_states)
	print("done.")

def finalize_and_check(matrixBuilder : RandomAccessSparseMatrixBuilder, satisfying_state_idxs : list, time_bound : int, crn : Crn = None):
	global state_ids
	# First, connect all terminal states to absorbing
	global all_states
	# NOTE: in the paper, we flush the queue, however here, we go through all states and connect all PERIMETER
	# states to the absorbing, which is the same thing.
	num_perim_satstates = 0
	next_available_idx = len(all_states)
	for state in all_states[1::]:
		if state.perimeter:
			if satisfies(state.vec, crn.boundary):
				# print(f"Found satisfying state {tuple(curr_state)}")
				num_perim_satstates += 1
				satisfying_state_idxs.append(state.idx)
				# We will create a self-loop later, so declare the total exit rate as 1.0
				# matrixBuilder.add_exit_rate(state.idx, 1.0)
				matrixBuilder.add_next_value(state.idx, state.idx, 1.0)
				# deadlock_idxs.append(state.idx)
				state.perimeter = False
				continue
			# Expand the state and create transitions ONLY TO EXISTING STATES
			successors, total_exit_rate = state.successors(True, all_successors=True)
			total_full_rate = state.get_total_outgoing_rate()
			# states not expanded will go to the absorbing state
			rate_to_abs = total_full_rate - total_exit_rate
			for s, rate in successors:
				stup = s  # tuple(s.vec)
				if stup in state_ids:
					next_idx = state_ids[stup]
					matrixBuilder.add_next_value(state.idx, next_idx, rate)
				elif satisfies(stup, crn.boundary):
					# We can easily just add this state.
					num_perim_satstates += 1
					next_state = State(np.matrix(stup).T, next_available_idx, need_compute_order=False)
					state_ids[stup] = next_state.idx
					next_available_idx += 1
					all_states.append(next_state)
					next_state.perimeter = True
					matrixBuilder.add_next_value(state.idx, next_state.idx, rate)
					matrixBuilder.add_next_value(next_state.idx, next_state.idx, 1.0)
					satisfying_state_idxs.append(next_state.idx)

				else:
					rate_to_abs += rate
			if rate_to_abs > 0.0:
				matrixBuilder.add_next_value(state.idx, 0, rate_to_abs)
			# matrixBuilder.add_exit_rate(state.idx, total_full_rate)
	if num_perim_satstates > 0:
		print(f"We found an additional {num_perim_satstates} satisfying states in the perimeter state indecies!")
	deadlock_idxs = matrixBuilder.deadlocks()
	matrix = matrixBuilder.build()
	matrixBuilder.assert_all_entries_correct()
	print(f"{matrixBuilder.size()} ?= {len(all_states)}")
	assert matrixBuilder.size() == len(all_states)
	assert matrix.nr_rows == matrixBuilder.size()
	print(f"Shape (Row, Col): {matrix.nr_rows}, {matrix.nr_columns}")
	assert matrix.nr_columns == matrix.nr_rows
	print(f"Number of rows in matrix: {matrix.nr_rows}")
	labeling = StateLabeling(matrixBuilder.size())
	label_file = [[] for _ in range(len(all_states))]
	# with open("model.lab", 'a') as lf:
		# lf.write("0=\"init\" 1=\"satisfy\" 2=\"absorbing\" 3=\"deadlock\"\n")
	# Add initial state labeling
	labeling.add_label("init")
	labeling.add_label_to_state("init", 1)
	label_file[1].append(0)
	labeling.add_label("satisfy")
	labeling.add_label("absorbing")
	labeling.add_label_to_state("absorbing", 0)
	label_file[0].append(2)
	for idx in satisfying_state_idxs:
		assert idx < matrix.nr_rows
		labeling.add_label_to_state("satisfy", idx)
		label_file[idx].append(1)
	# Add the deadlock state indexes
	labeling.add_label("deadlock")
	for idx in deadlock_idxs:
		assert idx < matrix.nr_rows
		labeling.add_label_to_state("deadlock", idx)
		label_file[idx].append(3)
	# with open("model.lab", 'a') as lf:
		# for i, r in enumerate(label_file):
			# lf.write(f"{i}: ")
			# lf.write(' '.join([str(ri) for ri in r]))
			# lf.write("\n")
	components = SparseModelComponents(matrix, labeling, {}, rate_transitions=True)
	prop_bound = "" if time_bound is None else f"[0, {time_bound}]"
	chk_property = f"P=? [ true U{prop_bound} \"satisfy\" ]"
	# [rate if rate is not None else 1.0 for rate in matrixBuilder.exit_rates]
	exit_rates = matrixBuilder.create_exit_rates()
	assert len(exit_rates) == matrix.nr_rows
	components.exit_rates = exit_rates
	# print(f"Exit rates size = {len(exit_rates)}. Model size = {matrixBuilder.size()}")
	model = stormpy.storage.SparseCtmc(components)
	print(model)
	print(f"Matrix built (size {matrixBuilder.size()})")
	print(f"Checking model with formula `{chk_property}`")
	prop = stormpy.parse_properties(chk_property)[0]  # stormpy.Property("Lower Bound", )
	env = stormpy.Environment()
	env.solver_environment.native_solver_environment.precision = stormpy.Rational(1e-25)
	start_time = time.time()
	result = stormpy.check_model_sparse(model, prop, only_initial_states=True)
	print(f"Model checking took {time.time() - start_time} seconds.")
	print(f"Pmin = {result.at(1)}")
	assert (result.min + 1e-6 >= 0.0 and result.max <= 1.0 + 1e-6)
	if SolverSettings.COMPUTE_UPPER_BOUND:
		# Upper bound propert
		chk_property_upper = f"P=? [ true U{prop_bound} \"satisfy\" | \"absorbing\" ]"
		# print(f"Checking upper bound with property `{chk_property_upper}`")
		prop_upper = stormpy.parse_properties(chk_property_upper)[0]
		result_upper = stormpy.check_model_sparse(model, prop_upper, only_initial_states=True)
		print(f"Pmax = {result_upper.at(1)}")
