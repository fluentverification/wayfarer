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
			assert (self.from_list[row] == 0.0)
		self.from_list[row].append(Entry(col, val))

	def has_entry(self, row: int, col: int) -> bool:
		if len(self.from_list) <= row:
			return False
		for entry in self.from_list[row]:
			if entry.col == col:
				return True
		return False

	def to_smb(self):
		'''
		Creates a stormpy.SparseMatrixBuilder
		'''
		matrix_builder = SparseMatrixBuilder()
		for row in range(len(self.from_list)):
			self.from_list[row].sort()
			if len(self.from_list[row]) == 0:
				matrix_builder.add_next_value(row, row, 1.0)
			for entry in self.from_list[row]:
				col = entry.col
				val = entry.val
				if row == col:
					assert (len(self.from_list[row]) == 1)
					matrix_builder.add_next_value(row, col, 1.0)
					break
				matrix_builder.add_next_value(row, col, val)
		return matrix_builder

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

	def assert_all_entries_correct(self):
		# assert(len(self.exit_rates) == len(self.from_list))
		for i in range(len(self.from_list)):
			# print(f"{i}: {self.exit_rates[i]}, {[str(entry) for entry in self.from_list[i]]}")
			if len(self.from_list[i]) == 0:
				self.exit_rates[i] = None
				assert (self.exit_rates[i] is None)
				continue
			max_entry = max(self.from_list[i])
			max_rate = max_entry.val
			if self.exit_rates[i] is not None and not self.exit_rates[i] >= max_rate:
				print(f"Error: {self.exit_rates[i]} < {max_rate} (state index {i})")
			assert (self.exit_rates[i] is None or (self.exit_rates[i] >=
			        max_rate or math.isclose(max_rate, self.exit_rates[i])))

def min_probability_subsp(crn, dep, number=1, print_when_done=False, write_when_done=False, time_bound=None, expand_all_states=False, single_order=False, cnc=False):
	global all_states
	global state_ids
	State.initialize_static_vars(crn, dep, single_order=single_order, cnc=cnc)
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
	deadlock_idxs = [0]
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
			matrixBuilder.add_exit_rate(curr_state_data.idx, 1.0)
			matrixBuilder.add_next_value(curr_state_data.idx, curr_state_data.idx, 1.0)
			deadlock_idxs.append(curr_state_data.idx)
			curr_state_data.perimeter = False
			continue

		# Total expanded rate: the rate of transitions we EXPANDED in the graph
		# Total full rate: the total rate of all POSSIBLE enabled transitions from this state.
		successors, total_expanded_rate = curr_state_data.successors(all_successors=expand_all_states)
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
			matrixBuilder.add_next_value(curr_state_data.idx, 0, total_full_rate - total_expanded_rate)
		matrixBuilder.add_exit_rate(curr_state_data.idx, total_full_rate)
		for s, rate in successors:
			next_state = s.vec
			# If the state is new, we explore it
			next_state_tuple = tuple(next_state)
			if next_state_tuple not in state_ids:
				next_index = len(all_states)
				all_states.append(s)
				assert (last_index == next_index)
				# Assign new index
				state_ids[next_state_tuple] = last_index
				s.idx = last_index
				last_index += 1
				next_state_tuple = tuple(next_state)
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
			matrixBuilder.add_next_value(curr_state_data.idx, s.idx, rate)
	if print_when_done:
		print(f"Explored {len(matrixBuilder.from_list)} states (expanded {
		      num_explored}). Found {num_satstates} satisfying states.")
	if num_satstates == 0:
		print(f"Could not find any satisfying states!")
		return
	apply_cycles(matrixBuilder, sat_states, deadlock_idxs, crn, last_index)
	sanity_check()
	finalize_and_check(matrixBuilder, sat_states, deadlock_idxs, time_bound, crn)

def apply_cycles(matrixBuilder, satisfying_state_idxs, deadlock_idxs, crn, next_available_idx):
	if len(State.cycles) == 0:
		print("Cannot apply cycles! No cycles exist!")
		return
	print(f"Applying {len(State.cycles)} cycles...", end="")
	global all_states
	global state_ids
	assert next_available_idx == len(all_states)
	new_cycle_states = []
	# First we will apply all of the cycles, creating internal connections, and then create connections
	# between states that have been newly created if there is a one-step transition between them.
	for state in all_states[1::]:
		# We need to iterate over the states first, then the cycles.
		for cycle in State.cycles:
			# We can apply the cycle forward and backward. We apply forward first.
			cur_state = state
			cur_state_idx = state.idx
			last_state_in_graph = True
			for t in cycle.ordered_reactions:
				if not t.enabled(cur_state.vecm):
					break
				# Get the next state and see if it's new or not
				next_state_vec = cur_state.vec + t.vector
				if np.any(next_state_vec < 0.0):
					break
				rate = t.rate_finder(cur_state.vec)
				nsvt = tuple(next_state_vec)
				if nsvt in state_ids:
					# State is not new, just add it to the matrix builder
					next_idx = state_ids[nsvt]
					if matrixBuilder.has_entry(cur_state_idx, next_idx):
						# No need to add the entry since it already exists
						continue
					if cur_state_idx == next_idx:
						print(f"Warning: Wanted self-loop on index {next_idx}")
					matrixBuilder.add_next_value(cur_state_idx, next_idx, rate)
					cur_state = all_states[next_idx]
					cur_state_idx = next_idx
					last_state_in_graph = True
				else:
					# State is new, so we have to add it. We only have to actually add the edge if the previous state was already
					# in the graph, since otherwise finalize_and_check() will take care of that
					if last_state_in_graph:
						matrixBuilder.add_next_value(cur_state_idx, next_available_idx, rate)
					next_state = State(next_state_vec, next_available_idx, need_compute_order=False)
					next_available_idx += 1
					next_state.perimeter = True
					all_states.append(next_state)
					new_cycle_states.append(next_state)
					last_state_in_graph = False
			cur_state = state
			cur_state_idx = state.idx
			last_state_in_graph = True
			for t in cycle.ordered_reactions[::-1]:
				if not t.enabled(cur_state.vecm):
					break
				# Get the next state and see if it's new or not
				next_state_vec = cur_state.vec + t.vector
				if np.any(next_state_vec < 0.0):
					break
				rate = t.rate_finder(cur_state.vec)
				nsvt = tuple(next_state_vec)
				if nsvt in state_ids:
					# State is not new, just add it to the matrix builder
					next_idx = state_ids[nsvt]
					if matrixBuilder.has_entry(cur_state_idx, next_idx):
						# No need to add the entry since it already exists
						continue
					if cur_state_idx == next_idx:
						print(f"Warning: Wanted self-loop on index {next_idx}")
					matrixBuilder.add_next_value(cur_state_idx, next_idx, rate)
					cur_state = all_states[next_idx]
					cur_state_idx = next_idx
					last_state_in_graph = True
				else:
					# State is new, so we have to add it. We only have to actually add the edge if the previous state was already
					# in the graph, since otherwise finalize_and_check() will take care of that
					if last_state_in_graph:
						matrixBuilder.add_next_value(cur_state_idx, next_available_idx, rate)

					next_state = State(next_state_vec, next_available_idx, need_compute_order=False)
					next_available_idx += 1
					next_state.perimeter = True
					all_states.append(next_state)
					new_cycle_states.append(next_state)
					last_state_in_graph = False
	print("...finished")
	# # Next, we need to expand successors for all of these newly created state
	# for state in new_cycle_states:
	# 	successors, total_expanded_rate = state.successors(all_successors=True)
	# 	for s, rate in successors:
	# 		next_state_vec = s.vec
	# 		ns = tuple(next_state_vec)
	# 		if ns not in state_ids:
	# 			# We should not create this state or give it an ID. However, we have to handle the rate going to it.
	# 			# The rate going to that state should go to the absorbing state, however, it is included in total_expanded_rate
	# 			# in finalize_and_check(). So, we have to account for it. The way we do that is kind of "hackey": we SUBTRACT the
	# 			# rate from state.total_cycle_exit
	# 			state.total_cycle_exit -= rate
	# 			continue
	# 		else:
	# 			next_state_id = state_ids[ns]

	# We have to account for the rate we have taken away from
	# TODO: Gotta figure out what to do. Should we just leave the states floating and have them be taken care of by finalize_and_check()

# This can become a lemma when we eventually use Nagini to verify this
def sanity_check():
	print("Performing sanity check...", end="")
	global all_states
	# Check our indecies
	idx = 0
	for state in all_states:
		assert (state is None or state.idx == idx)
		idx += 1
	print("done.")

def finalize_and_check(matrixBuilder : RandomAccessSparseMatrixBuilder, satisfying_state_idxs : list, deadlock_idxs : list, time_bound : int, crn : Crn = None):
	global state_ids
	# First, connect all terminal states to absorbing
	global all_states
	# NOTE: in the paper, we flush the queue, however here, we go through all states and connect all PERIMETER
	# states to the absorbing, which is the same thing.
	num_perim_satstates = 0
	for state in all_states[1::]:
		if state.perimeter:
			if satisfies(state.vec, crn.boundary):
				# print(f"Found satisfying state {tuple(curr_state)}")
				num_perim_satstates += 1
				satisfying_state_idxs.append(state.idx)
				# We will create a self-loop later, so declare the total exit rate as 1.0
				matrixBuilder.add_exit_rate(state.idx, 1.0)
				matrixBuilder.add_next_value(state.idx, state.idx, 1.0)
				deadlock_idxs.append(state.idx)
				state.perimeter = False
				continue
			# else:
			# 	if len(State.commutable_transitions) > 0:
			# 		pass  # TODO
			# 	if len(State.orthocycles) > 0:
			# 		pass  # TODO
			# 	if len(State.non_orthocycles) > 0:
			# 		pass  # TODO
			# Expand the state and create transitions ONLY TO EXISTING STATES
			successors, total_exit_rate = state.successors(True)
			total_full_rate = state.get_total_outgoing_rate()
			# states not expanded will go to the absorbing state
			rate_to_abs = total_full_rate - total_exit_rate
			for s, rate in successors:
				stup = s  # tuple(s.vec)
				if stup in state_ids:
					next_idx = state_ids[stup]
					matrixBuilder.add_next_value(state.idx, next_idx, rate)
				else:
					rate_to_abs += rate
			if rate_to_abs > 0.0:
				matrixBuilder.add_next_value(state.idx, 0, rate_to_abs)
			matrixBuilder.add_exit_rate(state.idx, total_full_rate)
	if num_perim_satstates > 0:
		print(f"We found an additional {
		      num_perim_satstates} satisfying states in the perimeter state indecies!")
	matrix = matrixBuilder.build()
	matrixBuilder.assert_all_entries_correct()
	labeling = StateLabeling(matrixBuilder.size())
	# Add initial state labeling
	labeling.add_label("init")
	labeling.add_label_to_state("init", 1)
	labeling.add_label("satisfy")
	labeling.add_label("absorbing")
	labeling.add_label_to_state("absorbing", 0)
	for idx in satisfying_state_idxs:
		labeling.add_label_to_state("satisfy", idx)
	# Add the deadlock state indexes
	labeling.add_label("deadlock")
	for idx in deadlock_idxs:
		labeling.add_label_to_state("deadlock", idx)
	components = SparseModelComponents(matrix, labeling, {}, rate_transitions=True)
	prop_bound = "" if time_bound is None else f"[0, {time_bound}]"
	chk_property = f"P=? [ true U{prop_bound} \"satisfy\" ]"
	exit_rates = [rate if rate is not None else 1.0 for rate in matrixBuilder.exit_rates]
	components.exit_rates = exit_rates
	# print(f"Exit rates size = {len(exit_rates)}. Model size = {matrixBuilder.size()}")
	model = stormpy.storage.SparseCtmc(components)
	print(model)
	print(f"Matrix built (size {matrixBuilder.size()})")
	print(f"Checking model with formula `{chk_property}`")
	prop = stormpy.parse_properties(chk_property)[0]  # stormpy.Property("Lower Bound", )
	env = stormpy.Environment()
	env.solver_environment.native_solver_environment.precision = stormpy.Rational(1e-100)
	result = stormpy.check_model_sparse(model, prop, only_initial_states=True)
	print(f"Pmin = {result.at(1)}")
	assert (result.min + 1e-6 >= 0.0 and result.max <= 1.0 + 1e-6)
	if SolverSettings.COMPUTE_UPPER_BOUND:
		# Upper bound propert
		chk_property_upper = f"P=? [ true U{prop_bound} \"satisfy\" | \"absorbing\" ]"
		# print(f"Checking upper bound with property `{chk_property_upper}`")
		prop_upper = stormpy.parse_properties(chk_property_upper)[0]
		result_upper = stormpy.check_model_sparse(model, prop_upper, only_initial_states=True)
		print(f"Pmax = {result_upper.at(1)}")
