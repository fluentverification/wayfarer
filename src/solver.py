from distance import *
from crn import *
from subspace import *

import sys

import queue
import random

from stormpy import SparseMatrixBuilder, StateLabeling, SparseModelComponents
import stormpy

sys.setrecursionlimit(sys.getrecursionlimit() * 50)

DESIRED_NUMBER_STATES=2
ABSORBING_INDEX=0

PRINT_FREQUENCY=100000

force_end_traceback=False
reaches = {}
counterexamples = []
num_satstates = 0
state_ids = {}
all_states = []

all_transitions = []

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
		self.from_list[row].append(Entry(col, val))

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
					assert(len(self.from_list[row]) == 1)
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
				assert(self.exit_rates[i] is None)
				continue
			max_entry = max(self.from_list[i])
			max_rate = max_entry.val
			assert(self.exit_rates[i] >= max_rate)

def reset():
	# global DESIRED_NUMBER_COUNTEREXAMPLES
	global reaches
	global counterexamples
	global num_satstates
	global state_ids
	global all_states
	num_satstates = 0
	counterexamples = []
	reaches = {}
	state_ids = {} # tuple : int
	all_states = []

def min_probability_subsp(crn, dep, number=1, print_when_done=False, write_when_done=False):
	reset()
	State.initialize_static_vars(crn, dep)
	global DESIRED_NUMBER_STATES
	global backward_pointers
	global num_satstates
	global force_end_traceback
	global ABSORBING_INDEX
	global all_states
	# Add the absorbing state
	matrixBuilder = RandomAccessSparseMatrixBuilder()
	last_index = ABSORBING_INDEX + 1
	all_states.append(None)
	# Other stuff
	DESIRED_NUMBER_STATES = number
	boundary = crn.boundary
	num_explored = 0
	sat_states = []
	# Min queue
	pq = queue.PriorityQueue()
	curr_state = None
	init_state = crn.init_state
	reaches[tuple(init_state)] = 1.0
	init_state = State(init_state, last_index)
	all_states.append(init_state)
	pq.put((init_state))
	last_index += 1
	deadlock_idxs = [0]
	while (not pq.empty()) and num_satstates < number:
		num_explored += 1
		if num_explored % PRINT_FREQUENCY == 0:
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
		successors, total_expanded_rate = curr_state_data.successors()
		total_full_rate = curr_state_data.get_total_outgoing_rate()
		assert(total_full_rate >= total_expanded_rate)
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
				assert(last_index == next_index)
				# Assign new index
				state_ids[next_state_tuple] = last_index
				s.idx = last_index
				last_index += 1
				next_state_tuple = tuple(next_state)
				# Only explore new states
				pq.put(s)
			# If this state already exists, use the state data we already have
			else:
				s = all_states[state_ids[next_state_tuple]]
			assert(s.idx is not None)
			# Place the transition in the matrix
			matrixBuilder.add_next_value(curr_state_data.idx, s.idx, rate)
	if print_when_done:
		print(f"Explored {len(matrixBuilder.from_list)} states (expanded {num_explored}). Found {num_satstates} satisfying states.")
	sanity_check()
	finalize_and_check(matrixBuilder, sat_states, deadlock_idxs)


# This can become a lemma when we eventually use Nagini to verify this
def sanity_check():
	print("Performing sanity check...", end="")
	global all_states
	# Check our indecies
	idx = 0
	for state in all_states:
		assert(state is None or state.idx == idx)
		idx += 1
	print("done.")

def finalize_and_check(matrixBuilder : RandomAccessSparseMatrixBuilder, satisfying_state_idxs : list, deadlock_idxs : list):
	# First, connect all terminal states to absorbing
	global all_states
	for state in all_states[1::]:
		if state.perimeter:
			matrixBuilder.add_next_value(state.idx, 0, 1.0)
			matrixBuilder.add_exit_rate(state.idx, 1.0)
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
	chk_property = "P=? [ true U \"satisfy\" ]"
	exit_rates = [rate if rate is not None else 1.0 for rate in matrixBuilder.exit_rates]
	components.exit_rates = exit_rates
	# print(f"Exit rates size = {len(exit_rates)}. Model size = {matrixBuilder.size()}")
	model = stormpy.storage.SparseCtmc(components)
	print(model)
	print(f"Matrix built (size {matrixBuilder.size()})")
	print("Checking model...")
	prop = stormpy.parse_properties(chk_property)[0] # stormpy.Property("Lower Bound", )
	result = stormpy.check_model_sparse(model, prop, only_initial_states=True)
	assert(result.min >= 0.0 and result.max <= 1.0)
	print(f"Pmin = {result.at(1)}")
