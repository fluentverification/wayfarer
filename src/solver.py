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

force_end_traceback=False
reaches = {}
counterexamples = []
num_satstates = 0
state_ids = {}
all_states = []

all_transitions = []

class Entry:
	def __init__(self, col : int, val : float):
		self.col = col
		self.val = val

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
			for entry in self.from_list[row]:
				col = entry.col
				val = entry.val
				matrix_builder.add_next_value(row, col, val)
		return matrix_builder

	def build(self):
		'''
		Creates a stormpy SparseMatrix
		'''
		matrix_builder = self.to_smb()
		return matrix_builder.build()

	def size(self):
		print(len(self.from_list), len(self.exit_rates))
		return len(self.from_list)

	def add_exit_rate(self, idx : int, rate : float):
		while len(self.exit_rates) <= idx:
			self.exit_rates.append(None)
		self.exit_rates[idx] = rate

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
	while (not pq.empty()) and num_satstates < number:
		num_explored += 1
		if num_explored % 20000 == 0:
			print(f"Explored {num_explored} states")
		curr_state_data = pq.get()
		# print(f"Exploring state with index {curr_state_data.idx}")
		curr_state = curr_state_data.vec
		if satisfies(curr_state, boundary):
			print(f"Found satisfying state {tuple(curr_state)}")
			num_satstates += 1
			sat_states.append(curr_state_data.idx)
			# We will create a self-loop later, so declare the total exit rate as 1.0
			matrixBuilder.add_exit_rate(curr_state_data.idx, 1.0)
			continue

		successors, total_rate = curr_state_data.successors()
		if len(successors) == 0:
			# Introduce a self-loop
			matrixBuilder.add_next_value(curr_state_data.idx, curr_state_data.idx, 1.0)
		matrixBuilder.add_exit_rate(curr_state_data.idx, total_rate)
		for s, rate in successors:
			next_state = s.vec
			# If the state is new, we explore it
			next_state_tuple = tuple(next_state)
			if next_state_tuple not in state_ids:
				all_states.append(s)
				# Assign new index
				state_ids[next_state_tuple] = last_index
				s.idx = last_index
				last_index += 1
				next_state_tuple = tuple(next_state)
				# Only explore new states
				pq.put(s)
			else:
				s = all_states[state_ids[next_state_tuple]]
			assert(s.idx is not None)
			# Place the transition in the matrix
			matrixBuilder.add_next_value(curr_state_data.idx, s.idx, rate)
	if print_when_done:
		print(f"Explored {num_explored} states. Found {num_satstates} satisfying states.")
	finalize_and_check(matrixBuilder, sat_states)

def finalize_and_check(matrixBuilder : RandomAccessSparseMatrixBuilder, satisfying_state_idxs : list):
	# First, connect all terminal states to absorbing
	global all_states
	for state in all_states[1::]:
		if state.perimeter:
			matrixBuilder.add_next_value(state.idx, 0, 1.0)
			matrixBuilder.add_exit_rate(state.idx, 1.0)
	matrix = matrixBuilder.build()
	labeling = StateLabeling(matrixBuilder.size())
	labeling.add_label("satisfy")
	for idx in satisfying_state_idxs:
		labeling.add_label_to_state("satisfy", idx)
	components = SparseModelComponents(matrix, labeling, {}, False)
	properties = [] # TODO
	exit_rates = [rate if rate is not None else 0.0 for rate in matrixBuilder.exit_rates]
	components.exit_rates = exit_rates
	print(f"Exit rates size = {len(exit_rates)}. Model size = {matrixBuilder.size()}")
	model = stormpy.storage.SparseCtmc(components)
	print(f"Matrix built (size {matrixBuilder.size()})")
