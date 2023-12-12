USE_CUDA=False
VERIFY=True

# if not USE_CUDA:
import numpy as np
# else:
	# import cupy as np

# if VERIFY:
# from nagini_contracts.contracts import *
# from typing import List

from distance import vass_distance
from crn import *

DONT_CARE = -1

class Subspace:
	mask = None # (target > DONT_CARE).astype(float)
	# Type of elements in transitions: crn.Transition
	def __init__(self, transitions, excluded_transitions):
		'''
		Creates a subspace with a projection matrix and all that fun stuff
		transitions : the transitions forming the basis of the subspace
		excluded_transitions : all other transitions in the VASS
		'''
		# Requires(Forall(int, lambda i : Implies(i > 0 && i < len(transitions), len(transitions[i].vector) == len(transitions[i - 1].vector))))
		# Requires(Forall(int, lambda i : Implies(i > 0 && i < len(excluded_transitions), len(excluded_transitions[i].vector) == len(excluded_transitions[i - 1].vector))))
		# Requires(len(transitions) >= 1 && len(excluded_transitions) >= 1)
		# Requires(len(transitions[0].vector) == len(excluded_transitions[0].vector))
		self.transitions = transitions
		self.excluded_transitions = excluded_transitions
		basis_vectors = [np.matrix(t.vector).T for t in transitions]
		# print(basis_vectors[0])
		# TODO: make sure we're appending to the right axis
		A = np.column_stack(basis_vectors) # np.matrix(basis_vectors[0].append(basis_vectors[1:], axis=1))
		print(A)
		# Use the pseudoinverse since with rectangular matrices,
		# A.T * A is not guaranteed to be invertible. This SHOULD
		# still produce a valid projection matrix
		self.P = A * np.linalg.pinv(A.T * A) * A.T
		self.rank = np.linalg.matrix_rank(self.P)

	# @Pure
	def contains(self, other, test_vec): # -> bool:
		# Check if contains
		# Requires(type(State.init) == np.matrix)
		# Requires(len(test_vec) == len(Subspace.mask))
		# Ensures(Implies(Result(), self.dist(test_vec) >= other.dist(test_vec)))
		if self.rank < other.rank:
			return False
		return np.linalg.matrix_rank(np.block([self.P, other.P])) == self.rank

	# @Pure
	def norm(self, vec): # -> float:
		'''
		Computes the norm of a vector, only accounting for species in the CRN that
		aren't "Don't care" (value -1)
		'''
		# Requires(len(vec) == len(Subspace.mask))
		# Ensures(Result() >= 0.0)
		return float(np.linalg.norm(np.multiply(vec, Subspace.mask)))

	# @Pure
	def dist(self, vec): # -> float:
		# Requires(len(vec) == len(Subspace.mask))
		# Ensures(Result() >= 0.0)
		return self.norm(self.P * vec - vec)


class State:
	subspaces : list = []
	target : np.matrix = None
	# A "mask" vector that gives us 1.0's in the species we care about and 0.0s in
	# the species we don't

	init : np.matrix = None
	# @staticmethod
	def initialize_static_vars(crn, dep):
		State.subspaces = dep.create_subspaces(crn)
		State.init = np.matrix(crn.init_state).T
		State.target = np.matrix([b.to_num() for b in crn.boundary]).T

	def __init__(self, vec):
		'''
		Constructor for a new State element. Members within the State class:
		1. vec (type: np.matrix) : the actual vector representing the state values
		for a particular state in the state space.
		2. adj (type: np.matrix) : the ADJUSTED state. I.e., the state minus the
		initial state.
		3. order (type : int) : The number of subspaces we must pass through in order
		to reach the smallest subspace.
		4. epsilon (type list(int)) : The list of nonzero subspace distances, starting
		with the largest and working towards the smallest
		'''
		# Requires(type(State.init) == np.matrix)
		# Requires(type(State.target) == np.matrix)
		# Requires(type(State.subspaces) == List[Subspace])
		# # Requires(Forall(int, lambda i : Implies(i > 0 and i < len(State.subspaces), State.subspaces[i].contains(State.subspaces[i - 1], State.init))))
		# Ensures(self.order >= -1)
		# Ensures(len(self.epsilon) == len(State.subspaces) + 1)
		self.vec = vec
		self.adj = vec - State.init
		self.order : int = 0
		self.__compute_order()

	def __compute_order(self):
		'''
		Computes the order and epsilon vector of the state
		'''
		# Requires(type(State.init) == np.matrix)
		# Requires(type(State.target) == np.matrix)
		# Ensures(type(self.order) == int and self.order >= -1)
		# Ensures(type(self.epsilon == list[float]))
		# Ensures(len(self.epsilon) >= 1)
		# Ensures(Forall(int, lambda i : (Implies(i > 0 and i < len(State.subspaces), self.epsilon[i] >= self.epsilon[i - 1]))))
		dist_to_target = State.subspaces[0].norm(self.vec - State.target)
		if dist_to_target == 0.0:
			self.epsilon : list = [0.0]
			self.order = -1
			return
		self.epsilon = [dist_to_target]
		self.order = 0
		for s in reversed(State.subspaces):
			ep = s.dist(self.adj)
			if ep == 0:
				return
			self.epsilon.insert(0, epsilon)
			self.order += 1

	# @Pure
	def successors(self): # -> tuple:
		'''
		Only returns the successors using the vectors in the dependency graph
		that get us closer to the target.
		'''
		# Requires(type(State.init) == np.matrix)
		# Requires(type(State.target) == np.matrix)
		# Ensures(type(Result()) == tuple)
		# Ensures(len(Result()) == 2)
		# Ensures(type(Result()[0]) == list)
		# Ensures(type(Result()[1]) == float)
		# Ensures(Result()[1] >= 0.0)
		succ = []
		total_outgoing_rate = 0.0
		subspace = State.subspaces[len(State.subspaces) - (self.order + 1)]
		for t in subspace.transitions:
			if t.enabled(self.vec):
				next_state = State(self.vec)
				rate = t.rate_finder(next_state.vec)
				total_outgoing_rate += rate
				# Due to the cycle-free nature of the dependency graph, we
				# can ignore successors with a higher distance if both the
				# current state and successor have order 0 (are in the last subspace)
				# Note: this only works if the last subspace has rank 1
				if subspace.rank == 1 and self.order == 0 and \
					next_state.epsilon[len(next_state.epsilon) - 1] > self.epsilon[len(self.epsilon) - 1]:
					continue
				succ.append((next_state, rate))
		# Compute this rate using ALL transitions, not just the ones we use for successors
		for t in subspace.excluded_transitions:
			if t.enabled(self.vec):
				rate = t.rate_finder(self.vec)
				total_outgoing_rate += rate
		return succ, total_outgoing_rate

	# Comparators. ONLY COMPARES THE ORDER AND THE LOWEST VALUE FOR EPSILON
	# @Pure
	def __gt__(self, other):
		# Requires(len(self.epsilon) == len(other.epsilon))
		# Requires(len(self.epsilon) > 0)
		return self.order > other.order or self.epsilon[0] > other.epsilon[0]

	# @Pure
	def __lt__(self, other):
		# Requires(len(self.epsilon) == len(other.epsilon))
		# Requires(len(self.epsilon) > 0)
		return self.order < other.order or self.epsilon[0] < other.epsilon[0]

	# @Pure
	def __le__(self, other):
		# Requires(len(self.epsilon) == len(other.epsilon))
		# Requires(len(self.epsilon) > 0)
		return not self > other

	# @Pure
	def __ge__(self, other):
		# Requires(len(self.epsilon) == len(other.epsilon))
		# Requires(len(self.epsilon) > 0)
		return not self < other

	# @Pure
	def __eq__(self, other):
		'''
		Weak equal means are we equally close to the next subspace.
		Operator overloading done here because of priority queue, which
		possibly may use `=` operator.
		'''
		# Requires(len(self.epsilon) == len(other.epsilon))
		# Requires(len(self.epsilon) > 0)
		return self.order == other.order or self.epsilon[0] == other.epsilon[0]

	# Strong equality means we are the same state
	# @Pure
	def strong_equal(self, other):
		# Requires(len(self.vec) == len(other.vec))
		return self.vec == other.vec

