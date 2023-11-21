USE_CUDA=False

if not USE_CUDA:
	import numpy as np
else:
	import cupy as np

from distance import vass_distance
from crn import *

DONT_CARE = -1

class Subspace:
	def __init__(self, transitions):
		self.transitions = transitions
		basis_vectors = [t.vector for t in transitions]
		A = np.matrix(basis_vectors[0].append(basis_vectors[1:]), axis=1)
		self.P = A * (A.T * A) ** -1 * A.T
		self.rank = np.linalg.rank(self.P)

class State:
	subspaces = None
	target = None
	# A "mask" vector that gives us 1.0's in the species we care about and 0.0s in
	# the species we don't
	mask = None # (target > DONT_CARE).astype(float)
	init = None
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
		self.vec = vec
		self.adj = vec - State.init
		self.__compute_order()

	def __compute_order(self):
		'''
		Computes the order and epsilon vector of the state
		'''
		dist_to_target = self.norm(self.vec - State.target)
		if dist_to_target == 0:
			self.epsilon = [0.0]
			self.order = -1
			return 0
		self.epsilon = [dist_to_target]
		self.order = 0
		for s in reversed(State.subspaces):
			epsilon = self.norm(s.P * self.adj - self.adj)
			if epsilon == 0:
				return
			self.epsilon.insert(0, epsilon)
			self.order += 1

	def norm(self, vec):
		'''
		Computes the norm of a vector, only accounting for species in the CRN that
		aren't "Don't care" (value -1)
		'''
		return np.linalg.norm(np.multiply(vec, mask))

	def successors(self):
		'''
		Only returns the successors using the vectors in the dependency graph
		that get us closer to the target.
		'''
		succ = []
		# TODO: actually compute this rate using ALL transitions, not just the ones we use for successors
		total_outgoing_rate = 0.0
		subspace = State.subspaces[len(State.subspaces) - (self.index + 1)]
		for t in subspace.transitions:
			if t.enabled(self.vec)
				next_state = State(self.vec + t.vec_as_mat)
				rate = t.rate_finder(next_state.vec)
				# Due to the cycle-free nature of the dependency graph, we
				# can ignore successors with a higher distance if both the
				# current state and successor have order 0 (are in the last subspace)
				# Note: this only works if the last subspace has rank 1
				if subspace.rank == 1 and self.order == 0 and \
					next_state.epsilon[len(next_state.epsilon) - 1] > self.epsilon[len(self.epsilon) - 1]:
					continue
				succ.append((next_state, rate))
		return succ, total_outgoing_rate

	# Comparators. ONLY COMPARES THE ORDER AND THE LOWEST VALUE FOR EPSILON
	def __gt__(self, other):
		return self.order > other.order or self.epsilon[0] > other.epsilon[0]

	def __lt__(self, other):
		return self.order < other.order or self.epsilon[0] < other.epsilon[0]

	def __le__(self, other):
		return not self > other

	def __ge__(self, other):
		return not self < other

	def __eq__(self, other):
		'''
		Weak equal means are we equally close to the next subspace.
		Operator overloading done here because of priority queue, which
		possibly may use `=` operator.
		'''
		return self.order == other.order or self.epsilon[0] == other.epsilon[0]

	# Strong equality means we are the same state
	def strong_equal(self, other):
		return self.vec == other.vec

