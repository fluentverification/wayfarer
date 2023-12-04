USE_CUDA=False
VERIFY=True

if not USE_CUDA:
	import numpy as np
else:
	import cupy as np

# if VERIFY:
from nagini_contracts.contracts import *

from distance import vass_distance
from crn import *

DONT_CARE = -1

class Subspace:
	mask = None # (target > DONT_CARE).astype(float)
	def __init__(self, transitions : Transition, excluded_transitions : list(Transition)):
		'''
		Creates a subspace with a projection matrix and all that fun stuff
		transitions : the transitions forming the basis of the subspace
		excluded_transitions : all other transitions in the VASS
		'''
		self.transitions = transitions
		self.excluded_transitions = excluded_transitions
		basis_vectors = [t.vector for t in transitions]
		A = np.matrix(basis_vectors[0].append(basis_vectors[1:]), axis=1)
		self.P = A * (A.T * A) ** -1 * A.T
		self.rank = np.linalg.rank(self.P)

	# @Pure
	def contains(self, other : Subspace, test_vec : np.matrix = State.init) -> bool:
		# Check if contains
		Requires(len(test_vec) == len(Subspace.mask))
		Ensures(Implies(Result(), self.dist(test_vec) >= other.dist(test_vec)))
		if self.rank < other.rank:
			return False
		return np.linalg.rank(np.block(self.P, other.P)) == self.rank

	# @Pure
	def norm(self, vec : np.matrix) -> float:
		'''
		Computes the norm of a vector, only accounting for species in the CRN that
		aren't "Don't care" (value -1)
		'''
		Requires(len(vec) == len(mask))
		Ensures(Result() >= 0.0)
		return np.linalg.norm(np.multiply(vec, mask))

	@Pure
	def dist(self, vec) -> float:
		Requires(len(vec) == len(mask))
		Ensures(Result() >= 0.0)
		return self.norm(s.P * vec - vec)


class State:
	subspaces = None
	target = None
	# A "mask" vector that gives us 1.0's in the species we care about and 0.0s in
	# the species we don't

	init = None
	def __init__(self, vec : np.matrix):
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
		Requires(type(State.subspaces) == list(Subspace))
		Requires(
			Forall(int, lambda i : (Implies(i > 0 and i < len(State.subspaces),
			State.subspaces[i].contains[State.subspaces[i - 1]))))
		Ensures(self.order >= -1)
		Ensures(len(self.epsilon) == len(subspaces) + 1)
		self.vec = vec
		self.adj = vec - State.init
		self.order = None
		self.__compute_order()

	def __compute_order(self):
		'''
		Computes the order and epsilon vector of the state
		'''
		Ensures(type(self.order) == int and self.order >= -1)
		Ensures(type(self.epsilon == list[float]))
		Ensures(len(self.epsilon) >= 1)
		Ensures(Forall(int, lambda i : (Implies(i > 0 and i < len(State.subspaces), self.epsilon[i] >= self.epsilon[i - 1]))))
		dist_to_target = State.subspaces[0].norm(self.vec - State.target)
		if dist_to_target == 0.0:
			self.epsilon = [0.0]
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
	def successors(self) -> (list, float):
		'''
		Only returns the successors using the vectors in the dependency graph
		that get us closer to the target.
		'''
		Ensures(type(Result()) == tuple)
		Ensures(len(Result()) == 2)
		Ensures(type(Result()[0]) == list)
		Ensures(type(Result()[1]) == float)
		Ensures(Result()[1] >= 0.0)
		succ = []
		total_outgoing_rate = 0.0
		subspace = State.subspaces[len(State.subspaces) - (self.index + 1)]
		for t in subspace.transitions:
			if t.enabled(self.vec)
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
	@Pure
	def __gt__(self, other):
		Requires(len(self.epsilon) == len(other.epsilon))
		Requires(len(self.epsilon) > 0)
		return self.order > other.order or self.epsilon[0] > other.epsilon[0]

	@Pure
	def __lt__(self, other):
		Requires(len(self.epsilon) == len(other.epsilon))
		Requires(len(self.epsilon) > 0)
		return self.order < other.order or self.epsilon[0] < other.epsilon[0]

	@Pure
	def __le__(self, other):
		Requires(len(self.epsilon) == len(other.epsilon))
		Requires(len(self.epsilon) > 0)
		return not self > other

	@Pure
	def __ge__(self, other):
		Requires(len(self.epsilon) == len(other.epsilon))
		Requires(len(self.epsilon) > 0)
		return not self < other

	@Pure
	def __eq__(self, other):
		'''
		Weak equal means are we equally close to the next subspace.
		Operator overloading done here because of priority queue, which
		possibly may use `=` operator.
		'''
		Requires(len(self.epsilon) == len(other.epsilon))
		Requires(len(self.epsilon) > 0)
		return self.order == other.order or self.epsilon[0] == other.epsilon[0]

	# Strong equality means we are the same state
	@Pure
	def strong_equal(self, other):
		Requires(len(self.vec) == len(other.vec))
		return self.vec == other.vec

