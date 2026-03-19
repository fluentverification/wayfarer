import numpy as np
from scipy.linalg import null_space

import pulp
from pulp import LpProblem, LpMinimize, LpVariable, lpSum, LpInteger, LpStatus, LpStatusOptimal, value, LpSolutionInfeasible, LpSolutionUnbounded, LpSolutionNoSolutionFound

import time
from copy import deepcopy
from itertools import permutations

from crn import *
# from sbspc import Subspace

# TODO: are we sure we need both?
from fractions import Fraction

class Cycle:
	commute_depth = 3
	def __init__(self, ordered_reactions : list, crn: Crn):
		'''
	Creates an object which represents a cycle in the abstract.
		'''
		self.__crn = crn
		self.ordered_reactions = ordered_reactions.copy()
		self.__check_cycle_valid()
		self.__in_s0 = [r.in_s0 for r in ordered_reactions]
		self.is_commute_cycle = False
		# if not np.any(self.__in_s0):
		# 	raise Exception("Cycle must leave S0 (else it may be useless)!")
		self.is_orthocycle = np.all(self.__in_s0)
		if len(self.ordered_reactions) == 2:
			COMMUTE_DEPTH=Cycle.commute_depth
			self.ordered_reactions = (COMMUTE_DEPTH * [self.ordered_reactions[0]]) + (COMMUTE_DEPTH * [self.ordered_reactions[1]])
			self.is_commute_cycle = True
	def order_by_rate(self, state_vec) -> list:
		return sorted(self.ordered_reactions, key=lambda t: t.rate_finder(state_vec), reverse=True)

	def __check_cycle_valid(self):
		sum = self.ordered_reactions[0].vec_as_mat
		for r in self.ordered_reactions[1::]:
			sum += r.vec_as_mat
		assert (np.all(np.isclose(sum, 0)))

	def apply_cycle(self, state : np.matrix) -> list | None:
		'''
	TODO: I would like to apply cycles to every state individually, and then go through the list of
	newly created states and then after all of that is done add edges only to states which exist/have
	been created in the model.
		'''
		s = state
		new_states = []
		for r in self.ordered_reactions:
			s += r.vec_as_mat
			if np.any(s < 0):
				return None
			new_states.append(s)
		return new_states

	def __str__(self):
		return ' <-> '.join([t.name for t in self.ordered_reactions])

	def perms(self):
		# if len(self.ordered_reactions) == 2 or self.is_commute_cycle:
		if self.is_commute_cycle:
			return [self.ordered_reactions, self.ordered_reactions[::-1]]
		return list(permutations(self.ordered_reactions))
		# ps = [(perm, perm[::-1]) for perm in list(permutations(self.ordered_reactions))[::4]]
		# return [perm for perm_and_reverse in ps for perm in perm_and_reverse]


def get_commutable_transitions(crn : Crn, s0, _ss) -> list:
	transitions = []
	for t in crn.transitions:
		# TODO: need offset?
		# print(t.vec_as_mat.T)
		# If the transition's update vector is orthogonal to S0, then it is very likely commutable
		# print((s0.P * t.vec_as_mat).T)
		# print(np.isclose(s0.P * t.vec_as_mat, 0).T)
		# print(np.all(np.isclose(s0.P * t.vec_as_mat, 0).T))
		if np.all(np.isclose(s0.P * t.vec_as_mat, 0)):
			transitions.append(t)
	print(f"[INFO] Got {len(transitions)} trivially commutable transitions: {
            ','.join([t.name for t in transitions])}")
	return transitions


primes = [2]

def is_prime(i : int) -> bool:
	for j in range(2, i):
		if i % j == 0:
			return False
	return True

def init_primes_to(n : int):
	last_prime = primes[len(primes) - 1]
	for i in range(last_prime, n):
		if is_prime(i + 1):
			primes.append(i + 1)

def get_prime_factors(n : int) -> list:
	factors = []
	init_primes_to(n)
	for p in primes:
		if p == n:
			return [(p, 1)]
		elif p > n:
			break
		exponent : int = 0
		while n % p == 0:
			exponent += 1
			n //= p
		if exponent > 0:
			factors.append((p, exponent))
	return factors

def get_fraction(f : float, mx_denom=100) -> tuple:
	'''
In order to get past floating point weirdness, we have to put the float
into a string first.
	'''
	# assert(type(f) == Fraction)
	# return f.as_integer_ratio()
	return Fraction(f).limit_denominator(max_denominator=mx_denom).as_integer_ratio()

def minimal_scaling_factor(vec : np.matrix) -> int:
	'''
This does not work if the data type of the matrix is `float`
	'''
	primes = {}
	for elem in vec:
		n, d = get_fraction(elem[0, 0])
		factors = get_prime_factors(d)
		for p, e in factors:
			primes[p] = max(primes[p], e) if p in primes else e
	scale_factor = 1
	for prime, exponent in primes.items():
		scale_factor *= prime ** exponent
	return scale_factor

def get_nullvectors(R : np.matrix, atol=1e-13, rtol=0) -> list:
	'''
Gets the nullvectors of matrix R (where R has all positive integers)
and ensures that all of the nullvectors are also of type int.
	'''
	A = np.atleast_2d(R)
	u, s, vh = np.linalg.svd(A)
	# v = vh.T # This is all real valued so the hermetian is just the transpose
	rank = np.linalg.matrix_rank(A)
	# We can get the columns of v by using vh.tolist() which provides the rows of vh (i.e., the columns of v) in a list
	vcols = vh.tolist()
	# The last n - r columns of v are the nullspace basis
	ns = vcols[rank::]
	return ns
	# return null_space(R) # Todo: turn into list of columns

def get_cycles(crn: Crn, transitions: list, num: int = 3) -> list:
	matrix = np.column_stack([t.vec_as_mat.copy() for _, t in transitions])
	tran_to_idx = [idx for idx, _ in transitions]
	max_rate_const = np.max([t.rate_constant for _, t in transitions])
	weights = [max_rate_const * 1 / max(t.rate_constant, 1) for _, t in transitions]
	vecs = get_cycle_vectors(matrix, num, weights=weights)
	return cycles_from_cycle_vectors(vecs, crn, tran_to_idx)

def get_cycle_vectors(R : np.matrix, num=5, weights: list | None = None):
	'''
Any positive integer linear combination of the nullvectors of R are the cycles
of the graph. This gives us a set of `num` *reasonably small* cycle vectors.

The weights are how important it is to minimize each variable
	'''
	# TODO: add support for combinations beyond that
	print("[WARNING] Wayfarer only supports \"first level\" cycle detection currently. This means only linear combinations of null vectors with coefficients equal to 1 or 0")
	cycles = []
	start_time = time.time()
	# utilize PuLP's ILP engine to get potential cycle Parikh vectors
	problem = LpProblem("Cycle_Detection_Problem", sense=LpMinimize)
	solver = pulp.PULP_CBC_CMD(msg=False)

	# Create variables for integer null vectors with bounds 0-2
	n = R.shape[1]  # Number of columns
	x = [LpVariable(f'x{i}', lowBound=0, upBound=2, cat='Integer') for i in range(n)]

	if weights is not None:
		assert len(weights) == n

	while len(cycles) < num:
		# clear previous problem
		problem += lpSum(0)
		# Create constraints for null vector
		for row in R.tolist():
			problem += lpSum(row[j] * x[j] for j in range(n)) == 0

		# Minimize either L1 norm or weighted L1 norm (smaller cycles) -- only on the first iteration
		# if len(cycles) == 0:
		if weights is not None:
			problem += lpSum(weights[j] * x[j] for j in range(n))
		else:
			# There are no weights so we just minimize L1 Norm
			problem += lpSum(x)

		# Add a constraint to ensure x is not the zero vector
		problem += lpSum(x) >= 1  # At least one component must be greater than zero

		# Exclude previously found solutions
		for cycle in cycles:
			# Add constraints to ensure not equal to any previously found vector
			# At least one component is different
			problem += lpSum([x[i] - cycle[i] for i in range(n)]
			                 ) >= 1 or lpSum([cycle[i] - x[i] for i in range(n)]) >= 1

		# Solve the problem
		problem.solve(solver)

		# Check results
		if problem.status == LpStatusOptimal:
			cycle = [value(var) for var in x]
			cycles.append(cycle)
		else:
			break
	print(f"[INFO] Found {len(cycles)} cycles after {time.time() - start_time} s.")

	return cycles

def get_ordered_reactions(crn : Crn) -> list:
	'''
Orders the reactions by rate constant, returns a list of their
indexes in the crn's `transitions` member list
	'''
	# Yeah this is messy, but whatever
	sortable_transitions = [SortableTransition(t, 0) for t in crn.transitions]
	for i in range(len(sortable_transitions)):
		sortable_transitions[i].index = i
	sortable_transitions.sort(reverse=True)
	return [st.index for st in sortable_transitions]

def expand_vec(vec: np.matrix, tran_to_idx: list, num_transitions: int) -> list:
	'''
When we generate cycle Parikh vectors, we only provide the transitions we want included in the cycle.
As a result, we have to convert this into a parikh vector that is actually useful to the CRN (i.e.,
has slots for all transitions in the S-VAS, not just for the ones we want included in the VAS.
	'''
	veclist = [0 for _ in range(num_transitions)]
	for old_idx, new_idx in enumerate(tran_to_idx):
		veclist[new_idx] = vec[old_idx]
	return veclist

def cycles_from_cycle_vectors(vecs : list, crn : Crn, tran_to_idx: list) -> list:
	'''
Creates cycles from cycle vectors
	'''
	cycles = []
	sorted_transitions = get_ordered_reactions(crn)
	n = len(crn.transitions)
	# print(f"Parikh vector dimension: {n}")
	for v in vecs:
		# Rather than combinatorially expand to all possible
		# just get those with the most probable transitions
		# first, since they are most likely to be probable
		v_counter = expand_vec(v, tran_to_idx, n)
		transitions = []
		while not np.all([vi == 0 for vi in v_counter]):
			for idx in sorted_transitions:
				if v_counter[idx] > 0:
					transitions.append(deepcopy(crn.transitions[idx]))
					v_counter[idx] -= 1
		try:
			cycles.append(Cycle(transitions, crn))
		except Exception:
			pass

	print(f"[INFO] Found {len(cycles)} usable cycles:")
	for c in cycles:
		print(c)

	return cycles

def is_cyclable(transition, s0) -> bool:
	P = s0.P
	r = transition.vec_as_mat
	return np.isclose(P * r, r).all()
