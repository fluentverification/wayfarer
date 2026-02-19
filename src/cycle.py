import numpy as np
from scipy.linalg import null_space

from pulp import LpProblem, LpMinimize, LpVariable, lpSum, LpInteger, LpStatus, LpStatusOptimal, value, LpSolutionInfeasible, LpSolutionUnbounded, LpSolutionNoSolutionFound
import time

from crn import *
# from sbspc import Subspace

# TODO: are we sure we need both?
from fractions import Fraction

class Cycle:
	def __init__(self, ordered_reactions : list):
		'''
	Creates an object
		'''
		self.__ordered_reactions = ordered_reactions
		self.__check_cycle_valid()
		self.__in_s0 = [r.in_s0 for r in ordered_reactions]
		if not np.any(self.__in_s0):
			raise Exception("Cycle must leave S0 (else it may be useless)!")
		self.is_orthocycle = np.all(self.__in_s0)

	def __check_cycle_valid(self) -> bool:
		sum = self.__ordered_reactions[0].vec_as_mat
		for r in self.__ordered_reactions[1::]:
			sum += r.vec_as_mat
		assert(np.all(np.is_close(sum, 0)))

	def apply_cycle(self, state : np.matrix) -> list | None:
		s = state
		new_states = []
		for r in self.__ordered_reactions:
			s += r.vec_as_mat
			if np.any(s < 0):
				return None
			new_states.append(s)
		return new_states

def get_commutable_transitions(crn : Crn, s0, ss) -> list:
	transitions = []
	for t in crn.transitions:
		# TODO: need offset?
		if np.all(np.isclose(s0.P * t.vec_as_mat, 0)).all() and np.all(np.isclose(ss.P * t.vec_as_mat, t.vec_as_mat)):
			transitions.append(t)
	return t

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
			n /= p
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

def get_cycles(crn: Crn, transitions: list, num: int = 5) -> list:
	matrix = np.column_stack([t.vec_as_mat for t in transitions])
	vecs = get_cycle_vectors(matrix, num)
	print(f"Cycle vectors: {vecs}")
	return cycles_from_cycle_vectors(vecs, crn)

def get_cycle_vectors(R : np.matrix, num=5):
	'''
Any positive integer linear combination of the nullvectors of R are the cycles
of the graph. This gives us a set of `num` *reasonably small* cycle vectors.
	'''
	print(R)
	# TODO: add support for combinations beyond that
	print("[WARNING] Wayfarer only supports \"first level\" cycle detection currently. This means only linear combinations of null vectors with coefficients equal to 1 or 0")
	cycles = []
	start_time = time.time()
	# utilize PuLP's ILP engine to get potential cycle Parikh vectors
	problem = LpProblem("Cycle_Detection_Problem", sense=LpMinimize)

	# Create variables for integer null vectors with bounds 0-1
	n = R.shape[1] # Number of columns
	x = [LpVariable(f'x{i}', lowBound=0, upBound=2, cat='Integer') for i in range(n)]

	while len(cycles) < num:
		# clear previous problem
		problem += lpSum(0)
		# Create constraints for null vector
		for row in R.tolist():
			problem += lpSum(row[j] * x[j] for j in range(n)) == 0

		# Minimize L1 norm (smaller cycles)
		# problem += lpSum(x)

		# Add a constraint to ensure x is not the zero vector
		problem += lpSum(x) >= 1  # At least one component must be greater than zero

		# Exclude previously found solutions
		for cycle in cycles:
			# Add constraints to ensure not equal to any previously found vector
			problem += lpSum([x[i] - cycle[i] for i in range(n)]) != 1  # At least one component is different
			problem += lpSum([cycle[i] - x[i] for i in range(n)]) != 1  # Ensure no matching components

		# Solve the problem
		problem.solve()

		# Check results
		print(LpStatus[problem.status])
		if problem.status == LpStatusOptimal:
			cycle = [value(var) for var in x]
			# print(f"Found cycle {cycle}")
			cycles.append(cycle)
		else:
			print("Could not find any more cycles!")
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

def cycles_from_cycle_vectors(vecs : list, crn : Crn) -> list:
	'''
Creates cycles from cycle vectors
	'''
	cycles = []
	sorted_transitions = get_ordered_reactions(crn)
	print(sorted_transitions)
	for v in vecs:
		print(v)
		# Rather than combinatorially expand to all possible
		# just get those with the most probable transitions
		# first, since they are most likely to be probable
		v_counter = v.copy()
		transitions = []
		while not np.all(v_counter == 0):
			for idx in sorted_transitions:
				transitions.append(crn.transitions[idx])
				v_counter[idx] -= 1
		cycles.append(Cycle(transitions, crn))

	print(f"Found cycles: {cycles}")

	return cycles

def is_cyclable(transition, s0) -> bool:
	P = s0.P
	r = transition.vec_as_mat
	return np.isclose(P * r, r).all()
