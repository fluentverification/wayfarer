import numpy as np

from crn import *
from subspace import *

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

	def apply_cycle(self, state : np.matrix) -> list:
		s = state
		new_states = []
		for r in self.__ordered_reactions:
			s += r.vec_as_mat
			new_states.append(s)
		return new_states

def get_commutable_transitions(crn : Crn, s0 : Subspace, ss : Subspace) -> list:
	transitions = []
	for t in crn.transitions:
		# TODO: need offset?
		if np.all(np.is_close(s0.P * t.vec_as_mat, 0)).all() and np.all(np.isclose(ss.P * t.vec_as_mat, t.vec_as_mat)):
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

def get_fraction(f : Fraction) -> tuple:
	'''
Requires type Fraction
	
Old: ~~In order to get past floating point weirdness, we have to put the float
into a string first.~~
	'''
	assert(type(f) == Fraction)
	return f.as_integer_ratio()
	# return Fraction(str(f)).as_integer_ratio()

def minimal_scaling_factor(vec : np.matrix) -> int:
	'''
This does not work if the data type of the matrix is `float`
	'''
	primes = {}
	for elem in vec:
		n, d = get_fraction(elem[0, 0])
		print(f"{n} / {d}")
		factors = get_prime_factors(d)
		for p, e in factors:
			primes[p] = max(primes[p], e) if p in primes else e
	scale_factor = 1
	for prime, exponent in primes.items():
		scale_factor *= prime ** exponent
	return scale_factor
	

