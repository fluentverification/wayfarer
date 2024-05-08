'''
Restrictive upper bound calculator

Requires:

 - No "don't-cares". We are looking for a specific state
 - No degridation transitions. All states cannot be reached twice.
'''

import numpy as np
# from stormpy import Rational

from crn import *

def get_rmatrix(crn):
	R = np.matrix(crn.transitions[0]).T
	for t in crn.transitions[1:]:
		col = np.matrix(t.vector).T
		R = np.append(R, col, axis=1)
	return R

def get_rhs(crn):
	s0 = np.matrix(crn.init_state).T
	s = []
	for b in crn.boundary:
		if b.bound_type != BoundTypes.EQUAL:
			raise Exception("All bound types must be BoundTypes.EQUAL")
		s.append(b.bound)
	s = np.matrix(s).T
	return np.append(s, -s0, axis=1)

def search_from_result(crn, sol, x=8):
	# If all of the transition rates are constant and all transitions are always enabled,
	# we can calculate the upper bound exactly
	if crn.all_rates_const and crn.all_trans_always_enabled:
		print(f"Good news! we can get an exact probability!")
		# TODO: calculate permutations without repetition and find the path probability of a single
		# path. Afterward, multiply the two to get the exact upper bound
	elif crn.all_rates_const:
		# TODO
	else:
		# Can we compute here?

def test_get_ubound(crn):
	R = get_rmatrix(crn)
	ss0 = get_rhs(crn)
	try:
		sol = np.linalg.solve(R, ss0)
	except LinAlgError:
		raise Exception("This transition system must not be singular!")
	# S = np.append(R, ss0, axis=1)


