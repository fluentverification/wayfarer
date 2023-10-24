#!/usr/bin/env python3

from crn import *
from counterexample import *

crn = Crn([
	# Transition system
	Transition(
		[-10, 0, 0]
		, lambda state : state[0] > 50
		, lambda state : 1.7
		)
	, Transition(
		[40, 0, 0]
		, lambda state : state[0] < 1000
		, lambda state : 5.7
		)
	, Transition(
		[0, 1, 0]
		, lambda state : state[1] < 100
		, lambda state : 3
		)
	, Transition(
		[0, -1, 0]
		, lambda state : state[1] > 1
		, lambda state : 0.5
		)
	, Transition(
		[0, 0, 10]
		, lambda state : state[2] < 100
		, lambda state : 4
		)
	, Transition(
		[0, 0, -2]
		, lambda state : state[2] > 2
		, lambda state : 0.3
		)
	]
	# The satisfying condition
	, [Bound(650, BoundTypes.GREATER_THAN), Bound(20, BoundTypes.GREATER_THAN), Bound(3, BoundTypes.EQUAL)]
	# , lambda state : return state[0] > 650 and state[1] > 20 and state[3] == 8
	# The initial state
	, np.matrix([0, 0, 0])
)

if __name__=="__main__":
	find_counterexamples(crn, number=3)
