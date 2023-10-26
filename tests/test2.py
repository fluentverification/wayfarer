#!/usr/bin/env python3

from crn import *
from counterexample import *

import time

crn = Crn([
	# Transition system
	Transition(
		[-1, 0, 0, 0]
		, lambda state : state[0] > 50
		, lambda state : 0.01
		)
	, Transition(
		[10, 0, 0, 0]
		, lambda state : state[0] < 1000
		, lambda state : 3.3
		)
	, Transition(
		[0, 1, 0, 0]
		, lambda state : state[1] < 100
		, lambda state : 1
		)
	, Transition(
		[0, -1, 0, 0]
		, lambda state : state[1] > 1
		, lambda state : 0.05
		)
	, Transition(
		[0, 0, 10, 0]
		, lambda state : state[2] < 100
		, lambda state : 4
		)
	, Transition(
		[0, 0, -2, 0]
		, lambda state : state[2] > 2
		, lambda state : 0.03
		)
	, Transition(
		[0, 0, 0, 2]
		, lambda state : state[2] < 100
		, lambda state : 4
		)
	, Transition(
		[0, 0, 0, -1]
		, lambda state : state[2] > 2
		, lambda state : 0.003
		)
	]
	# The satisfying condition
	# At least 6.852600255650947e-75
	, [Bound(200, BoundTypes.GREATER_THAN), Bound(2, BoundTypes.GREATER_THAN), Bound(98, BoundTypes.EQUAL), Bound(98, BoundTypes.EQUAL)]
	# At least 1.7743981340125083e-91
	# , [Bound(400, BoundTypes.GREATER_THAN), Bound(10, BoundTypes.GREATER_THAN), Bound(98, BoundTypes.EQUAL), Bound(98, BoundTypes.EQUAL)]
	# At least 6.518817323914444e-103
	# , [Bound(600, BoundTypes.GREATER_THAN), Bound(10, BoundTypes.GREATER_THAN), Bound(98, BoundTypes.EQUAL), Bound(98, BoundTypes.EQUAL)]
	# , lambda state : return state[0] > 650 and state[1] > 20 and state[3] == 8
	# The initial state
	, np.array([0, 0, 0, 0])
)

if __name__=="__main__":
	print("========================================================")
	print("Targeted Exploration")
	print("========================================================")
	start_time = time.time()
	find_counterexamples(crn, number=10, print_when_done=True)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")
	print("========================================================")
	print("Random Exploration")
	print("========================================================")
	start_time = time.time()
	find_counterexamples_randomly(crn, number=10, print_when_done=True)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")
