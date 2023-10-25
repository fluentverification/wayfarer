from distance import *
from crn import *

import queue

backward_pointers = {}
counterexamples = []

all_transitions = []

def traceback(c_state, tail=[]):
	'''
recursive algorithm to get all tracebacks to init state from current state
	'''
	# TODO: Dynamic programming version of this
	global backward_pointers
	global counterexamples
	# We can afford a copy of the tail on every call, since traceback is only called from
	# satisfying states and works BACKWARDS
	new_tail = tail.copy()
	new_tail.append(c_state)
	for state in backward_pointers[c_state]:
		if not state in backward_pointers:
			print("Found counterexample!")
			counterexamples.append(new_tail)
		else:
			traceback(state, new_tail)

def find_counterexamples(crn, number=1, print_when_done=False):
	global backward_pointers
	boundary = crn.boundary
	init_state = crn.init_state
	# According to https://docs.python.org/3/library/queue.html, the priority queue
	# is, by default, a min queue
	pq = queue.PriorityQueue()
	# TODO: according to the same docs, these queues are threadsafe, so
	# we can eventually multithread this
	global counterexamples
	curr_state = None
	state_priority = vass_priority(init_state, boundary, crn)
	print(f"State {init_state} has priority {state_priority}")
	pq.put((state_priority, tuple(init_state)))
	num_explored = 0
	while (not pq.empty()) and len(counterexamples) < number:
		num_explored += 1
		curr_state = pq.get()[1]
		if satisfies(curr_state, boundary):
			traceback(curr_state)
		# else:
			# print(f"State {curr_state} does not satisfy condition")
		for _, next_state in get_transitions(curr_state, crn):
			if not tuple(next_state) in backward_pointers:
				priority = vass_priority(next_state, boundary, crn)
				# print(f"State {next_state} has priority {priority}")
				pq.put((priority, tuple(next_state)))
				backward_pointers[tuple(next_state)] = [curr_state]
			else:
				backward_pointers[tuple(next_state)].append(curr_state)
	if print_when_done:
		print(f"Explored {num_explored} states")
		print_counterexamples()

def print_counterexamples():
	global counterexamples
	print(f"Finished finding {len(counterexamples)} counterexamples")
	for ce in counterexamples:
		print(ce)
