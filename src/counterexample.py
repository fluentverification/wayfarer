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
			counterexamples.append(new_tail)
		else:
			traceback(state, new_tail)

def find_counterexamples(crn, number=1):
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
	pq.put((vass_priority(init_state, boundary, crn), init_state))
	while (not pq.empty()) and len(counterexamples) < number:
		curr_state = pq.get()
		if satisfies(curr_state, boundary):
			traceback(curr_state)
		priority = vass_priority(curr_state, boundary, crn)
		for _, next_state in get_transitions(curr_state, crn):
			if not next_state in backward_pointers:
				backward_pointers[next_state] = [curr_state]
			else:
				backward_pointers[next_state].append(curr_state)

