from distance import *
from crn import *

import queue
import random

TRACEBACK_BOUND=10
DESIRED_NUMBER_COUNTEREXAMPLES=2

force_end_traceback=False
backward_pointers = {}
reaches = {}
counterexamples = []
num_counterexamples = 0

all_transitions = []

def reset():
	# global DESIRED_NUMBER_COUNTEREXAMPLES
	global backward_pointers
	global reaches
	global counterexamples
	global num_counterexamples
	num_counterexamples = 0
	counterexamples = []
	reaches = {}
	backward_pointers = {}

def traceback(c_state, tail=[], tail_probability = 1.0):
	'''
recursive algorithm to get SOME (not all) tracebacks to init state from current state
	'''
	# TODO: Dynamic programming version of this
	global backward_pointers
	global counterexamples
	global DESIRED_NUMBER_COUNTEREXAMPLES
	global num_counterexamples
	global force_end_traceback
	if force_end_traceback or num_counterexamples >= DESIRED_NUMBER_COUNTEREXAMPLES + TRACEBACK_BOUND:
		return
	# If our tail probability is too small to handle or is zero we shouldn't continue exploring
	if tail_probability == 0.0:
		# We've probably found all of the counterexamples we can
		force_end_traceback=True
		return
	# We can afford a copy of the tail on every call, since traceback is only called from
	# satisfying states and works BACKWARDS
	new_tail = tail.copy()
	new_tail.append(c_state)
	for normalized_rate, state in backward_pointers[c_state]:
		if not state in backward_pointers:
			new_tail.append(state)
			counterexamples.append((tail_probability * normalized_rate, new_tail))
			num_counterexamples = len(counterexamples)
			# print(f"Found counterexample! Now we have {len(counterexamples)}")
		else:
			traceback(state, new_tail, tail_probability * normalized_rate)
	# print(f"ERROR: found no counterexample in traceback!")

def find_counterexamples(crn, number=1, print_when_done=False, include_flow_angle=False):
	reset()
	global DESIRED_NUMBER_COUNTEREXAMPLES
	global backward_pointers
	global num_counterexamples
	global force_end_traceback
	DESIRED_NUMBER_COUNTEREXAMPLES = number
	boundary = crn.boundary
	init_state = crn.init_state
	# According to https://docs.python.org/3/library/queue.html, the priority queue
	# is, by default, a min queue
	pq = queue.PriorityQueue()
	# TODO: according to the same docs, these queues are threadsafe, so
	# we can eventually multithread this
	curr_state = None
	state_priority = vass_priority(init_state, boundary, crn, include_flow_angle=include_flow_angle)
	# print(f"State {init_state} has priority {state_priority}")
	reaches[tuple(init_state)] = 1.0
	pq.put((state_priority, tuple(init_state)))
	num_explored = 0
	while (not pq.empty()) and num_counterexamples < number and not force_end_traceback:
		# print(num_counterexamples)
		num_explored += 1
		curr_priority, curr_state = pq.get()
		# print(f"Got state {curr_state} at priority {curr_priority}")
		if satisfies(curr_state, boundary):
			print(f"Found satisfying state {curr_state}")
			force_end_traceback = False
			traceback(curr_state)
		# else:
			# print(f"State {curr_state} does not satisfy condition")
		transitions = get_transitions(curr_state, crn)
		total_rate = 0.0
		# reach = reaches[curr_state]
		for rate, _ in transitions:
			total_rate += rate
		for rate, next_state in transitions:
			if not tuple(next_state) in backward_pointers:
				# print(f"State {next_state} has priority {priority}")
				next_state_tuple = tuple(next_state)
				# Only explore new states
				# if not next_state_tuple in reaches:
				# curr_reach = rate / total_rate * reach
				priority = vass_priority(next_state, boundary, crn, include_flow_angle=include_flow_angle) # , curr_reach)
				# reaches[next_state_tuple] = curr_reach
				pq.put((priority, next_state_tuple))
				backward_pointers[next_state_tuple] = [(rate / total_rate, curr_state)]
			else:
				backward_pointers[tuple(next_state)].append((rate / total_rate, curr_state))
	if print_when_done:
		print(f"Explored {num_explored} states")
		print_counterexamples()

def find_counterexamples_randomly(crn, number=1, print_when_done=False, trace_length=100):
	reset()
	global DESIRED_NUMBER_COUNTEREXAMPLES
	global backward_pointers
	global num_counterexamples
	global counterexamples
	ces_set = {}
	DESIRED_NUMBER_COUNTEREXAMPLES = number
	boundary = crn.boundary
	init_state = crn.init_state
	while num_counterexamples < number:
		ce = [tuple(init_state)]
		curr_state = init_state
		prob = 1.0
		for _ in range(trace_length):
			transitions = get_transitions(curr_state, crn)
			total_rate = 0.0
			for rate, _ in transitions:
				total_rate += rate
			chosen_transition = transitions[random.randint(0, len(transitions) - 1)]
			next_state = tuple(chosen_transition[1])
			next_rate = chosen_transition[0]
			prob *= next_rate / total_rate
			ce.append(next_state)
			if satisfies(next_state, boundary):
				# print(next_state, " satisfies")
				num_counterexamples += 1
				ce_tup = tuple(ce)
				# print(ce_tup)
				if not ce_tup in ces_set:
					counterexamples.append((prob, ce))
					ces_set[ce_tup] = True
					break
			else:
				# print(next_state, " does not satisfy")
				curr_state = next_state
	if print_when_done:
		print_counterexamples()


def print_counterexamples(show_entire_trace=False):
	global counterexamples
	print(f"Finished finding {len(counterexamples)} counterexamples")
	lower_bound = 0.0
	for est_prob, ce in counterexamples:
		lower_bound += est_prob
		print(f"Counterexample size {len(ce)} (esimated probability {est_prob})")
		# print(ce)
		if show_entire_trace:
			for i in range(len(ce)):
				state = ce[i]
				additional_message = ""
				if i == 0:
					additional_message = " (satisfying state)"
				elif i == len(ce) - 1:
					additional_message = " (initial state)"
				print(f"\tState: {state}{additional_message}")
	print(f"Total lower bound probability: {lower_bound}")
