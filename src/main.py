#!/usr/bin/env python3

from distance import *
from counterexample import *
from parser import *
from solver import *
from subspace import Subspace

import argparse
import time

store_traces = False

def basic_priority(filename, num):
	crn = parse_ragtimer(filename)

	print("========================================================")
	print("Targeted Exploration (just distance)")
	print("========================================================")
	start_time = time.time()
	find_counterexamples(crn, number=num, print_when_done=True)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")
	# print("========================================================")
	# print("Targeted Exploration (Distance and angle)")
	# print("========================================================")
	# start_time = time.time()
	# find_counterexamples(crn, number=num, print_when_done=True, include_flow_angle=True)
	# end_time = time.time()
	# print(f"Total time {end_time - start_time} s")

def random(filename, num):
	crn = parse_ragtimer(filename)
	print("========================================================")
	print("Random Exploration")
	print("========================================================")
	start_time = time.time()
	find_counterexamples_randomly(crn, number=num, print_when_done=True)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")

def subspace_priority(filename, num):
	dep, crn = parse_dependency_ragtimer(filename)
	print("========================================================")
	print("Targeted Exploration (Subspace)")
	print("========================================================")
	start_time = time.time()
	find_counterexamples_subsp(crn, dep, number=num, print_when_done=True, write_when_done=store_traces)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")

def subspace_priority_solver(filename, num, time_bound, agnostic=False, piped=False, all_expand=False, single_order=False, use_rate_const=False):
	dep, crn = parse_dependency_ragtimer(filename, agnostic=agnostic)
	print("========================================================")
	print("Targeted Exploration (Subspace - With Solver)")
	print("========================================================")
	if piped:
		print("Info: Norms during on-the-fly exploration will be scaled by the \"flow\" of the CRN via a \"piped\" matrix.")
		piped_matrix = create_piped(crn, use_rate_const)
		Subspace.initialize_piped(piped_matrix)
	start_time = time.time()
	min_probability_subsp(crn, dep, number=num, print_when_done=True, write_when_done=store_traces, time_bound=time_bound, expand_all_states=all_expand, single_order=single_order)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")

if __name__=="__main__":
	parser = argparse.ArgumentParser(
		prog="wayfarer"
		, description="wayfarer is a heuristic partial state space generator for lower bounds compatible with RAGTIMER files"
		, epilog="Developed at USU")
	parser.add_argument("-r", "--ragtimer", default=None,
			help="The name of the .ragtimer file to check")
	parser.add_argument("-n", "--number", default=3,
			help="The number of either counterexamples or satisfying states to find, if without or with solver respectively")
	parser.add_argument("-p", "--primitive", action="store_true",
			help="Primitive, single order distance with traceback")
	parser.add_argument("-s", "--subspace", action="store_true",
			help="Iterative subspace reduction with traceback")
	parser.add_argument("-S", "--subspace_with_solver", action="store_true",
			help="Iterative subspace reduction with CTMC analysis in StormPy")
	parser.add_argument("-V", "--solver", action="store_true",
			help="Single order priority with CTMC analysis in StormPy")
	parser.add_argument("-m", "--random", action="store_true",
			help="Random exploration")
	parser.add_argument("-t", "--traces", action="store_true",
			help="Store traces to a file") # Store traces to a file
	parser.add_argument("-T", "--time", default=None,
			help="If doing CTMC analysis, the time bound on the eventually property. "
			+ "If this is not provided, checks `P=? [ true U \"satisfy\" ]`")
	parser.add_argument("-a", "--agnostic", action="store_true",
			help="Initial-state agnostic. Will construct dependency graph with ALL reactions that produce or consume a species, even if enough of that species exists in the initial state (this may help increase probability bounds).")
	parser.add_argument("-R", "--random_freq", default=None,
			help="When set, the frequency of iterations in which a state's ENTIRE catalog of enabled reactions are expanded, rather than just those in the smallest subspace which contains it. This may cause it to \"jump\" to a path with either more or less probability than the current one.")
	parser.add_argument("-P", "--piped", action="store_true",
			help="Use a \"piped\" matrix based on the normalized update vectors of the available transitions scaled by their rates. Assumes constant rates as rates are evaluated at the initial state for the matrix")
	parser.add_argument("-e", "--expand_all", action="store_true",
			help="When expanding states, expand ALL transitions rather than just those which the dependency graph specifies will get the shortest traces. This may help find higher probability traces.")
	parser.add_argument("-Q", "--rate_constant", action="store_true",
			help="Use rate constant in piped method.")
	parser.add_argument("-F", "--rate_finder", default=None,
			help="If reaction rate ought to be calculated by something other than the standard rate constant method, a Python file may be provided here in order to perform these custom calculations. Will look for a function called rate_finder(state : arraylike, rate_constant : float, reaction_name : str) -> float")
	parser.add_argument("-U", "--upper", action="store_true",
			help="Also compute upper bound")
	parser.add_argument("-A", "--abstract", action="store_true",
			help="Abstract states from trivially commutable reactions")
	args = parser.parse_args()
	store_traces = args.traces
	State.abstract_states = args.abstract
	if args.ragtimer is None:
		print("Missing args.")
		sys.exit(1)
	num = int(args.number)

	if args.rate_finder is not None:
		parse_custom_rate_finder(args.rate_finder)

	SolverSettings.COMPUTE_UPPER_BOUND = args.upper
	if args.subspace:
		subspace_priority(args.ragtimer, num)

	if args.subspace_with_solver:
		t = None
		if args.time is not None and not args.time.isnumeric():
			print(f"Time bound {args.time} is invalid. Will ignore.")
		elif args.time is not None:
			t = int(args.time)
		subspace_priority_solver(args.ragtimer
						, num
						, time_bound=t
						, agnostic=args.agnostic
						, piped=args.piped
						, all_expand=args.expand_all
						, use_rate_const=args.rate_constant)

	if args.solver:
		t = None
		if args.time is not None and not args.time.isnumeric():
			print(f"Time bound {args.time} is invalid. Will ignore.")
		elif args.time is not None:
			t = int(args.time)
		subspace_priority_solver(args.ragtimer
						, num
						, time_bound=t
						, agnostic=args.agnostic
						, piped=args.piped
						, all_expand=args.expand_all
						, single_order=True
						, use_rate_const=args.rate_constant)

	if args.primitive:
		basic_priority(args.ragtimer, num)

	if args.random:
		random(args.ragtimer, num)

