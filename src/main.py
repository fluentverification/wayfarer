#!/usr/bin/env python3

from distance import *
from counterexample import *
from parser import *
from solver import *

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

def subspace_priority_solver(filename, num, time_bound):
	dep, crn = parse_dependency_ragtimer(filename)
	print("========================================================")
	print("Targeted Exploration (Subspace - With Solver)")
	print("========================================================")
	start_time = time.time()
	min_probability_subsp(crn, dep, number=num, print_when_done=True, write_when_done=store_traces, time_bound=time_bound)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")

if __name__=="__main__":
	parser = argparse.ArgumentParser(
		prog="wayfarer"
		, description="Wayfarer is a heuristic counterexample generator compatible with RAGTIMER files"
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
	parser.add_argument("-m", "--random", action="store_true",
			help="Random exploration")
	parser.add_argument("-t", "--traces", action="store_true",
			help="Store traces to a file") # Store traces to a file
	parser.add_argument("-T", "--time", default=None,
			help="If doing CTMC analysis, the time bound on the eventually property. "
			+ "If this is not provided, checks `P=? [ true U \"satisfy\" ]`")
	args = parser.parse_args()
	store_traces = args.traces
	if args.ragtimer is None:
		print("Missing args.")
		sys.exit(1)
	num = int(args.number)
	if args.subspace:
		subspace_priority(args.ragtimer, num)

	if args.subspace_with_solver:
		t = None
		if args.time is not None and not args.time.isnumeric():
			print(f"Time bound {args.time} is invalid. Will ignore.")
		elif args.time is not None:
			t = int(args.time)
		subspace_priority_solver(args.ragtimer, num, time_bound=t)

	if args.primitive:
		basic_priority(args.ragtimer, num)

	if args.random:
		random(args.ragtimer, num)

