#!/usr/bin/env python3

from distance import *
from counterexample import *
from parser import *

import argparse
import time

def basic_priority(filename, num):
	crn = parse_ragtimer(filename)

	# print("========================================================")
	# print("Targeted Exploration (just distance)")
	# print("========================================================")
	# start_time = time.time()
	# find_counterexamples(crn, number=num, print_when_done=True)
	# end_time = time.time()
	# print(f"Total time {end_time - start_time} s")
	print("========================================================")
	print("Targeted Exploration (Distance and angle)")
	print("========================================================")
	start_time = time.time()
	find_counterexamples(crn, number=num, print_when_done=True, include_flow_angle=True)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")

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
	find_counterexamples_subsp(crn, dep, number=num, print_when_done=True)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")

if __name__=="__main__":
	parser = argparse.ArgumentParser(
		prog="wayfarer"
		, description="Wayfarer is a heuristic counterexample generator compatible with RAGTIMER files"
		, epilog="Developed at USU")
	parser.add_argument("-r", "--ragtimer", default=None)
	parser.add_argument("-n", "--number", default=3)
	parser.add_argument("-p", "--primitive", action="store_true")
	parser.add_argument("-s", "--subspace", action="store_true")
	parser.add_argument("-m", "--random", action="store_true")
	args = parser.parse_args()
	if args.ragtimer is None:
		print("Missing args.")
		sys.exit(1)
	num = int(args.number)
	if args.subspace:
		subspace_priority(args.ragtimer, num)

	if args.primitive:
		basic_priority(args.ragtimer, num)

	if args.random:
		random(args.ragtimer, num)

