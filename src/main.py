#!/usr/bin/env python3

from distance import *
from counterexample import *
from parser import *

import argparse
import time

if __name__=="__main__":
	parser = argparse.ArgumentParser(
		prog="wayfarer"
		, description="Wayfarer is a heuristic counterexample generator compatible with RAGTIMER files"
		, epilog="Developed at USU")
	parser.add_argument("-r", "--ragtimer", default=None)
	parser.add_argument("-n", "--number", default=3)
	args = parser.parse_args()
	if args.ragtimer is None:
		print("Missing args.")
		sys.exit(1)
	crn = parse_ragtimer(args.ragtimer)
	num = int(args.number)
	print("========================================================")
	print("Targeted Exploration (just distance)")
	print("========================================================")
	start_time = time.time()
	find_counterexamples(crn, number=num, print_when_done=True)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")
	print("========================================================")
	print("Targeted Exploration (Distance and angle)")
	print("========================================================")
	start_time = time.time()
	find_counterexamples(crn, number=num, print_when_done=True, include_flow_angle=True)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")
	print("========================================================")
	print("Random Exploration")
	print("========================================================")
	start_time = time.time()
	find_counterexamples_randomly(crn, number=num, print_when_done=True)
	end_time = time.time()
	print(f"Total time {end_time - start_time} s")
