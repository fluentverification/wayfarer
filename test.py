#!/usr/bin/env python3
import time
import subprocess
import re
import sys
import os

models = [
	("models/2react/2react", 100)
	, ("models/6react/6react", 100)
	, ("models/8react/8react", 20)
	, ("models/12react/12react", 10)
]

def getResultFromOutput(output : str) -> float:
	pMin : float = 0.0
	for line in output.split("\n"):
		if line.startswith("Pmin = "):
			pMin = float(line.replace("Pmin = ", ""))
			break
	return pMin

def getFoundStatesFromOutput(output : str) -> int:
	foundStates : int = 0
	matches = re.findall(r'.+Found ([0-9]+) satisfying states.', output)
	print(len(matches))
	print(output)
	assert(len(matches) == 1)
	return matches[0]

def runTest(fileBaseName : str, timeBound : int, num : int = 2, agnostic : bool = False, piped : bool = False, sop : bool = False) -> None:
	print(f"Testing Wayfarer with file {fileBaseName} searching for {num} sat states")
	modelFile = f"{fileBaseName}.sm"
	ragtimerFile = f"{fileBaseName}.ragtimer"
	propFile = f"{fileBaseName}.csl"
	# Test wayfarer
	startTime = time.time()
	cmdArray = None
	if sop:
		print("Testing single order priority")
		cmdArray = ["./wayfarer", "-r", ragtimerFile, "-V", "-T", f"{timeBound}", "-n", f"{num}"]
	else:
		cmdArray = ["./wayfarer", "-r", ragtimerFile, "-S", "-T", f"{timeBound}", "-n", f"{num}"]
	if agnostic:
		cmdArray.append("--agnostic")
	if piped:
		cmdArray.append("--piped")

	p = subprocess.run(cmdArray, capture_output=True)
	endTime = time.time()
	results = p.stdout.decode("utf-8")
	return getFoundStatesFromOutput(results), getResultFromOutput(results), endTime - startTime

if __name__=="__main__":
	postfix="/"
	if "WAYFARER_TEST_POSTFIX" in os.environ:
		postfix = f"{os.environ['WAYFARER_TEST_POSTFIX']}/"
		print(f"Putting results in results/{postfix}")
	agnostic = "-a" in sys.argv or "--agnostic" in sys.argv
	sop = "-s" in sys.argv or "--sop" in sys.argv
	if agnostic:
		print("Initial state agnostic")
	if sop:
		print("Single order priority")
	piped = "--piped" in sys.argv
	if piped:
		print("Piped method")
	for baseName, tBound in models:
		with open(f"results/{postfix}{baseName.split('/').pop()}.csv", "w") as f:
			f.write("Number Sat Requested,Number Sat Found,Probability Min,Duration\n")
			for p in range(14):
				num = 2 ** p
			#for p in range(4):
			#	num = 10 ** p
				numFound, pMin, duration = runTest(baseName, tBound, num=num, agnostic=agnostic, piped=piped, sop=sop)
				print(f"Took {duration} s and found pMin={pMin} ({numFound} actual found satisfying states)")
				f.write(f"{num},{numFound},{pMin},{duration}\n")
