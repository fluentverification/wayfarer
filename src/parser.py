from crn import *
from distance import Bound, BoundTypes
from dependency import *

import sys

def check_valid_identifier(identifier):
	'''
Checks to see if the reactant identifier is valid.
Does not yet cover all cases
	'''
	if ' ' in identifier or '.' in identifier or identifier.isdigit():
		print(f"[PARSE ERROR] '{identifier}' is an invalid identifier")
		sys.exit(1)

def create_species_idxes(reactants):
	species_idxes = {}
	for i in range(len(reactants)):
		reactant = reactants[i].replace("\n", "")
		# print(f"Adding reactant {reactant}")
		check_valid_identifier(reactant)
		species_idxes[reactant] = i
	return species_idxes

def create_bound(bound_text):
	bound_int = int(bound_text)
	if bound_int == -1:
		return Bound(0, BoundTypes.DONT_CARE)
	return Bound(bound_int, BoundTypes.EQUAL)

def create_transition(transition_line, species_idxes):
	transition_info = transition_line.split("\t")
	sep_idx = transition_info.index(">")
	tname = transition_info[0]
	reactants = transition_info[1:sep_idx]
	products = transition_info[sep_idx + 1:len(transition_info) - 1]
	rate = float(transition_info[len(transition_info) - 1])
	transition_vector = np.array([0 for _ in range(len(species_idxes))])
	always_enabled = False
	is_consumer = False
	for reactant in reactants:
		if reactant == "0":
			always_enabled = True
			break
		transition_vector[species_idxes[reactant]] -= 1
	for product in products:
		if product == "0":
			is_consumer = True
			break
		transition_vector[species_idxes[product]] += 1
	if always_enabled:
		return Transition(transition_vector, lambda state : True, lambda state : rate, tname)
	reactant_idxes = [species_idxes[reactant] for reactant in reactants]
	return Transition(transition_vector, lambda state : np.all([state[i] > 0 for i in reactant_idxes]), lambda state : rate, tname)

def parse_ragtimer(filename):
	with open(filename, 'r') as rag:
		lines = rag.readlines()
		assert(len(lines) >= 4)
		# The first line is the list of reactants
		species_idxes = create_species_idxes(lines[0].split("\t"))
		# Get the initial condition
		init_state = np.array([int(val) for val in lines[1].split("\t")])
		# Get the boundary condition (exact equal)
		bound_vals = lines[2].split("\t")
		boundary = [create_bound(bound_val) for bound_val in bound_vals]
		transitions = [create_transition(line, species_idxes) for line in lines[3:]]
		return Crn(transitions, boundary, init_state)

def parse_dependency_ragtimer(filename):
	with open(filename, 'r') as rag:
		lines = rag.readlines()
		assert(len(lines) >= 4)
		return DepGraph(lines), parse_ragtimer(filename)
