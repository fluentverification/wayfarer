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
	rate_const = float(transition_info[len(transition_info) - 1])
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
	rate_mul_vector = np.array([float(elem < 0) for elem in transition_vector])
	# The rate, from rate constant k and reactants A, B, is k * A^count(A) * B^count(B)
	rate_finder = lambda state : rate_const * np.prod([state[i] ** rate_mul_vector[i] for i in range(len(rate_mul_vector))])
	if always_enabled:
		return Transition(transition_vector
					, lambda state : True
					, rate_finder
					, tname)
	reactant_idxes = [species_idxes[reactant] for reactant in reactants]
	# Require all reactants to be strictly greater than zero
	return Transition(transition_vector
				   , lambda state : np.all([state[i] > 0 for i in reactant_idxes])
				   , rate_finder
				   , tname)

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

def parse_dependency_ragtimer(filename, agnostic=False):
	with open(filename, 'r') as rag:
		lines = rag.readlines()
		assert(len(lines) >= 4)
		return DepGraph(lines, agnostic=agnostic), parse_ragtimer(filename)

def create_piped(crn : Crn):
	'''
	Creates a piped matrix. Assumes that the Crn has constant rates since rates are gleaned
	from the rate finder at the initial state.
	'''
	matrix = np.column_stack([
		# Normalize the vector
		t.vec_as_mat / np.linalg.norm(t.vec_as_mat)
		# Scale by the transition rate
		* t.rate_finder(crn.init_state)
		for t in crn.transitions])
	# Warn when rank is not equal to the number of transitions and species.
	# In this case, the pseudoinverse must be used.
	rank = np.linalg.matrix_rank(matrix)
	if not (rank == len(crn.transitions) and rank == len(crn.transitions[0].vector)):
		print(f"Warning: must use pseudoinverse for piped!")
	return matrix
