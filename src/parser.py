from crn import *
from distance import Bound, BoundTypes
from dependency import *

# from stormpy import Rational

import sys
# Cursed
import importlib

# A dynamically loaded function pointer from a user file to customize how rates are found
custom_rate_finder = None

def parse_custom_rate_finder(filename : str):
	'''
Loads `filename` (a python file) and looks for a function called rate_finder which is then
used to generate rates for a CRN
	'''
	global custom_rate_finder # this should be the ONLY function which MODIFIES this
	assert(filename.endswith(".py"))
	module_name = filename.replace(".py", "").split("/").pop()
	# mod = importlib.import_module(module_name)
	module_spec = importlib.util.spec_from_file_location(module_name, filename)
	mod = importlib.util.module_from_spec(module_spec)
	module_spec.loader.exec_module(mod)
	try:
		custom_rate_finder = mod.rate_finder
	except Exception as e:
		print(f"[WARNING] Could not load custom rate finder. Reason: {e}")

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
	# Is 1.0 iff is reactant
	rate_mul_vector = np.array([0.0 for _ in transition_vector])
	catalysts = np.array([0.0 for _ in transition_vector])
	for reactant in reactants:
		if reactant == "0":
			always_enabled = True
			break
		idx = species_idxes[reactant]
		transition_vector[idx] -= 1
		rate_mul_vector[idx] = 1.0
	for product in products:
		if product == "0":
			is_consumer = True
			break
		elif product in reactants:
			idx = species_idxes[product]
			catalysts[idx] = 1.0
		transition_vector[species_idxes[product]] += 1
	# The rate, from rate constant k and reactants A, B, is k * A^count(A) * B^count(B)
	rate_finder = None
	if custom_rate_finder is None:
		rate_finder = lambda state : rate_const * np.prod([state[i] ** rate_mul_vector[i] for i in range(len(rate_mul_vector))])
	else:
		rate_finder = lambda state : custom_rate_finder(state, rate_const, tname)
	if always_enabled:
		return Transition(transition_vector
					, lambda state : True
					, rate_finder
					, tname
					, rate_const
					, catalysts)
	reactant_idxes = [species_idxes[reactant] for reactant in reactants]
	# Require all reactants to be strictly greater than zero
	return Transition(transition_vector
				, lambda state : np.all([state[i] > 0 for i in reactant_idxes])
				, rate_finder
				, tname
				, rate_const
				, catalysts)

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

def parse_dependency_ragtimer(filename: str, agnostic : bool =False):
	with open(filename, 'r') as rag:
		lines = rag.readlines()
		assert(len(lines) >= 4)
		return DepGraph(lines, agnostic=agnostic), parse_ragtimer(filename)

def create_piped(crn : Crn, use_rc : bool = False):
	'''
	Creates a piped matrix. Assumes that the Crn has constant rates since rates are gleaned
	from the rate finder at the initial state.
	'''
	matrix = None
	if use_rc:
		matrix = np.column_stack([
			# Normalize the vector
			t.vec_as_mat / np.linalg.norm(t.vec_as_mat)
			# Scale by the transition rate
			* t.rate_constant
			for t in crn.transitions])
	else:
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
