import numpy as np

from crn import *

# from counterexample import get_transitions

def get_transitions(state, crn : Crn):
	'''
Use properties of a VASS to get the enabled transitions of a particular
	'''
	# global all_transitions
	transitions = []
	for transition in crn.transitions:
		if transition.enabled(state):
			next_state = state + transition.vector
			transitions.append((transition.rate_finder(state), next_state))
	return transitions

def species_distance(value, bound, bound_type=BoundTypes.EQUAL, normalize=True):
	norm_factor = 1 if normalize else abs(value - bound)
	if bound_type == BoundTypes.DONT_CARE:
		return 0.0
	elif bound_type == BoundTypes.EQUAL:
		return abs(value - bound) / norm_factor
	elif bound_type == BoundTypes.LESS_THAN or bound_type == BoundTypes.LESS_THAN_EQ: # LESS_THAN_EQ is counted here too
		return max(value - bound, 0) / norm_factor
	elif bound_type == BoundTypes.GREATER_THAN or bound_type == BoundTypes.GREATER_THAN:
		return max(bound - value, 0) / norm_factor
	else:
		raise Exception(f"bound_type not supported: '{bound_type}'!")

def species_satisfies(value, bound, bound_type=BoundTypes.EQUAL):
	if bound_type == BoundTypes.DONT_CARE:
		return True
	elif bound_type == BoundTypes.EQUAL:
		return value == bound
	elif bound_type == BoundTypes.LESS_THAN:
		return value < bound
	elif bound_type == BoundTypes.LESS_THAN_EQ:
		return value <= bound
	elif bound_type == BoundTypes.GREATER_THAN:
		return value > bound
	elif bound_type == BoundTypes.GREATER_THAN:
		return value >= bound
	else:
		raise Exception(f"bound_type not supported: '{bound_type}'!")

def satisfies(state, boundary):
	'''
Checks to see if a state satisfies (or is within) the boundary
	'''
	assert(len(state) == len(boundary))
	for i in range(len(state)):
		# dist = species_distance(state[i], boundary[i].bound, boundary[i].bound_type)
		# print(f"Distance = {dist}")
		if not species_satisfies(state[i], boundary[i].bound, boundary[i].bound_type):# dist != 0:
			# print(f"failed on index {i}: val {state[i]}, bound type: {boundary[i].bound_type}")
			return False
	return True

def averaging_total_distance(state, bounds, bound_types, weights=None, normalize=True):
	assert(len(state) == len(bounds) == len(bound_types))
	if weights is None:
		weights = [1.0 for _ in range(len(state))]
	else:
		assert(len(weights) == len(state))

	return np.average([
	species_distance(
		state[i]
		, bounds[i]
		, bound_types[i]
		, normalize)
	* weights[i] for i in range(len(state))])

def direction_vector(state, transitions, dist_weight=0.5, rate_weight=0.5, normalize=True):
	'''
Parameters:
state: the current state we're in
transitions: The transitions from this state (enabled), including the rate and target_state
dist_weight: The amount of weight to put in the distance from the current state to the next state in the direction vector
rate_weight: The amount of weight to put in the CTMC rate
	'''
	vecs = np.zeros(len(state))
	normFactor = 1.0
	if normalize:
		norm_factor = 0.0
		for rate, target_state in transitions:
			norm_factor += rate
	# Now actually create the vector
	for rate, target_state in transitions:
		direction_vec = target_state - state
		distance = np.linalg.norm(direction_vec)
		normal_vec = direction_vec / distance
		vecs += normal_vec * (dist_weight * distance + rate_weight / norm_factor * rate)
	return vecs

def angle(v1, v2):
	return np.rad2deg(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))

def vass_distance(state, boundary, exact_equal=False):
	'''
Parameters:
state: the current state we're in
boundary: The variable boundaries
	'''
	# The shortest distance to the boundary is a superposition of the individual species'
	# boundary distance vectors, and since each species is an element in the vector, we just
	# need a vector with each species' distance in its index
	assert(len(state) == len(boundary))
	return np.array([
			species_distance(
				state[i]
				, boundary[i].bound
				, boundary[i].bound_type
			) for i in range(len(state))
		])

def vass_priority(state, boundary, crn, reach=1.0, include_flow_angle=False, include_flow_mag=False, deprioritize_dont_cares=False):
	'''
Creates a priority from a state based on a boundary. The lower the priority, the sooner we should explore the state.
This means you should use a min queue
	'''
	dist_vector = vass_distance(state, boundary)
	v_dist = None
	if deprioritize_dont_cares:
		v_dist = np.linalg.norm(np.multiply(dist_vector, crn.react_depriority))
	else:
		v_dist = np.linalg.norm(dist_vector)
	priority = v_dist

	if include_flow_angle:
		transitions = get_transitions(state, crn)
		flow_vector = direction_vector(state, transitions)
		flow_dist_angle = np.abs(angle(flow_vector, dist_vector))
		# print(f"Angle = {flow_dist_angle}")
		# TODO: figure out better weight
		WEIGHT = 1
		# We want to minimize flow_dist_angle
		priority += WEIGHT * flow_dist_angle / 180
		if include_flow_mag:
			flow_magnitude = np.linalg.norm(flow_vector)
			# Maximize the magnitude of the flow vector
			WEIGHT = 1
			priority += 1 / flow_magnitude

	# TODO: we want to minimize flow_dist_angle, and also minimize v_dist
	# if possible, also maximizing the magnitude of the flow vector, assuming
	# the angle is low enough
	# return 1 / ((flow_dist_angle) / flow_magnitude * reach + v_dist)
	return priority
