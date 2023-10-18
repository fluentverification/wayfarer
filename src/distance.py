import numpy as np

class BoundTypes:
	LESS_THAN=0
	LESS_THAN_EQ=0
	EQUAL=1
	GREATER_THAN=2
	GREATER_THAN_EQ=2

def species_distance(value, bound, bound_type=BoundTypes.EQUAL, normalize=True):
	norm_factor = 1 if normalize else abs(value - bound)
	if bound_type == BoundTypes.EQUAL:
		return abs(value - bound) / norm_factor
	elif bound_type == BoundTypes.LESS_THAN: # LESS_THAN_EQ is counted here too
		return min(value - bound, 0) / norm_factor
	elif bound_type == BoundTypes.GREATER_THAN:
		return min(bound - value, 0) / norm_factor
	else:
		raise Exception(f"bound_type not supported: '{bound_type}'!")

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

def vass_distance(state,
