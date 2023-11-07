import numpy as np

from distance import vass_distance

class RVecLevel:
	def __init__(self, vector, level):
		self.vector = vector
		self.level  = level

# TODO: multiple paths to the target?
# Maybe store projection matrices? P = A(A^T A)^-1 A

class SubspacePriority:
	# These two must be set prior to using
	# ordered_reaction_vectors must be in order and must be at RVecLevel
	ordered_reaction_vectors = []
	boundary = None
	def __init__(self, state):
		self.__distance = vass_distance(state, SubspacePriority.boundary)
		# TODO: Optimize if we are already in one space
		self.__subspace_distances = [
			self.__dist_to_subspace(
				state
				# TODO
				, SubspacePriority.ordered_reaction_vectors[i:]
			) for i in range(len(SubspacePriority.ordered_reaction_vectors) - 1)
		]

	def __dist_to_subspace(self, state, B):
		P = A * (A.T * A) ** -1 * A
		dist_vect = P * state
		return np.linalg.norm(dist_vect)

	# Optimized NEQ with early termination
	def __ne__(self, other):
		if self.__distance != other.__distance:
			return True
		for i in range(len(self.__subspace_distances)):
			if self.__subspace_distances != other.__subspace_distances:
				return True
		return False

	def __eq__(self, other):
		return not self != other

	# We must implement at least one element from the groups
	# of [<, >=] and [>, <=] to accurately allow short-circuiting
	def __gt__(self, other):
		for i in range(len(self.__subspace_distances)):
			if self.__subspace_distances > other.__subspace_distances:
				return True
		if self.__distance > other.__distance:
			return True
		return False

	def __ge__(self, other):
		for i in range(len(self.__subspace_distances)):
			if self.__subspace_distances >= other.__subspace_distances:
				return True
		if self.__distance >= other.__distance:
			return True
		return False

	def __lt__(self, other)
		return not self >= other

	def __le__(self, other):
		return not self > other


