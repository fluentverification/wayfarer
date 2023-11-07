import numpy as np

from distance import vass_distance

class SubspacePriority:
	# These two must be set prior to using
	ordered_reaction_vectors = []
	boundary = None
	def __init__(self, state):
		self.__distance = vass_distance(state, SubspacePriority.boundary)
		# TODO: Optimize if we are already in one space
		self.__subspace_distances = [
			self.__dist_to_subspace(
				state
				, SubspacePriority.ordered_reaction_vectors[i:]
			) for i in range(len(SubspacePriority.ordered_reaction_vectors) - 1)
		]

	def __dist_to_subspace(self, state, basis_vectors):
		# TODO
		pass

	def __eq__(self, other):
		return self.__distance == other.__distance and np.all(self.__subspace_distances == other.__subspace_distances)

	# Optimized
	def __ne__(self, other):
		if self.__distance != other.__distance:
			return True
		for i in range(len(self.__subspace_distances)):
			if self.__subspace_distances != other.__subspace_distances:
				return True
		return False

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
