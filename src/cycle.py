import numpy as np

class Cycle:
	def __init__(self, ordered_reactions : list):
		'''
	Creates an object 
		'''
		self.__ordered_reactions = ordered_reactions
		self.__check_cycle_valid()
		self.__in_s0 = [r.in_s0 for r in ordered_reactions]
		if not np.any(self.__in_s0):
			raise Exception("Cycle must leave S0 (else it may be useless)!")
		self.is_orthocycle = np.all(self.__in_s0)
	
	def __check_cycle_valid(self) -> bool:
		sum = self.__ordered_reactions[0].vec_as_mat
		for r in self.__ordered_reactions[1::]:
			sum += r.vec_as_mat
		assert(np.all(np.is_close(sum, 0)))

	def apply_cycle(self, state : np.vector) -> list:
		s = state
		new_states = []
		for r in self.__ordered_reactions:
			s += r.vec_as_mat
			new_states.append(s)
		return new_states

def get_commutable_transitions(crn : Crn, s0 : Subspace, ss : Subspace) -> list:
	transitions = []
	for t in crn.transitions:
		# TODO: need offset?
		if np.all(np.is_close(s0.P * t.vec_as_mat, 0)).all() and np.all(np.isclose(ss.P * t.vec_as_mat, t.vec_as_mat)):
			transitions.append(t)
	return t
