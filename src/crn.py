import numpy as np

class BoundTypes:
	LESS_THAN=0
	LESS_THAN_EQ=1
	EQUAL=2
	GREATER_THAN=3
	GREATER_THAN_EQ=4
	DONT_CARE=5

class Bound:
	def __init__(self, bound, bound_type):
		self.bound = bound
		self.bound_type = bound_type

class Transition:
	# vector      #: np.vector
	# enabled     #: lambda
	# rate_finder #: lambda
	def __init__(self, vector, enabled, rate_finder):
		self.vector = np.array(vector)
		self.vec_as_mat = np.matrix(vector).T
		self.enabled_lambda = enabled
		self.rate_finder = rate_finder
		self.can_eliminate = False

	def enabled(self, state):
		return self.enabled_lambda(state)

class Crn:
	def __init__(self, transitions=None, boundary=None, init_state=None, all_trans_always_enabled=False, all_rates_const=False):
		self.boundary = boundary
		if transitions is None:
			self.transitions = []
		else:
			self.transitions = transitions
		self.init_state = init_state
		# self.prune_transitions()
		self.react_depriority = np.array([0.0 if self.boundary[i].bound_type == BoundTypes.DONT_CARE else 1.0 for i in range(len(boundary))])
		self.all_trans_always_enabled = all_trans_always_enabled
		self.all_rates_const = all_rates_const

	def add_transition(self, transition : Transition):
		self.transitions.append(transition)

	def prune_transitions(self):
		'''
Unused: prunes transitions
		'''
		for i in range(len(self.boundary)):
			bound = self.boundary[i]
			if bound.bound_type == BoundTypes.DONT_CARE:
				self.__prune_transitions_at_index(i)

	def __prune_transitions_at_index(self, index):
		'''
Unused
		'''
		for transition in self.transitions:
			can_eliminate = True
			for i in range(len(transition.vector)):
				if i != index and transition.vector[i] != 0.0:
					can_eliminate = False
			if can_eliminate:
				transition.can_eliminate = True
			# TODO: prune it

