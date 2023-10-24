import numpy as np

class Transition:
	# vector      #: np.vector
	# enabled     #: lambda
	# rate_finder #: lambda
	def __init__(self, vector, enabled, rate_finder):
		self.vector = np.matrix(vector)
		self.enabled_lambda = enabled
		self.rate_finder = rate_finder

	def enabled(self, state):
		return self.enabled_lambda(state)

class Crn:
	def __init__(self, transitions=None, boundary=None, init_state=None):
		self.boundary = boundary
		if transitions is None:
			self.transitions = []
		else:
			self.transitions = transitions
		self.init_state = init_state

	def add_transition(self, transition : Transition):
		self.transitions.append(transition)

