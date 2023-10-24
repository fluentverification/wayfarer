import numpy as np

class Transition:
	vector      #: np.vector
	enabled     #: lambda
	rate_finder #: lambda
	def __init__(self, vector, enabled, rate_finder):
		self.vector = np.vector(vector)
		self.enabled = enabled
		self.rate_finder = rate_finder

class Crn:
	def __init__(self, transitions=None : list, boundary : list, init_state):
		self.boundary = boundary
		if transitions is None:
			self.transitions = []
		else:
			self.transitions = transitions
		self.init_state = init_state

	def add_transition(self, transition : Transition):
		self.transitions.append(transition)

