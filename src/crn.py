class Transition:
	vector      : np.vector
	enabled     : lambda
	rate_finder : lambda
	def __init__(self, vector : np.vector, enabled : lambda, rate_finder : lambda):
		self.vector = vector
		self.enabled = enabled
		self.rate_finder = rate_finder

class Crn:
	def __init__(self, transitions=None : list):
		if transitions is None:
			self.transitions = []
		else:
			self.transitions = transitions

	def add_transition(self, transition : Transition):
		self.transitions.append(transition)
