# TODO: re-implement Landon's dependency graph

class Reaction:
	def __init__(self, name, in_species, out_species):
		self.name = name
		self.in_species = in_species
		self.out_species = out_species

	def __init__(self, line):
		line_elems = line.replace("\n", "").split("\n")
		self.name = line_elems
		sep_idx = transition_info.index(">")
		self.in_species = transition_info[1:sep_idx]
		self.out_species = transition_info[sep_idx + 1:len(transition_info) - 1]

class DepGraph:
	def __init__(self, reactions, desired
