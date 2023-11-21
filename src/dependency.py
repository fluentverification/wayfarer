# TODO: re-implement Landon's dependency graph

class Node:
	def __init__(self, reaction):
		self.reaction = reaction
		self.children = []

	def __init__(self, reaction, children):
		self.reaction = reaction
		self.children = children

	def add_child(self, child):
		self.children.append(child)

class Reaction:
	def __init__(self, name, in_species, out_species):
		self.name = name
		self.in_species = in_species
		self.out_species = out_species

	def __init__(self, line, dep_graph):
		line_elems = line.replace("\n", "").split("\n")
		self.name = line_elems
		sep_idx = transition_info.index(">")
		self.in_species = transition_info[1:sep_idx]
		self.out_species = transition_info[sep_idx + 1:len(transition_info) - 1]
		for spec in self.out_species:
			if not spec in self.in_species:
				dep_graph.declare_producer(self, spec)
		for spec in self.in_species:
			if not spec in self.out_species:
				dep_graph.declare_consumer(self, spec)

class DepGraph:
	def __init__(self, reactions, species_names, desired_values, init_state):
		change = desired_values - init_state
		self.producers = {}
		self.consumers = {}
		self.reactions = []
		self.root = []
		for reaction in reactions:
			r = Reaction(reaction, self)
			self.reactions.append(r)
		for c in change:
			if c > 0:
				# TODO: recurse down until satisfied
				pass
			elif c < 0:
				# TODO: also recurse down
				pass

	def declare_producer(self, reaction, species):
		if not species in self.producers:
			self.producers[species] = [reaction]
		else:
			self.producers[species].append(reaction)

	def declare_consumer(self, reaction, species):
		if not species in self.consumers:
			self.consumers[species] = [reaction]
		else:
			self.consumers[species].append(reaction)
