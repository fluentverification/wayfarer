# TODO: re-implement Landon's dependency graph

class Node:
	def __init__(self, reaction, level = 0):
		self.reaction = reaction
		self.children = []
		self.level = level

	def __init__(self, reaction, children, level = 0):
		self.reaction = reaction
		self.children = children
		self.level = level

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
		sep_idx = line.index(">")
		self.in_species = line[1:sep_idx]
		self.out_species = line[sep_idx + 1:len(line) - 1]
		for spec in self.out_species:
			if not spec in self.in_species:
				dep_graph.declare_producer(self, spec)
		for spec in self.in_species:
			if not spec in self.out_species:
				dep_graph.declare_consumer(self, spec)

class DepGraph:
	def __init__(self, reactions, species_names, desired_values, init_state):
		self.producers = {}
		self.consumers = {}
		self.reactions = []
		self.reaction_levels = []
		self.species_names = species_names
		self.root = []
		for reaction in reactions:
			r = Reaction(reaction, self)
			self.reactions.append(r)
		species_idx = 0
		self.graph_root = self.create_graph(desired_values, init_state)

	def create_graph(self, desired_values, state, level = 0):
		successors = []
		change = desired_values - state
		for c in change:
			species = species_names[species_idx]
			if c > 0:
				# We need a producer reaction
				spec_producers = self.producers[species]
				for producer in spec_producers:
					new_state = state.copy()
					new_state[species_idx] += 1
					self.declare_reaction_at_level(producer, level)
					# recurse down until satisfied
					next_reactions = self.create_graph(desired_values, new_state)
					successors.append(Node(producer, next_reactions, level + 1))
			elif c < 0:
				# We need a consumer reaction
				spec_consumers = self.consumers[species]
				for consumer in spec_consumers:
					new_state = state.copy()
					new_state[species_idx] -= 1
					self.declare_reaction_at_level(consumer, level)
					next_reactions = self.create_graph(desired_values, new_state)
					# recurse down until satisfied
					successors.append(Node(consumer, next_reactions, level + 1))
			species_idx += 1
		return successors

	def declare_reaction_at_level(self, reaction, level);
		while len(self.reaction_levels) <= level:
			self.reaction_levels.append([])
		self.reaction_levels[level].append(reaction)

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
