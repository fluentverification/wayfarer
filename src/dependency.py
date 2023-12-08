# TODO: re-implement Landon's dependency graph

class Node:
	'''
	A Node in the dependency graph
	'''
	def __init__(self, reaction, level = 0, count = 1):
		'''
		A constructor for a node with no successors
		'''
		self.reaction = reaction
		self.children = []
		self.level = level
		self.count = count

	def __init__(self, reaction, children, level = 0, count = 1):
		'''
		A constructor for a node with successors (children)
		'''
		self.reaction = reaction
		self.children = children
		self.level = level
		self.count = count

	def add_child(self, child):
		'''
		Appends a successor to the node's children
		'''
		self.children.append(child)

class Reaction:
	'''
	A reaction (can read in a line from a ragtimer file)
	'''
	def __init__(self, name, in_species, out_species):
		'''
		Basic constructor if we know species used and produced
		'''
		self.name = name
		self.in_species = in_species
		self.out_species = out_species

	def __init__(self, line, dep_graph):
		'''
		A constructor which parses a RAGTIMER line
		'''
		line_elems = line.replace("\n", "").split("\n")
		self.name = line_elems
		sep_idx = line.index(">")
		# ===============================
		# actual names of species
		self.in_species = line[1:sep_idx]
		self.out_species = line[sep_idx + 1:len(line) - 1]
		# ===============================
		for spec in self.out_species:
			if not spec in self.in_species:
				dep_graph.declare_producer(self, spec)
		for spec in self.in_species:
			if not spec in self.out_species:
				dep_graph.declare_consumer(self, spec)

	def __str__(self):
		return self.name

class DepGraph:
	'''
	A basic graph class for the entire dependency graph
	'''
	def __init__(self, reactions, species_names, desired_values, init_state):
		'''
		Builds the dependency graph
		'''
		self.producers = {}
		self.consumers = {}
		self.reactions = []
		self.reaction_levels = []
		self.species_names = species_names
		self.root = []
		self.mask = np.matrix([1.0 for _ in range(len(init_state))]).T
		# Creates the reactions
		for reaction in reactions:
			r = Reaction(reaction, self)
			self.reactions.append(r)
		self.graph_root = self.create_graph(desired_values, init_state)

	def __init__(self, ragtimer_lines):
		self.producers = {}
		self.consumers = {}
		self.reaction_levels = []
		self.root = []
		self.species_names = [spec_name.split() for spec_name in lines[0].split("\t")]
		self.reactions = [Reaction(line, self) for line in ragtimer_lines[3:]]
		# We want values of zero when we find a -1 and a value of 1 for all others
		self.mask = np.matrix([int(int(val.strip()) != -1) for val in ragtimer_lines[2].split("\t")]).T
		desired_values = np.matrix([int(val.strip) for val in ragtimer_lines[2].split("\t")]).T
		init_state = np.matrix([int(val) for val in ragtimer_lines[1].split("\t")]).T
		self.graph_root = self.create_graph(desired_values, init_state)

	def create_graph(self, desired_values, state, level = 0):
		'''
		Recursive function that creates the graph
		'''
		# Decreases the values in the change vector
		species_idx = 0
		successors = []
		change = np.multiply(desired_values - state, self.mask)
		# If c is all zeros, we don't recurse any further
		for c in change:
			species = species_names[species_idx]
			if c > 0:
				# We need a producer reaction
				spec_producers = self.producers[species]
				for producer in spec_producers:
					new_state = state.copy()
					new_state[species_idx] = 0
					count = abs(c)
					self.declare_reaction_at_level(producer, level)
					# recurse down until satisfied
					next_reactions = self.create_graph(desired_values, new_state)
					successors.append(Node(producer, next_reactions, level + 1, count))
			elif c < 0:
				# We need a consumer reaction
				spec_consumers = self.consumers[species]
				for consumer in spec_consumers:
					new_state = state.copy()
					new_state[species_idx] = 0
					count = abs(c)
					self.declare_reaction_at_level(consumer, level)
					next_reactions = self.create_graph(desired_values, new_state)
					# recurse down until satisfied
					successors.append(Node(consumer, next_reactions, level + 1, count))
			species_idx += 1
		return successors

	def declare_reaction_at_level(self, reaction, level):
		'''
		Declares the reaction at a specific level in the graph (for easy parsing later without having to
		traverse the actual graph)
		'''
		while len(self.reaction_levels) <= level:
			self.reaction_levels.append([])
		self.reaction_levels[level].append(reaction)

	def declare_producer(self, reaction, species):
		'''
		Declares a reaction as a producer (called by the reaction's constructor)
		'''
		if not species in self.producers:
			self.producers[species] = [reaction]
		else:
			self.producers[species].append(reaction)

	def declare_consumer(self, reaction, species):
		'''
		Declares a reaction as a consumer (called by the reaction's constructor)
		'''
		if not species in self.consumers:
			self.consumers[species] = [reaction]
		else:
			self.consumers[species].append(reaction)

	def __str__(self):
		'''
		Produces a string representation of the graph
		'''
		s = f"Dependency Graph:\nTotal Reactions: {len(self.reactions)}\n"
		s = f"{s}Producing Reactions: {len(self.producers)}\nConsuming Reactions: {len(self.consumers)}\n"
		for lev_idx in range(len(self.reaction_levels)):
			level = self.reaction_levels[lev_idx]
			s = f"{s}Level {level_idx}:"
			for reaction in level:
				s = f"{s} {reaction}"
			s = f"{s}\n"
		return s
