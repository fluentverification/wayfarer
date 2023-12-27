# TODO: re-implement Landon's dependency graph

import numpy as np

from crn import *
from subspace import *

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
		line_elems = line.replace("\n", "").split("\t")
		# print(line_elems)
		self.name = line_elems[0]
		# print(f"Creating reaction name {self.name}")
		sep_idx = line_elems.index(">")
		# ===============================
		# actual names of species
		self.in_species = line_elems[1:sep_idx]
		if "0" in self.in_species:
			self.in_species.remove("0")
		# print("in", self.in_species)
		self.out_species = line_elems[sep_idx + 1:len(line_elems) - 1]
		if "0" in self.out_species:
			self.out_species.remove("0")
		# print("out", self.out_species)
		# ===============================
		for spec in self.out_species:
			if not spec in self.in_species:
				dep_graph.declare_producer(self, spec)
			# else:
				# self.out_species.remove(spec)
				# self.in_species.remove(spec)
		for spec in self.in_species:
			if not spec in self.out_species:
				dep_graph.declare_consumer(self, spec)

	def __str__(self):
		return str(self.name)

	def vec(self, species_names):
		v = []
		for n in species_names:
			if n in self.in_species:
				v.append(-1.0)
			elif n in self.out_species:
				v.append(1.0)
			else:
				v.append(0.0)
		return np.matrix(v).T

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
		self.desired_values = desired_values
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
		self.species_names = [spec_name.strip() for spec_name in ragtimer_lines[0].split("\t")]
		self.reactions = [Reaction(line, self) for line in ragtimer_lines[3:]]
		# We want values of zero when we find a -1 and a value of 1 for all others
		self.mask = np.matrix([int(int(val.strip()) != -1) for val in ragtimer_lines[2].split("\t")]).T
		desired_values = np.matrix([int(val.strip()) for val in ragtimer_lines[2].split("\t")]).T
		self.desired_values = desired_values
		init_state = np.matrix([int(val) for val in ragtimer_lines[1].split("\t")]).T
		change = np.multiply(desired_values - init_state, self.mask)
		self.visited_changes = {}
		self.used_reactions = {}
		self.used_species = {}
		self.init_state = init_state
		# self.produced_species = {}
		self.graph_root = self.create_graph(change)
		self.create_reaction_levels()

	def create_graph(self, change, level = 0):
		'''
		Recursive function that creates the graph
		'''
		# Decreases the values in the change vector
		successors = []
		hashable_change = tuple([float(f) for f in change])
		# print(hashable_change)
		if hashable_change in self.visited_changes:
			return
		self.visited_changes[hashable_change] = True
		# Set each reaction as visited
		species_idx = 0
		allowed_reactions = []
		for c in change:
			species = self.species_names[species_idx]
			if c > 0:
				self.used_species[species] = True
				if not species in self.producers:
					print(f"Unable to continue this path to satisfiability! (This is not an error): produce {species}")
					continue
				if self.init_state[species_idx] > 0:
					continue
				spec_producers = self.producers[species]
				for producer in spec_producers:
					if not producer.name in self.used_reactions:
						allowed_reactions.append(producer.name)
					self.used_reactions[producer.name] = True
			elif c < 0:
				self.used_species[species] = True
				if not species in self.consumers:
					print(f"Unable to continue this path to satisfiability! (This is not an error): consume {species}")
					continue
				spec_consumers = self.consumers[species]
				for consumer in spec_consumers:
					if not consumer.name in self.used_reactions:
						allowed_reactions.append(consumer.name)
					self.used_reactions[consumer.name] = True
			species_idx += 1
		# Actually visit new ones
		species_idx = 0
		# If c is all zeros, we don't recurse any further
		for c in change:
			species = self.species_names[species_idx]
			# print(species)
			if c > 0:
				# TODO: create a new mask where the current index is zero and the (current) producer index is one
				# We need a producer reaction
				# print(self.producers)
				if not species in self.producers:
					continue
				if self.init_state[species_idx] > max(self.desired_values[species_idx], 0):
					continue
				spec_producers = self.producers[species]
				for producer in spec_producers:
					if producer.name in self.used_reactions and not producer.name in allowed_reactions:
						continue
					producer_idxs = [self.species_names.index(producer_species) for producer_species in producer.in_species if producer_species not in self.used_species]
					new_change = np.matrix([float(i in producer_idxs) for i in range(len(self.mask))]).T
					count = abs(c)
					# self.declare_reaction_at_level(producer, level)
					# recurse down until satisfied
					next_reactions = self.create_graph(new_change, level + 1)
					successors.append(Node(producer, next_reactions, level + 1, count))
			elif c < 0:
				# We need a consumer reaction
				# print(self.consumers)
				if not species in self.consumers:
					continue
				if self.init_state[species_idx] < max(self.desired_values[species_idx], 0):
					continue

				spec_consumers = self.consumers[species]

				for consumer in spec_consumers:
					if consumer.name in self.used_reactions and not consumer.name in allowed_reactions:
						continue
					consumer_idxs = [self.species_names.index(consumer_species) for consumer_species in consumer.in_species if consumer_species not in self.used_species]
					new_change = np.matrix([0.0 - float(i in consumer_idxs) for i in range(len(self.mask))]).T
					count = abs(c)
					# self.declare_reaction_at_level(consumer, level)
					next_reactions = self.create_graph(new_change, level + 1)
					# recurse down until satisfied
					successors.append(Node(consumer, next_reactions, level + 1, count))
			species_idx += 1

		return successors

	def create_reaction_levels(self):
		for node in self.graph_root:
			self.create_reaction_level(node)
		self.reaction_levels.reverse()

	def create_reaction_level(self, node : Node):
		if node is None:
			return -1
		if node.children is None or len(node.children) == 0:
			self.declare_reaction_at_level(node.reaction, 0)
			return 0
		successor_levels = [self.create_reaction_level(n) for n in node.children]
		level = 1 + min(successor_levels)
		self.declare_reaction_at_level(node.reaction, level)
		return level

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
		# print(f"Declaring {reaction} as producer of {species}")
		if not species in self.producers:
			self.producers[species] = [reaction]
		else:
			self.producers[species].append(reaction)

	def declare_consumer(self, reaction, species):
		'''
		Declares a reaction as a consumer (called by the reaction's constructor)
		'''
		# print(f"Declaring {reaction} as consumer of {species}")
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
			s = f"{s}Level {lev_idx}:"
			for reaction in level:
				s = f"{s} {reaction}"
			s = f"{s}\n"
		return s

	def create_subspaces(self, crn):
		available_reactions = []
		indecies = []
		for lev_idx in range(len(self.reaction_levels)):
			level = self.reaction_levels[lev_idx]
			for reaction in level:
				transition = crn.find_transition_by_name(reaction.name)
				available_reactions.append(transition)
			indecies.append(len(available_reactions) - 1)
		# The last index is unnecessary as it is just the last available index
		# indecies.pop()
		subspaces = []
		# print(indecies)
		# Perhaps restrict the reactions only to the current level?
		# print(indecies, "and ", [t.name for t in available_reactions])
		for i in range(len(indecies)):
			idx = indecies[i]
			last_idx = indecies[i - 1] if i > 1 else 0
			# TODO: should be idx or idx + 1
			used_transitions = available_reactions[:idx + 1]
			unused_transitions = available_reactions[idx + 1:]
			last_layer = available_reactions[last_idx:idx]
			subspace = Subspace(used_transitions, unused_transitions) #, last_layer)
			# print(f"Last layer is {[t.name for t in last_layer]}")
			subspaces.append(subspace)
		subspaces.reverse()
		# print([str(subspace) for subspace in subspaces])
		return subspaces
