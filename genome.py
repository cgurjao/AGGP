import networkx as nx
import numpy as np
import matplotlib as np
import numpy.random as r



class genome:
	"""docstring for genome"""
	def __init__(self, number_gene, p_ponctual_mutation, p_inversion_mutation, p_Xover_mutation, p_mutation):
		
		# The number of the genes we want
		self.nb_gene = number_gene

		# to create the random matrix for the initialisation, full of 0 and 1
		self.genome = r.randint(2,size=(nb_gene,nb_gene))
		# print genome

		# the differente probability for the mutations
		self.ponctual = p_ponctual_mutation
		self.inversion = p_inversion_mutation
		self.Xover = p_Xover_mutation

		# the probability to mutate
		self.mutation = p_mutation


	####################################################
	# function to mutate the genome
	####################################################
	def MutateGenome(self, fitness):

		# define the formula to associate the fitness to the number of mutation
		nb_mutation = int(fitness * 5)

		for i in xrange(nb_mutation):

			# choose one gene among the others genes (will it be one mutation on one gene or can we allow several mutations
			# on one genes ?)
			n_gene = r.randint(self.nb_gene)
			gene = self.genome[n_gene,:]

			