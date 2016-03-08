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
		self.genome = r.randint(2,size=(self.nb_gene,self.nb_gene))
		# print genome

		# the differente probability for the mutations
		self.ponctual = p_ponctual_mutation
		self.inversion = p_inversion_mutation
		self.Xover = p_Xover_mutation

		# to count the nunmber of differente mutation
		self.nb_ponctual = 0
		self.nb_inversion = 0
		self.nb_Xover = 0

		# the probability to mutate
		self.mutation = p_mutation



	####################################################
	# function to do ponctual mutations
	####################################################
	def PonctualMutation(self, gene):
		a = r.randint(self.nb_gene)
		gene[a] = abs(gene[a]-1)

		return gene


	####################################################
	# function to do inversion mutations
	####################################################
	def InversionMutation(self, gene):
		b, c  = r.randint(self.nb_gene), r.randint(self.nb_gene)
		while b == c:
			c = r.randint(self.nb_gene)

		gene[min(b,c):max(b,c)] = abs(gene[min(b,c):max(b,c)]-1)
		return gene



	####################################################
	# function to do Xover mutations
	####################################################
	def XoverMutation(self, gene, gene2):
		b, c = r.randint(self.nb_gene), r.randint(self.nb_gene)
		while b == c:
			c = r.randint(self.nb_gene)

		d = r.randint(self.nb_gene)
		e = d + r.choice([-1,1])*(max(b,c)-min(b,c))
		while e < 0 or e > self.nb_gene:
			d = r.randint(self.nb_gene)
			e = d + r.choice([-1,1])*(max(b,c)-min(b,c))

		gene_int = gene[min(b,c):max(b,c)]
		print b, c, d, e
		print len(gene_int), len(gene2[min(d,e):max(d,e)])
		gene[min(b,c):max(b,c)] = gene2[min(d,e):max(d,e)]
		gene2[min(d,e):max(d,e)] = gene_int
		return gene, gene2


	####################################################
	# function to mutate the genome
	####################################################
	def MutateGenome(self, fitness):

		# define the formula to associate the fitness to the number of mutation
		nb_mutation = int(fitness * 5)

		numero_mutation = []
		gene_mod = []

		for i in xrange(nb_mutation):
			while True:

				# choose one gene among the others genes (will it be one mutation on one gene or can we allow several mutations
				# on one genes ?)
				n_gene = r.randint(self.nb_gene)
				gene = self.genome[n_gene,:]

				# if the mutation is ponctual
				if r.random() < self.ponctual:
					gene_mod.append(self.PonctualMutation(gene))
					self.nb_ponctual += 1
					numero_mutation.append(n_gene)
					break
					
				# if the mutation is an inversion
				if r.random() < self.inversion:
					gene_mod.append(self.InversionMutation(gene))
					self.nb_inversion += 1
					numero_mutation.append(n_gene)
					break

				# if the mutation is an Xover
				if r.random() < self.Xover:
					n_gene2 = n_gene
					while n_gene2 == n_gene:
						n_gene2 = r.randint(self.nb_gene)
					gene2 = self.genome[n_gene2,:]

					gene, gene2 = self.XoverMutation(gene, gene2)
					self.nb_Xover += 1
					numero_mutation.append([n_gene, n_gene2])
					gene_mod.append([gene, gene2])
					break

		return numero_mutation, gene_mod


# genome(number_gene, p_ponctual_mutation, p_inversion_mutation, p_Xover_mutation, p_mutation)

a = genome(10, 0, 0, 1, 0)

print a.genome

b, c = a.MutateGenome(2)

print b
print c