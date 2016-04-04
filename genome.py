import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as r
import fitness
import copy
import random
import bisect
import collections
from paramsGlob import *
from operator import add
import csv


#=============================================================================
#					GENOME CLASS
#=============================================================================


class genome:
	def __init__(self, number_genes):
		
		# The number of the genes we want
		self.nb_genes = number_genes

		# to create the random matrix for the initialisation, full of 0 and 1
		self.genome = r.randint(2,size=(self.nb_genes,self.nb_genes))

		# to symmetrize the matrix
		for i in xrange(self.nb_genes):
			for j in xrange(i):
				self.genome[i,j] = 0
			self.genome[i,i] = 1

		self.nb_fig = 1

	#___________________________________________________
	# function to update the matrix

	def UpdateMatrix(self):
		for i in xrange(self.nb_genes):
			for j in xrange(i):
				if r.random() <  mutation_rate:
					self.genome[j,i] = abs(self.genome[j,i]-1)


	#___________________________________________________
	# function to do crossing over
	
	def CrossingOver(self, M1, M2):        

		n = self.nb_genes
		# to do the crossing over with a line
		line = r.randint(n)

		#lineint = np.array(M1[line,:])
		#M1[line,:] = np.array(M2[line, :])
		#M2[line,:] = lineint

		lineint = np.copy(M1[line,:])
		M1[line,:] = M2[line,:]
		M2[line,:] = lineint
		return M1, M2



	#___________________________________________________	
	# functions to draw the graph
	def graph(self):
		# to construct the graph
		G = nx.Graph()
		G.add_nodes_from(range(self.nb_genes))

		for i in xrange(self.nb_genes):
			for j in xrange(i+1):
				if self.genome[j,i] == 1:
					G.add_edge(i,j)
		return G


	def draw(self, G, indiv = 0, iteration=0, save=True):
		# if you use the option save = True, you will save also all the past graph but not the future graph
		if indiv==0:
			pos = nx.spring_layout(G)
			nx.draw(G, pos)
			node_labels = nx.get_node_attributes(G,'state')
			nx.draw_networkx_labels(G, pos, labels = node_labels)

			name = "Graph_Indiv_%d _Iteration_%d .png" %(indiv, iteration)
			plt.savefig(name)
			plt.close()


'''
		position = nx.circular_layout(G)

		fig = plt.figure(self.nb_fig)
		self.nb_fig += 1

		list_deg = np.array(nx.degree(G).values())

		nx.draw_networkx_nodes(G, pos=position, node_size=300, alpha=0.5, node_color=list_deg, cmap=plt.cm.YlOrRd)# RdPu, BuPu
		# http://matplotlib.org/examples/color/colormaps_reference.html
		nx.draw_networkx_edges(G, pos=position, alpha=0.2)
		nx.draw_networkx_labels(G,pos=position)
		plt.draw()
		if iteration != 9999 and save == True:

			name = "Graph_Indiv_%d _Iteration_%d .png" %(indiv, iteration)
			plt.savefig(name)
			plt.close()

		if save == True and iteration == 9999:
			plt.title('Graph')
			plt.savefig('graph')
			plt.close()'''


#=============================================================================
#								POPULATION CLASS
#=============================================================================


class population:
	def __init__(self, number_genomes, number_genes):
		
		# The number of individuals/genomes
		self.nb_genomes = number_genomes

		# List of genomes
		self.pop = [genome(number_genes) for i in xrange(number_genomes)]

		# Generation number...
		self.gen = 0



	#__________________________________________________________________________
	# Updating the population (mutations, crossing over, reproduction)

	def new_generation(self, compt):
		F = []		
		temp_avg = []
		temp_RSS = []	
		temp_hier = []	
		res = []
		avg = 0
		RSS = 0
		hier = 0								# list of fitnesses
		for i,g in enumerate(self.pop):					# For each individual
			g.UpdateMatrix()							# Do random mutation
			rd = r.randint(self.nb_genomes) 					

			#if rd!=i:
				#g.genome, self.pop[rd].genome = g.CrossingOver(g.genome, self.pop[rd].genome)  # Do random crossing over
			res = fitness.fitness_score(g.graph())
			temp_avg.append(res[0])
			temp_RSS.append(res[1])	# Calculate fitness for each ind
			temp_hier.append(res[2])

			if rd!=i:
				g.genome, self.pop[rd].genome = g.CrossingOver(g.genome, self.pop[rd].genome)   # Do random crossing over
			print "Iteration: %d, Individu: %d" %(compt, i)
			if i == 0:
				#fitness.draw_figure_scalefree(g.graph(), i, compt)
				#fitness.draw_figure_hierarchical(g.graph(), i, compt)
				avg = res[0]
				RSS = res[1] 
				hier = res[2]
				#g.draw(g.graph(), compt, i)

		ratio = 0
		for i in xrange(len(temp_avg)):
			ratio = ratio + (temp_avg[i]/temp_RSS[i])
		ratio = ratio/len(temp_avg)
		for i in xrange(len(temp_avg)):
			temp_avg[i] = temp_avg[i]/ratio
		for i in xrange(len(temp_avg)):
			ratio = ratio + (temp_hier[i]/temp_RSS[i])
		ratio = ratio/len(temp_avg)
		for i in xrange(len(temp_avg)):
			temp_hier[i] = temp_hier[i]/ratio

		for i in xrange(len(temp_avg)):
			F = map(add,temp_avg, temp_RSS)
			F = map(add,F, temp_hier)

		proba_rep = self.selection(F,0.7)
		#bprint sum(F)
		new_pop = []
		for j in xrange(self.nb_genomes):
			p = r.random()
			i = len(proba_rep)-1
			while p-proba_rep[i]>0 and i>0:
				p -= proba_rep[i]
				i -= 1
			new_pop.append(copy.deepcopy(self.pop[i]))
		self.pop = new_pop
		self.gen += 1
		return sum(F), avg, RSS, hier
			
	# Generates array of reproduction probabilities (one per genome)
	def selection(self,fit_table,c):
		order = np.argsort(fit_table)
		ranks = np.argsort(order)
		#print fit_table
		#print order
		reprod = []
		for i in xrange(len(fit_table)):
			N = len(fit_table)
			r = ranks[i]+1
			reprod.append(((c-1)/((c**N)-1))*c**(N-r))
		return reprod



#==============================================================================
#									TESTS


P = population(nb_genomes,nb_genes)
fitness_vect = []
iter_vect = []
avg_vect = []
RSS_vect = []
avg_vect = []
hier_vect = []
for i in xrange(100):
	fitness_vect.append(P.new_generation(i)[0])
	avg_vect.append(P.new_generation(i)[1])
	RSS_vect.append(P.new_generation(i)[1])
	hier_vect.append(P.new_generation(i)[1])
	iter_vect.append(i)

	csvfile = "Matrix"
	with open(csvfile,"w") as output:
		writer = csv.writer(output, lineterminator='\n')
		writer.writerows(P.pop[0].genome)

plt.plot(iter_vect,fitness_vect,marker='o')
plt.savefig("Evolution de la fitness")
plt.close()

plt.plot(iter_vect,avg_vect,marker='o')
plt.plot(iter_vect,math.log(nb_genes),marker='o')


plt.plot(iter_vect,RSS_vect,marker='o')
plt.plot(iter_vect,math.log(nb_genes),marker='o')


plt.plot(iter_vect,hier_vect,marker='o')
plt.plot(iter_vect,math.log(nb_genes),marker='o')

plt.savefig("Evolution de avg_shortest")
plt.close()

n = 50
a = r.rand(n,n) * r.choice([-1, 1], size=(n,n))
Gen = genome(n)



# for i in xrange(1):
# 	S0 = fitness.matrix_score(Gen.graph())
# 	print S0
# 	Gen.UpdateMatrix(S0)
# 	compteur +=1
# 	fitness.draw_figure_scalefree(Gen.graph(), fitness.score_matrix_scale_free(Gen.graph())[1],fitness. score_matrix_scale_free(Gen.graph())[2], compteur)
# 	#print compteur

