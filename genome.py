# coding=utf-8
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

open('fitness.txt', 'w').close()
open('small-world.txt', 'w').close()
open('scale-free.txt', 'w').close()
open('cluster.txt', 'w').close()


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
			res = fitness.fitness_score(g.graph())
			temp_avg.append(float(res[0]))
			temp_RSS.append(res[1])
			temp_hier.append(res[2])

			if rd!=i:
				g.genome, self.pop[rd].genome = g.CrossingOver(g.genome, self.pop[rd].genome)   # Do random crossing over
			print "Iteration: %d, Individu: %d" %(compt, i)
			if i == 0:
				fitness.draw_figure_scalefree(g.graph(), i, compt)
				fitness.draw_figure_hierarchical(g.graph(), i, compt)
				#g.draw(g.graph(), compt, i)

		ratio = 0.0
		for i in xrange(nb_genomes):
			if temp_RSS[i] != 0:
				ratio = ratio + float(temp_avg[i]/temp_RSS[i])
		ratio = ratio/float(nb_genomes)
		for i in xrange(nb_genomes):
			if ratio != 0:
				temp_avg[i] = float(temp_avg[i])/ratio
		ratio = 0
		for i in xrange(nb_genomes):
			if temp_RSS[i] != 0:
				ratio = ratio + (temp_hier[i]/temp_RSS[i])
		ratio = ratio/float(nb_genomes)
		for i in xrange(nb_genomes):
			if ratio != 0:
				temp_hier[i] = temp_hier[i]/ratio
		for i in xrange(nb_genomes):
			F = map(add,temp_avg, temp_RSS)
			F = map(add,F, temp_hier)

		with open("fitness.txt", "a") as myfile:
			for i in range(len(F)):
				myfile.write(str(F[i]))
				myfile.write("\t")
			myfile.write("\n")

		with open("small-world.txt", "a") as myfile:
			for i in range(len(F)):
				myfile.write(str(temp_avg[i]))
				myfile.write("\t")
			myfile.write("\n")

		with open("scale-free.txt", "a") as myfile:
			for i in range(len(F)):
				myfile.write(str(temp_RSS[i]))
				myfile.write("\t")
			myfile.write("\n")


		with open("cluster.txt", "a") as myfile:
			for i in range(len(F)):
				myfile.write(str(temp_hier[i]))
				myfile.write("\t")
			myfile.write("\n")

		proba_rep = self.selection(F,selection_intensity)
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
		return sum(F), temp_avg, temp_RSS, temp_hier
			
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
	res = P.new_generation(i)
	fitness_vect.append(res[0])
	if res[1] != 0:
		avg_vect.append(1.0/res[1][0])
	if res[2] != 0:
		RSS_vect.append(1.0/res[2][0])
	if res[3] != 0:
		hier_vect.append(1.0/res[3][0])
	iter_vect.append(i)

	#csvfile = "Matrix"
	#with open(csvfile,"w") as output:
		#writer = csv.writer(output, lineterminator='\n')
		#writer.writerows(P.pop[0].genome)

if len(iter_vect) == len(fitness_vect):
	plt.plot(iter_vect,fitness_vect,marker='o')
	plt.savefig("Evolution de la fitness")
	plt.ylabel("Fitness")
	plt.xlabel("Itération")
	plt.close()

if len(iter_vect) == len(avg_vect):
	plt.plot(iter_vect,avg_vect,marker='o')
	plt.savefig("Evolution du score petit-monde")
	plt.ylabel("Score petit-monde")
	plt.xlabel("Itération")
	plt.close()

if len(iter_vect) == len(RSS_vect):
	plt.plot(iter_vect,RSS_vect,marker='o')
	plt.savefig("Evolution du score d'invariance d'échelle")
	plt.ylabel("Score d'invariance d'échelle")
	plt.xlabel("Itération")
	plt.close()


if len(iter_vect) == len(hier_vect):
	plt.plot(iter_vect,hier_vect,marker='o')
	plt.savefig("Evolution du score de clustering")
	plt.ylabel("Score de clustering")
	plt.xlabel("Itération")
	plt.close()


table = np.loadtxt('fitness.txt', skiprows=1)
fig, ax = plt.subplots(1)
p = ax.pcolormesh(table)
fig.colorbar(p)
fig.savefig('Fitness.png')


table = np.loadtxt('cluster.txt', skiprows=1)
fig, ax = plt.subplots(1)
p = ax.pcolormesh(table)
fig.colorbar(p)
fig.savefig('Cluster.png')


table = np.loadtxt('small-world.txt', skiprows=1)
fig, ax = plt.subplots(1)
p = ax.pcolormesh(table)
fig.colorbar(p)
fig.savefig('Small-World.png')


table = np.loadtxt('scale-free.txt', skiprows=1)
fig, ax = plt.subplots(1)
p = ax.pcolormesh(table)
fig.colorbar(p)
fig.savefig('Scale-Free.png')
# for i in xrange(1):
# 	S0 = fitness.matrix_score(Gen.graph())
# 	print S0
# 	Gen.UpdateMatrix(S0)
# 	compteur +=1
# 	fitness.draw_figure_scalefree(Gen.graph(), fitness.score_matrix_scale_free(Gen.graph())[1],fitness. score_matrix_scale_free(Gen.graph())[2], compteur)
# 	#print compteur

