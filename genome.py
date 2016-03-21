import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as r
import fitness

class genome:
	def __init__(self, number_gene):
		
		# The number of the genes we want
		self.nb_gene = number_gene

		# to create the random matrix for the initialisation, full of 0 and 1
		self.genome = r.randint(2,size=(self.nb_gene,self.nb_gene))

		# to symmetrize the matrix
		for i in xrange(self.nb_gene):
			for j in xrange(i):
				self.genome[i,j] = 0
			self.genome[i,i] = 1
		print self.genome

		self.nb_fig = 1

	####################################################
	# function to update the matrix
	####################################################
	def UpdateMatrix(self, score_matrix):
		for i in xrange(self.nb_gene):
			for j in xrange(i):
				if score_matrix[j,i] < 0:						# If score is negative (bad) -> mutation!
					self.genome[j,i] = abs(self.genome[j,i]-1)

				if score_matrix[j,i] == 0:                      # If score is null-> mutation half of the time
					if r.rand() < 0.5:
						self.genome[j,i] = abs(self.genome[i,j]-1)

				# If score is positive (good) -> no mutation

	####################################################
	# function to draw the graph
	####################################################
	def graph(self):
		# to construct the graph
		G = nx.Graph()
		G.add_nodes_from(range(self.nb_gene))

		for i in xrange(self.nb_gene):
			for j in xrange(i+1):
				if self.genome[j,i] == 1:
					G.add_edge(i,j)
		return G


	def draw(self, graph, itteration=False, save=False):
		# if you use the option save = True, you will save also all the past graph but not the future graph

		
		position = nx.circular_layout(G)


		fig = plt.figure(self.nb_fig)
		self.nb_fig += 1

		list_deg = np.array(nx.degree(G).values())

		nx.draw_networkx_nodes(G, pos=position, node_size=300, alpha=0.5, node_color=list_deg, cmap=plt.cm.YlOrRd)# RdPu, BuPu
		# http://matplotlib.org/examples/color/colormaps_reference.html
		nx.draw_networkx_edges(G, pos=position, alpha=0.2)
		nx.draw_networkx_labels(G,pos=position)
		plt.draw()

		if itteration != False and save == True:

			plt.title('graph for itteration : %d'%itteration)
			plt.savefig('graph_itteration_%d'%itteration)

		if save == True and itteration == False:
			plt.title('Graph')
			plt.savefig('graph')




n = 10
a = r.rand(n,n) * r.choice([-1, 1], size=(n,n))
Gen = genome(n)

#plt.show()

compteur = 0

for i in xrange(10):
	S0 = fitness.matrix_score(Gen.graph())
	Gen.UpdateMatrix(S0)
	compteur +=1

print "\n",Gen.genome