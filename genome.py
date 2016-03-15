import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as r

class genome:
	def __init__(self, number_gene):
		
		# The number of the genes we want
		self.nb_gene = number_gene

		# to create the random matrix for the initialisation, full of 0 and 1
		self.genome = r.randint(2,size=(self.nb_gene,self.nb_gene))

		# to symmetrize the matrix
		for i in xrange(self.nb_gene):
			for j in xrange(i):
				self.genome[i,j] = self.genome[j,i]


		# to construct the graph
		self.G = nx.Graph()
		self.G.add_nodes_from(range(self.nb_gene))

		for i in xrange(self.nb_gene):
			for j in xrange(i+1):
				if self.genome[i,j] == 1:
					self.G.add_edge(i,j)

		self.position = nx.circular_layout(self.G)
		

	####################################################
	# function to update the matrix
	####################################################
	def UpdateMatrix(self, proba_matrix):
		for i in xrange(self.nb_gene):
			for j in xrange(i):

				if proba_matrix[i,j] > 0:
					if r.random() < proba_matrix[i,j]:
						self.genome[i,j] = 1
						self.genome[j,i] = 1
						self.G.add_edge(i,j)
		
				if proba_matrix[i,j] < 0:
					if r.random() < -proba_matrix[i,j]:
						self.genome[i,j] = 0
						self.genome[i,j] = 0
						if (i,j) in self.G.edges() or (j,i) in self.G.edges():
							self.G.remove_edge(i,j)


	####################################################
	# function to draw the graph and show the graph
	####################################################
	def draw(self):
		plt.figure()
		list_deg = np.array(nx.degree(self.G).values())

		nx.draw_networkx_nodes(self.G, pos=self.position, node_size=300, alpha=0.5, node_color=list_deg, cmap=plt.cm.YlOrRd)# RdPu, BuPu
		# http://matplotlib.org/examples/color/colormaps_reference.html
		nx.draw_networkx_edges(self.G, pos=self.position, alpha=0.2)
		nx.draw_networkx_labels(self.G,pos=self.position)
		plt.draw()


n = 100
a = r.rand(n,n) * r.choice([-1, 1], size=(n,n))
G = genome(n)
G.draw()
G.UpdateMatrix(a)

G.draw()
plt.show()