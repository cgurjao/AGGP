import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as r

class genome:
	"""docstring for genome"""
	def __init__(self, number_gene):
		
		# The number of the genes we want
		self.nb_gene = number_gene

		# to create the random matrix for the initialisation, full of 0 and 1
		self.genome = r.randint(2,size=(self.nb_gene,self.nb_gene))

		# to symmetrize the matrix
		for i in xrange(self.nb_gene):
			for j in xrange(i):
				self.genome[i,j] = self.genome[j,i]

		

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
		
				if proba_matrix[i,j] < 0:
					if r.random() < -proba_matrix[i,j]:
						self.genome[i,j] = 0
						self.genome[i,j] = 0



	####################################################
	# function to draw the graph and show the graph
	####################################################
	def draw(self):
		G = nx.Graph()
		G.add_nodes_from(range(self.nb_gene))

		for i in xrange(self.nb_gene):
			for j in xrange(i+1):
				if self.genome[i,j] == 1:
					G.add_edge(i,j)



		G0 = sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]



		list_deg = np.array(nx.degree(G).values())
		position = nx.circular_layout(G)

		nx.draw_networkx_nodes(G, pos=position, node_size=100*list_deg, alpha=0.4)
		nx.draw_networkx_edges(G, pos=position, alpha=0.2)
		nx.draw_networkx_labels(G,pos=position)
		plt.draw()
		plt.show()



genome(10).draw()

