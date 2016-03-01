import networkx as nx
import random
import numpy
import matplotlib.pyplot as plt

a = numpy.reshape(numpy.random.random_integers(0,1,size=100),(10,10))



print a

G = nx.from_numpy_matrix(a)

# print liste_edges

# the list of edges must be a tuple with the two labels of the nodes inside


pos = nx.spring_layout(G) # position of all nodes


# draw everything on the graph without any labels
plt.figure(1)
nx.draw(G,pos=nx.spring_layout(G))
plt.draw()



# to draw step by step the graph
plt.figure(2)
nx.draw_networkx_nodes(G,pos, node_size=700, node_color=range(10))

nx.draw_networkx_edges(G,pos, width=6, alpha=0.5, edge_color='b', style='dashed')
nx.draw_networkx_labels(G,pos, font_size=20, font_family='sans-serif')



# plt.show()




# if wants a weigthed graph you can add in the tuple a third number witch is the weigth of the edge
# and with the function             add_weighted_edges_from(list_edges_weighted)


# draw the histogramm of the degree of nodes
# plt.figure(3)
print nx.degree(G).values()
# plt.hist(nx.degree(G).values())



list_deg = nx.degree(G).values()

freq = {}

for i in list_deg:

	if i in freq.keys():
		freq[i] += 1

	else:
		freq[i] = 1

print freq

plt.show()
