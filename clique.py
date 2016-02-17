import networkx as nx
import random
import matplotlib.pyplot as plt

genome = [[random.randint(0,1) for x in range(10)] for x in range(10)]

for i in genome:
	print i

G = nx.Graph()
G.add_nodes_from(range(10))

liste_edges = []
for i, vi in enumerate(genome[0]):
	for j, vj in enumerate(genome[0]):
		if i != j:
			if vi == 1 or vj == 1:
				liste_edges.append( (i,j) )

print liste_edges


G.add_edges_from(liste_edges)


pos = nx.spring_layout(G) # position of all nodes


# draw everything on the graph without any labels
plt.figure(1)
nx.draw(G,pos=nx.spring_layout(G))
plt.draw()

# to draw step by step the graph
plt.figure(2)
nx.draw_networkx_nodes(G,pos)

nx.draw_networkx_edges(G,pos, width=6, alpha=0.5, edge_color='b', style='dashed')
nx.draw_networkx_labels(G,pos, font_size=20, font_family='sans-serif')



plt.show()