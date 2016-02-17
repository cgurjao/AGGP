import networkx as nx
import random
genome = [[random.randint(0,1) for x in range(10)] for x in range(10)]
for i in genome:
	print i

G = nx.graph()
G.add_nodes_from(range(10))

