import networkx as nx
from random import *
import numpy as numpy
try:
    import matplotlib.pyplot as plt
except:
    raise
###########################################
# Generate random genome
###########################################

genome = [[randint(0,1) for x in range(5)] for x in range(5)]
for i in xrange(len(genome[0])):
	print (genome[i])

###########################################
# Generate nx-compatible graph
###########################################
genome = numpy.array(genome)
G = nx.Graph(genome)

###########################################
# Average of shortest path
###########################################
D = nx.average_shortest_path_length(G)
print("Average of all path is")
print(D)

###########################################
# Manipulate degrees
###########################################

degree_sequence = list(G.degree().values())
print(degree_sequence)


#############################################
# Draw graph
#############################################
nx.draw(G)
plt.savefig("simple_path.png") # save as png
plt.show() 


#degree_sequence=sorted(nx.degree(G).values(),reverse=True) # degree sequence
#print "Degree sequence", degree_sequence
dmax=max(degree_sequence)

plt.plot(degree_sequence,'b-',marker='o')
plt.title("Degree rank plot")
plt.ylabel("degree")
plt.xlabel("rank")

# draw graph in inset
plt.axes([0.45,0.45,0.45,0.45])
Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
pos=nx.spring_layout(Gcc)
plt.axis('off')
nx.draw_networkx_nodes(Gcc,pos,node_size=20)
nx.draw_networkx_edges(Gcc,pos,alpha=0.4)

#plt.savefig("degree_histogram.png")
plt.show()


