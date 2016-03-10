import networkx as nx
from random import *
import numpy as numpy

N = 5#Define number of nodes

try:
    import matplotlib.pyplot as plt
except:
    raise
###########################################
######### Generate random genome ##########
###########################################
##Genome with random numbers
genome = numpy.reshape(numpy.random.random_integers(0,1,size=N*N),(N,N))
##To symmetrize the matrix
#genome = genome + genome.T - numpy.diag(genome.diagonal())
genome = (genome + genome.T)/2
##Put all non value to 1 and fill the diagonal
for i in xrange(N):
	for j in xrange(N):
		if genome[i][j] != 0:
			genome[i][j] = 1
		if i == j:
			genome[i][j] = 1
		
#Custom genome for tests
#genome = [[1, 0, 0, 0, 0],[0, 1, 0, 1, 1],[0, 0, 1, 1, 1],[0, 1, 1, 1, 0],[0, 1, 1, 0, 1]]


###########################################
###### Generate nx-compatible graph #######
###########################################
##Create empty graph
G = nx.Graph()
##Create nodes and label them (here, their column number)
for i in xrange(N):
	G.add_node(i)
	G.node[i]['state']= i
##Add edge when two nodes are connected
for i in xrange(N):
	for j in xrange(i+1):
		if genome[i][j] != 0:
			G.add_edge(i,j)

############################################
# Draw graph
############################################
#print genome
#print G.number_of_edges()
#print G.edges()

pos = nx.spring_layout(G)
nx.draw(G, pos)
node_labels = nx.get_node_attributes(G,'state')
nx.draw_networkx_labels(G, pos, labels = node_labels)
#plt.savefig("simple_path.png")
plt.show() 			

###########################################
######## Average of shortest path #########
###########################################
## It is the sum of all shortest paths divided by the number of paths (which is n!)

def average_shortest(g):
	avg = 0.0
	for node in g:
        	    path_length=nx.single_source_dijkstra_path_length(g, node)
        	    avg += sum(path_length.values())
    	n=len(g)
    	if n>1:
		    return float(avg)/float(n*(n-1))
    	else:
		    return 0.0

avg = 0.0
for g in nx.connected_component_subgraphs(G):
	avg = float(avg) + float(average_shortest(g))

print "Average of all shortest paths is %f" %avg

###########################################
########## Manipulate degrees #############
###########################################

print G.edges()
print G.degree(0)
degree_sequence = list(G.degree().values()) - numpy.ones(len(G.degree().values()))
print(degree_sequence)

#degree_sequence=sorted(nx.degree(G).values(),reverse=True) # degree sequence
#print "Sorted degree sequence", degree_sequence
#dmax=max(degree_sequence)

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


############################################
# To have the frequence list
############################################
list_deg = nx.degree(G).values()
freq = {}
for i in set(list_deg):
	freq[i] = list_deg.count(i)
print freq

