import networkx as nx
from random import *
import numpy as numpy
import math
try:
    import matplotlib.pyplot as plt
except:
    raise

##Define number of nodes
N = 500
##Scale-free network parameters
gamma = 2.2

###########################################
######### Generate random genome ##########
###########################################
##Genome with random numbers
genome = numpy.reshape(numpy.random.random_integers(0,1,size=N*N),(N,N))
##To symmetrize the matrix
genome = (genome + genome.T)/2
##Put all non value to 1 and fill the diagonal
for i in xrange(N):
	for j in xrange(N):
		if genome[i][j] != 0:
			genome[i][j] = 1
		if i == j:
			genome[i][j] = 1
		
##Custom genome for tests
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

##Uncomment to draw graph
'''
pos = nx.spring_layout(G)
nx.draw(G, pos)
node_labels = nx.get_node_attributes(G,'state')
nx.draw_networkx_labels(G, pos, labels = node_labels)
##Uncomment to save figure
#plt.savefig("simple_path.png")
plt.show() 			
'''

###########################################
######## Average of shortest path #########
###########################################
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
print "\n Score small-world!"
print "Average of all shortest paths is %f \n" %avg

###########################################
########## Manipulate degrees #############
###########################################
degree_sequence = list(G.degree().values()) - numpy.ones(len(G.degree().values()))

## Check-up for tests
#print G.edges()
#print G.degree(0)
#print(degree_sequence)


##Sort by degrees
degree_sequence_sorted=sorted(degree_sequence,reverse=True)
#print "Degree sequence"
#print degree_sequence_sorted

##Uncomment to draw histogram of degrees with graph on the same figure
"""
plt.plot(degree_sequence,'b-',marker='o')
plt.title("Degree rank plot")
plt.ylabel("degree")
plt.xlabel("rank")
## draw graph in inset
plt.axes([0.45,0.45,0.45,0.45])
Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
pos=nx.spring_layout(Gcc)
plt.axis('off')
nx.draw_networkx_nodes(Gcc,pos,node_size=20)
nx.draw_networkx_edges(Gcc,pos,alpha=0.4)
##Uncomment to save figure
#plt.savefig("degree_histogram.png")
plt.show()
"""

##Uncomment to draw figure with expected and observed curve of proportions of each degree
'''
plt.plot(degree_sequence,'b-',marker='o')
plt.title("Degree rank plot")
plt.ylabel("degree")
plt.xlabel("rank")
## draw graph in inset
plt.axes([0.45,0.45,0.45,0.45])
Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
pos=nx.spring_layout(Gcc)
plt.axis('off')
nx.draw_networkx_nodes(Gcc,pos,node_size=20)
nx.draw_networkx_edges(Gcc,pos,alpha=0.4)
##Uncomment to save figure
#plt.savefig("degree_histogram.png")
plt.show()
'''

RSS = 0.0
observed = []
expected = []
##Sort by unique degrees
unique_degrees=sorted(set(degree_sequence),reverse=True)
#print unique_degrees

## Calculate observed and expected values of proportions of each degree
for i in xrange(len(unique_degrees)):
	observed.append(float(float(degree_sequence_sorted.count(unique_degrees[i]))))
	expected.append(math.pow(unique_degrees[i], -gamma))
	RSS = RSS + abs(observed[i] - expected[i])


##To calculate proportionnality coefficient between observed and expected values
ratio = 0.0
for i in xrange(len(unique_degrees)):
	ratio = ratio + (observed[i]/expected[i])

ratio = ratio/len(unique_degrees)
print "Ratio"
print ratio
	
for i in xrange(len(unique_degrees)):
	expected[i] = expected[i]*ratio

##Uncomment to draw figure with expected and observed curve of proportions of each degree
plt.plot(unique_degrees, observed,marker='o')
plt.plot(unique_degrees, expected,marker='o')
plt.title("Degree proportion")
plt.ylabel("Proportion")
plt.xlabel("Degree")
plt.show()
## draw graph in inset
'''plt.axes([0.45,0.45,0.45,0.45])
Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
pos=nx.spring_layout(Gcc)
plt.axis('off')
nx.draw_networkx_nodes(Gcc,pos,node_size=20)
nx.draw_networkx_edges(Gcc,pos,alpha=0.4)
##Uncomment to save figure
#plt.savefig("degree_histogram.png")
'''

print "Score Scale-free!"
print "Root Sum Square is %f \n" %RSS

###########################################
########## Manipulate cliques #############
###########################################

