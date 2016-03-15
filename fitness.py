###########################################
#########  Import all libraries  ##########
###########################################
import networkx as nx
from random import *
import numpy as numpy
import math
import matplotlib.pyplot as plt
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    raise

from collections import deque
from itertools import chain, islice
try:
    from itertools import ifilter as filter
except ImportError:
    pass
import networkx
from networkx.utils.decorators import *

###########################################
#######  Define global parameters  ########
###########################################
##Define number of nodes
N = 100
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
#############  Draw graph ################## 
############################################

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

print "\nScore small-world!"
print "Average of all shortest paths is %f \n" %avg

##Calculate deviance of path length for each gene
deviance_path = []
gene_number = []
for g in nx.connected_component_subgraphs(G):
	for node in g:
        	    path_length=nx.single_source_dijkstra_path_length(g, node)
		    n=len(g)
    		    if n>1:
		    	deviance_path.append(sum(path_length.values())/float(n-1)-avg)
    		    else:
		    	deviance_path.append(0)
		    gene_number.append(node)

gene_number, deviance_path = zip(*sorted(zip(gene_number, deviance_path)))

##Score matrix according to deviance
deviance_matrix1 = numpy.zeros((N,N))
for i in range(0,N):
    for j in range(0,N):
		deviance_matrix1[i][j]=deviance_path[i] + deviance_path[j]

##Normalize data
if deviance_matrix1.max() != 0:
	norm_deviance_matrix1 = deviance_matrix1/deviance_matrix1.max()
else:
	norm_deviance_matrix1 = deviance_matrix1

##Create heatmap
plt.pcolor(norm_deviance_matrix1,cmap=plt.cm.Reds)
plt.show()
plt.close()


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
	
for i in xrange(len(unique_degrees)):
	observed[i] = observed[i]/ratio

##Associate deviance for each gene
deviance_degree = []
for i in xrange(len(degree_sequence)):
	for j in xrange(len(unique_degrees)):
		if (degree_sequence[i] == unique_degrees[j]):
			deviance_degree.append(observed[j]-expected[j])

deviance_matrix2 = numpy.zeros((N,N))
for i in range(0,N):
    for j in range(0,N):
		deviance_matrix2[i][j]=deviance_degree[i] + deviance_degree[j]

##Uncomment to draw figure with expected and observed curve of proportions of each degree
'''
plt.plot(unique_degrees, observed,marker='o')
plt.plot(unique_degrees, expected,marker='o')
plt.title("Degree proportion")
plt.ylabel("Proportion")
plt.xlabel("Degree")
plt.show()
'''
## draw graph in inset
'''
plt.axes([0.45,0.45,0.45,0.45])
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

##Normalize data
if deviance_matrix2.max() != 0:
	norm_deviance_matrix2 = deviance_matrix2/deviance_matrix2.max()
else:
	norm_deviance_matrix2 = deviance_matrix2

##Create heatmap
plt.pcolor(norm_deviance_matrix2,cmap=plt.cm.Greens)
plt.show()
plt.close()


#######################################################
# ----- Calculate the number of cliques in graph -----#
#######################################################
def enumerate_all_cliques(G): # list all cliques in G
  index = {}
  nbrs = {}
  for u in G:
      index[u] = len(index)
      # Neighbors of u that appear after u in the iteration order of G.
      nbrs[u] = {v for v in G[u] if v not in index}

  queue = deque(([u], sorted(nbrs[u], key=index.__getitem__)) for u in G)
  # Loop invariants:
  # 1. len(base) is nondecreasing.
  # 2. (base + cnbrs) is sorted with respect to the iteration order of G.
  # 3. cnbrs is a set of common neighbors of nodes in base.
  while queue:
      base, cnbrs = map(list, queue.popleft())
      yield base
      for i, u in enumerate(cnbrs):
          # Use generators to reduce memory consumption.
          queue.append((chain(base, [u]),
                        filter(nbrs[u].__contains__,
                               islice(cnbrs, i + 1, None))))


def nb_cliques(G) : # return the nomber of clique with minimum size k
  clique_list = []
  for clique in enumerate_all_cliques(G) :
    clique_list.append(clique)
  return len(clique_list)



###########################################
#---------- Clique score graph -----------#
###########################################
def Max_num_clique(N) : # Maximum number of cliques
  max_num_clique = 0
  for i in xrange(N):
    max_num_clique = max_num_clique + (math.factorial(N)/(math.factorial(i)*math.factorial(N-i)))
  return max_num_clique

def clique_score(G) :  # Clique Score 
  return float(nb_cliques(G)/float(Max_num_clique(N)))



###########################################
#---------- Matrix score graph -----------#
###########################################
# if m[i,j] < 0 : the interaction between i and j must be created if unexisting,
# or deleted if existing
# if m[i,j] > : the (un)interaction between i and j must be kept as it is

def matrix_score(G) :
  deviance_matrix3 = numpy.zeros((N,N))    # the matrix scores

  for i in range(0,N):
    for j in range(0,N) :
      score_obs = clique_score(G)   # score of the actual graph

      if G.has_edge(i,j):
        G.remove_edge(i,j)
      else :
        G.add_edge(i,j)

      score_inv = clique_score(G)   # score of the modified graph (interaction i-j added or deleted)

      deviance_matrix3[i,j] = float(score_obs - score_inv)

      if G.has_edge(i,j):
        G.remove_edge(i,j)
      else :
        G.add_edge(i,j)      
  return deviance_matrix3

deviance_matrix3 = matrix_score(G)
##Normalize data
if deviance_matrix3.max() != 0:
	norm_deviance_matrix3 = deviance_matrix3/deviance_matrix3.max()
else:
	norm_deviance_matrix3 = deviance_matrix3

##Create heatmap
plt.pcolor(norm_deviance_matrix3,cmap=plt.cm.Reds)
plt.show()
plt.close()


print "Score Clique!"
print "Average clique-number is %f \n" %clique_score(G)

