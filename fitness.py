###########################################
#########  Import all libraries  ##########
###########################################
import networkx as nx
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import deque
from itertools import chain, islice
try:
    from itertools import ifilter as filter
except ImportError:
    pass
from networkx.utils.decorators import *

###########################################
#######  Define global parameters  ########
###########################################
##Define number of nodes
N = 50
##Scale-free network parameters
gamma = 2.2

###########################################
#-----------General functions-------------#
###########################################
def generate_genome():
	##Genome with random numbers
	genome = np.reshape(np.random.random_integers(0,1,size=N*N),(N,N))
	##To symmetrize the matrix
	genome = (genome + genome.T)/2
	##Put all non value to 1 and fill the diagonal
	for i in xrange(N):
		for j in xrange(N):
			if genome[i][j] != 0:
				genome[i][j] = 1
			if i == j:
				genome[i][j] = 1
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
	return G, genome
	
def draw_heatmaps(matrix, save=False, name = "noname"):
	plt.pcolor(matrix,cmap=plt.cm.Reds)
	plt.show()	
	if save == True:	
		plt.savefig("Heatmaps.png")
	plt.close()

def draw_graph(save=False, name = "noname"):
	pos = nx.spring_layout(G)
	nx.draw(G, pos)
	node_labels = nx.get_node_attributes(G,'state')
	nx.draw_networkx_labels(G, pos, labels = node_labels)
	if save == True:	
		plt.savefig("simple_path.png")
	plt.show() 		

def normalize_data(matrix1):
	if matrix1.max() != 0:
		norm_matrix1 = matrix1/matrix1.max()
	else:
		norm_matrix1 = matrix1
	return norm_matrix1

###########################################
############### Small world ###############
###########################################

def overall_average_shortest(G):
	avg = 0.0
	for g in nx.connected_component_subgraphs(G):
		for node in g:
        		    path_length=nx.single_source_dijkstra_path_length(g, node)
        		    avg += sum(path_length.values())
    		n=len(g)
    		if n>1:
			    return float(avg)/float(n*(n-1))
    		else:
			    return 0.0
		avg = float(avg) + float(average_shortest(g))
	return avg			

def score_matrix_small_world(G):
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
	deviance_matrix1 = np.zeros((N,N))
	for i in range(0,N):
	    for j in range(0,N):
		deviance_matrix1[i][j]=deviance_path[i] + deviance_path[j]
	return deviance_matrix1

###########################################
############## Scale-free #################
###########################################

def overall_RSS(G):
	'''degree_sequence = [i-1 for i in G.degree().values()]
	##Sort by degrees
	RSS = 0.

	##Sort by unique degrees
	unique_degrees=set(degree_sequence)

	degree_sequence_sorted=sorted(degree_sequence,reverse=True)
	#print unique_degrees
	## Calculate observed and expected values of proportions of each degree
	for i in unique_degrees:
		observed = degree_sequence.count(i)
		expected = math.pow(i, -gamma)
		RSS += abs(observed - expected)

	for i in xrange(len(unique_degrees)):
		observed[i] = observed[i]/ratio'''

	degree_sequence = list(G.degree().values()) - np.ones(len(G.degree().values()))
	##Sort by degrees
	degree_sequence_sorted=sorted(degree_sequence,reverse=True)
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
	##To calculate proportionnality coefficient between observed and expected values
	ratio = 0.0
	for i in xrange(len(unique_degrees)):
		ratio = ratio + (observed[i]/expected[i])
	ratio = ratio/len(unique_degrees)
	for i in xrange(len(unique_degrees)):
		observed[i] = observed[i]/ratio

	for i in xrange(len(unique_degrees)):
		RSS = RSS + abs(observed[i] - expected[i])
	return RSS

def score_matrix_scale_free(G):
	degree_sequence = list(G.degree().values()) - np.ones(len(G.degree().values()))
	##Sort by degrees
	degree_sequence_sorted=sorted(degree_sequence,reverse=True)
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
	degree_sequence = list(G.degree().values()) - np.ones(len(G.degree().values()))
	unique_degrees=sorted(set(degree_sequence),reverse=True)
	##Associate deviance for each gene
	deviance_degree = []
	for i in xrange(len(degree_sequence)):
		for j in xrange(len(unique_degrees)):
			if (degree_sequence[i] == unique_degrees[j]):
				deviance_degree.append(observed[j]-expected[j])
	deviance_matrix2 = np.zeros((N,N))
	for i in range(0,N):
	    for j in range(0,N):
		deviance_matrix2[i][j]=deviance_degree[i] + deviance_degree[j]
	return deviance_matrix2, observed, expected

def draw_figure_scalefree(G,indiv = 0,compt = 0):
	degree_sequence = list(G.degree().values()) - np.ones(len(G.degree().values()))
	##Sort by degrees
	degree_sequence_sorted=sorted(degree_sequence,reverse=True)
	RSS = 0.0
	observed = []
	expected = []
	expected_curve = []
	##Sort by unique degrees
	unique_degrees=sorted(set(degree_sequence),reverse=True)
	#print unique_degrees
	## Calculate observed and expected values of proportions of each degree
	for i in xrange(len(unique_degrees)):
		observed.append(float(float(degree_sequence_sorted.count(unique_degrees[i]))))
		expected.append(math.pow(unique_degrees[i], -gamma))
		RSS = RSS + abs(observed[i] - expected[i])
	for i in range(1,100):
		expected_curve.append(math.pow(i, -gamma))

	##To calculate proportionnality coefficient between observed and expected values
	ratio = 0.0
	for i in xrange(len(unique_degrees)):
		ratio = ratio + (observed[i]/expected[i])
	ratio = ratio/len(unique_degrees)
	for i in xrange(len(unique_degrees)):
		observed[i] = observed[i]/ratio
	plt.plot(unique_degrees, observed,marker='o')
	plt.plot(np.arange(1, 100, 1), expected_curve,marker='o')
	plt.title("Individual %d" %indiv)
	plt.ylabel("Proportion")
	plt.xlabel("Degree")
	plt.xlim(10,50)
	plt.ylim(0,0.010)
	## draw graph in inset
	plt.axes([0.45,0.45,0.45,0.45])
	Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
	pos=nx.spring_layout(Gcc)
	plt.axis('off')
	nx.draw_networkx_nodes(Gcc,pos,node_size=20)
	nx.draw_networkx_edges(Gcc,pos,alpha=0.4)
	name = "Indiv_%d _Iteration_%d .png" %(indiv, compt)
	plt.savefig(name)
	plt.close()

###########################################
##############   Clique   #################
###########################################
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

def Max_num_clique(N) : # Maximum number of cliques
  max_num_clique = 0
  for i in xrange(N):
    max_num_clique = max_num_clique + (math.factorial(N)/(math.factorial(i)*math.factorial(N-i)))
  return max_num_clique

def overall_clique_score(G) :  # Clique Score 
  return float(nb_cliques(G)/float(Max_num_clique(N)))

# if m[i,j] < 0 : the interaction between i and j must be created if unexisting,
# or deleted if existing
# if m[i,j] > : the (un)interaction between i and j must be kept as it is
def matrix_score(G) :
  deviance_matrix3 = np.zeros((N,N))    # the matrix scores

  for i in range(0,N):
    for j in range(0,N) :
      score_obs = overall_clique_score(G)   # score of the actual graph

      if G.has_edge(i,j):
        G.remove_edge(i,j)
      else :
        G.add_edge(i,j)

      score_inv = overall_clique_score(G)   # score of the modified graph (interaction i-j added or deleted)

      deviance_matrix3[i,j] = float(score_obs - score_inv)

      if G.has_edge(i,j):
        G.remove_edge(i,j)
      else :
        G.add_edge(i,j)      
  return deviance_matrix3

###########################################
###############   TEST   ##################
###########################################
#Generate graph
G = generate_genome()[0]
#Generate genome
#genome = generate_genome()[1]

############ Small-World ##################
#Average score for small-world parameter
#avg = overall_average_shortest(G)
#print "\nScore small-world!"
#print "Average of all shortest paths is %f \n" %avg
#Score matrix for small-world parameter
#deviance_matrix1 = score_matrix_small_world(G)
#norm_deviance_matrix1 = normalize_data(deviance_matrix1[0])
############ Scale-free ###################
#Average score for scale-free parameter
#RSS_score = overall_RSS(G)
#print "Score Scale-free!"
#print "Root Sum Square is %f \n" %RSS_score
#Score matrix for small-world parameter
#deviance_matrix2 = score_matrix_scale_free(G)[0]
#norm_deviance_matrix2 = normalize_data(deviance_matrix2)
#Draw figure for scale-free
#observed = score_matrix_scale_free(G)[1]
#expected = score_matrix_scale_free(G)[2]
#draw_figure_scalefree(G, observed, expected)

############### Clique ###################
#Average score for scale-free parameter
#score_clique = overall_clique_score(G)
#print "Score Clique!"
#print "Average clique-number is %f \n" %score_clique
#Score matrix for small-world parameter
#deviance_matrix3 = matrix_score(G)
#norm_deviance_matrix3 = normalize_data(deviance_matrix3)

###### Output of this algorithm ##########
#Overall_score_matrix = norm_deviance_matrix1/3 +norm_deviance_matrix2/3 + norm_deviance_matrix3/3
#draw_heatmaps(Overall_score_matrix)

def fitness_score(G):
	avg = 0#overall_average_shortest(G)
	#print "\nScore small-world!"
	#print "Average of all shortest paths is %f \n" %avg
	RSS_score = overall_RSS(G)
	#print "Score Scale-free!"
	#print "Root Sum Square is %f \n" %RSS_score
	#score_clique = overall_clique_score(G)
	#print "Score Clique!"
	#print "Average clique-number is %f \n" %score_clique
	fitness = 0
	if avg+RSS_score != 0:
		fitness = 1/(avg+RSS_score)
	return fitness

