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
N = 100
##Scale-free network parameters
gamma = 2.2

###########################################
#-----------General functions-------------#
###########################################
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



###########################################
############## Scale-free #################
###########################################

def overall_RSS(G):
	degree_sequence = list(G.degree().values() - np.ones(len(G.degree().values())))
	##Sort by degrees
	RSS = 0.0

	##Sort by unique degrees
	unique_degrees=set(degree_sequence)
	#print unique_degrees
	## Calculate observed and expected values of proportions of each degree
	for i in unique_degrees:
		observed = degree_sequence.count(i)
		expected = math.pow(i, -gamma)
		RSS += abs(observed - expected)
	return RSS






def draw_figure_scalefree(G, observed, expected, compt = 0):
	degree_sequence = list(G.degree().values()) - np.ones(len(G.degree().values()))
	##Sort by degrees
	degree_sequence_sorted=sorted(degree_sequence,reverse=True)
	##Sort by unique degrees
	unique_degrees=sorted(set(degree_sequence),reverse=True)
	plt.plot(unique_degrees, observed,marker='o')
	plt.plot(unique_degrees, expected,marker='o')
	plt.title("Degree proportion")
	plt.ylabel("Proportion")
	plt.xlabel("Degree")
	## draw graph in inset
	plt.axes([0.45,0.45,0.45,0.45])
	Gcc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
	pos=nx.spring_layout(Gcc)
	plt.axis('off')
	nx.draw_networkx_nodes(Gcc,pos,node_size=20)
	nx.draw_networkx_edges(Gcc,pos,alpha=0.4)
	name = "Iteration %d .png" %compt
	plt.savefig(name)
	plt.close()







def fitness_score(G):
	avg = overall_average_shortest(G)
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
	#print fitness
	return fitness




# fitness_score(G)

