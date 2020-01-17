import math
import networkx as nx
import numpy as np
from . import variant_functions as vf
### The functions here perform variant-to-variant similarity analysis for pairs of nodes.
# Individual node analyses are in variant_functions.py

##Functions here should *only* take in n1 data, n2 data, and keyword args
# and must return a single numerical score


## Scoring functions
# Functions take in a dictionary for each node, and return a corresponding score

def score_iou(n1, n2, attr_name):
	'''Given two nodes, finds the intersection over union for a specific attribute
	The attributes here must be of a class type that is iterable (i.e., a list)
	'''
	set1 = set(n1['all'][attr_name])
	set2 = set(n2['all'][attr_name])

	intersect = len(set1.intersection(set2))
	union = len(set1.union(set2))
	score = 0

	if not union == 0:
		score = intersect/union	

	#find better way of creating weight name
	return score

def variant_score_domains(n1, n2, domain):
	'''Calculates similarity between variants based on whether they affect the
	same protein domains
	'''
	score = 0

	#if the inputted domain is in both variants
	if (domain in n1['all']['affected_domains']) and (domain in n2['all']['affected_domains']):
		score = 1
	#if the domain is in neither variant
	# elif not ((domain in n1['all']['affected_domains']) or (domain in n2['all']['affected_domains'])):
	# 	score = 0.1

	return score



def variant_aa_distance(n1, n2, sigma=5):
	'''Calculates similarity between snvs based on how far (in aa) they are
	Currently uses a gaussian distribution
	'''
	#TODO: get these from the nodes
	aa1 = 5
	aa2 = 10

	score = math.exp(-0.5 * (((aa1-aa2)/sigma )**2) )
	return score


def variant_nt_distance(n1, n2, sigma=200):
	'''Calculates similarity between snvs based on how far (in base pairs) they are
	Currently uses a gaussian distribution
	'''
	#TODO: get these from the nodes
	nt1 = 5
	nt2 = 10

	score = math.exp(-0.5 * (((nt1-nt2)/sigma )**2) )
	return score

def graph_distance(id1,id2, obo):
	"""Given two node ids and an obo, find their shortes distance
	"""
	distance = None
	if id1 == id2:
		distance = 0
	
	else:
		try:
			distance = nx.shortest_path_length(obo, source=id1, target=id2)
		except nx.exception.NetworkXNoPath as e:
			# print(repr(e))
			pass
	if isinstance(distance, int):
		distance += 1

	return distance 

# Note: an undirected graph is needed here for distance calculations
def variant_obo_distance(ids1, ids2, obo, sigma = 1 ):
	score = 0
	all_distances = []
	geo_mean = None

	longer_list = ids1 if len(ids1) >= len(ids2) else ids2
	shorter_list = ids1 if len(ids1) < len(ids2) else ids2
	for id1 in longer_list:
		distances = []
		for id2 in shorter_list:
			d = graph_distance(id1, id2, obo)
			if d is not None:
				distances.append(d)

		if len(distances) > 0:
			all_distances.append(min(distances))


	if len(all_distances) > 0:
		a = np.log(all_distances)
		geo_mean = np.exp(a.sum()/len(a))	
	
	if geo_mean is not None:
		score = math.exp(-0.5 * (((geo_mean-1)/sigma )**2) )

	return score

def variant_hpo_distance(n1, n2, sigma = 1):
	list1 = n1['all']['associatedPhenotypes']
	list2 = n2['all']['associatedPhenotypes']
	return variant_obo_distance(list1, list2, vf.OBONET_UD, sigma = 1)

def variant_so_distance(n1, n2, sigma = 1):
	list1 = n1['all']['variantTypes']
	list2 = n2['all']['variantTypes']
	return variant_obo_distance(list1, list2, vf.OBONET_UD, sigma = 1)
	