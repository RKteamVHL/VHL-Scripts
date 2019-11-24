from .VariantGraph import VariantGraph
import math
import matplotlib.pyplot as plt
import networkx as nx
import argparse
from snf import compute
import time
import json
import csv
import os 

import pprint

COLORMAP= {
	0:(0,0,0,0.5), 
	1:(0,0,1,0.5),
	2:(0,1,0,0.5),
	3:(0,1,1,0.5),
	4:(0,1,0,0.5),
	5:(0,1,1,0.5),
	6:(1,0,0,0.5),
	7:(1,0,1,0.5),

	
}

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('similarity_type', help = '''The type of similarity to base clustering on''')
	parser.add_argument('-d', '--directory', help = '''Local directory to save graph''', default="")
	parser.add_argument('-rc', '--pre_cache', help = '''Use the local graph raw node cache''', action="store_true")
	parser.add_argument('-oc', '--post_cache', help = '''Use the local graph calulated node cache''', action="store_true")
	parser.add_argument('-is', '--ignore_submitted', help = '''If set, ignores unreviewed variants (civic)''', action="store_true")

	args = parser.parse_args()

	VG = VariantGraph()

	#fetch/process all relevant data from all sources
	if not args.pre_cache:

		VG.add_nodes_from_db('Civic')
		VG.add_nodes_from_db('KimStudents2019')
		VG.add_nodes_from_db('Gnomad')
		VG.add_nodes_from_db('ClinVar')
		VG.merge_nodes()
		VG.save_to_json_file("variant_nodes.json")

	else:
		VG.load_from_json_file("variant_nodes.json")

	PG = VariantGraph()
	PG.add_nodes_from([(n, d) for n, d in VG.nodes(data=True) if d['all']['associatedPhenotypes']])
	PG.add_edges_from([(n1, n2, d) for n1, n2, d in VG.edges(data=True) if n1 in PG and n2 in PG])


	PG.save_to_json_file("pheno_variant_nodes.json")

	if not args.post_cache:
		PG.calculate_node_attributes()
		PG.calculate_similarities()
		PG.save_to_json_file("calculated_variant_nodes.json")	
	else:
		PG.load_from_json_file("calculated_variant_nodes.json")

	adj_mat = PG.get_adjacency_mats(types=[args.similarity_type])

	clust_count1, clust_count2 = compute.get_n_clusters(adj_mat[0])
	print("Cluster estimates: {}, {}".format(clust_count1, clust_count2))

	PG.cluster_by(args.similarity_type, num_clusters=4)
	PG.save_to_json_file("labeled_variant_nodes.json", nodes_only=True)
	#VG.remove_isolates()



	nx.draw(PG, 
		node_size=100, 
		with_labels=True,		
		pos=nx.spring_layout(PG),
		width=0.001,
		node_color=[COLORMAP[d[f'{args.similarity_type}_label']] for (u,d) in PG.nodes(data=True)]
	)
	plt.show()

	


#this code is for finding a summary of al variant types
# 	type_dict = {}
# 	for node in VG.nodes(data=True):
# 		for vtype in node[1]["variantTypes"]:
# 			type_dict[vtype] =  type_dict.get(vtype, 0)+1

# 	print(type_dict)

# QUESTIONS:
# - where do you find the keys each attribute type has?
# - is there a list of ids for mutation types (Civic Attribute)