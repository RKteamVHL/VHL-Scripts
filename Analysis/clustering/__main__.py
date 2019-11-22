from .VariantGraph import VariantGraph
import math
import matplotlib.pyplot as plt
import networkx as nx
import argparse
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


	parser.add_argument('-d', '--directory', help = '''Local directory to save graph''', default="")
	parser.add_argument('-c', '--cache', help = '''Use the local graph cache''', action="store_true")

	parser.add_argument('-is', '--ignore_submitted', help = '''If set, ignores unreviewed variants (civic)''', action="store_true")

	args = parser.parse_args()

	VG = VariantGraph()

	#fetch/process all relevant data from all sources
	if not args.cache:

		# VG.add_nodes_from_db('Civic')
		VG.add_nodes_from_db('KimStudents2019')
		VG.add_nodes_from_db('Gnomad')
		VG.add_nodes_from_db('ClinVar')
		VG.merge_nodes()
		VG.save_to_json_file("variant_nodes.json")

	else:
		VG.load_from_json_file("variant_nodes.json")

	VG.calculate_node_attributes()
	VG.calculate_similarities()	
	VG.calculate_snf()
	#VG.remove_isolates()



	nx.draw(VG, 
		node_size=100, 
		labels=nx.get_node_attributes(VG, 'variantName'), 
		pos=nx.spring_layout(VG),
		width=0,
		node_color=[COLORMAP[d['spectral_label']] for (u,d) in VG.nodes(data=True)]
	)
	plt.show()

	VG.save_to_json_file("calculated_variant_nodes.json")


#this code is for finding a summary of al variant types
# 	type_dict = {}
# 	for node in VG.nodes(data=True):
# 		for vtype in node[1]["variantTypes"]:
# 			type_dict[vtype] =  type_dict.get(vtype, 0)+1

# 	print(type_dict)

# QUESTIONS:
# - where do you find the keys each attribute type has?
# - is there a list of ids for mutation types (Civic Attribute)