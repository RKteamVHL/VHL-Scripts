from .VariantGraph import VariantGraph
from . import constants
import math
import matplotlib.pyplot as plt
import networkx as nx
import obonet
import argparse
import time
import json
import csv
import os 

import pprint


if __name__ == '__main__':
	parser = argparse.ArgumentParser()


	parser.add_argument('-d', '--directory', help = '''Local directory to save graph''', default="")
	parser.add_argument('-u', '--update', help = '''Update the local graph cache''', action="store_true")

	parser.add_argument('-is', '--ignore_submitted', help = '''If set, ignores unreviewed variants (civic)''', action="store_true")

	args = parser.parse_args()

	VG = VariantGraph()

	#fetch/process all relevant data from all sources
	if args.update:

		VG.add_nodes_from_civic_by_gene(constants.VHL_GENE_ID, ignore_submitted=args.ignore_submitted)	

		VG_json = nx.readwrite.json_graph.node_link_data(VG)

		#TODO: make these methods of VariantGraph		
		with open(os.path.join(args.directory, "variant_"+constants.GRAPH_FILENAME+ ".json"), "w") as file:
			json.dump(VG_json, file)

		print("# of variants saved to file: {}".format(len(VG)))

	else:
		#TODO: make these methods of VariantGraph
		with open(os.path.join(args.directory, "variant_"+constants.GRAPH_FILENAME+ ".json"), "r") as file:
			VG_json = json.load(file)
			VG = nx.readwrite.json_graph.node_link_graph(VG_json)	

		print("# of variants loaded from file: {}".format(len(VG)))


	VG_nodes = list(VG.nodes(data=True))

	#assuming undirected graphs with no self connections
	for i in range(0,len(VG_nodes)):
		for j in range(i+1, len(VG_nodes)):

			similarities = {}
			for metric in constants.SIMILARITY_METRICS:
				(name, score) = metric["function"](VG_nodes[i][1], VG_nodes[j][1], **metric["kwargs"])
				similarities[name] = score

			#TODO: dont add edge if all scores are 0
			VG.add_edge(VG_nodes[i][0], VG_nodes[j][0], **similarities)

			
	VG.save_to_json_file()
	# filt_edges = [(u,v,d) for (u,v,d) in variant_graph.edges(data=True) if d["weight"] !=0] 
	# edge_widths = [ d["weight"]/10 for (u,v,d) in filt_edges ]	
	# edge_colors = [[d["phenotype_weight"], d["vartype_weight"],0, 0.5] for (u,v,d) in filt_edges]



	# #graph drawing
	# options = {
	# 	'node_color': 'blue',
	# 	'node_size': 50,
	# 	'line_color': 'grey',
	# 	'font_size': 6,
	# 	'width':edge_widths,
	# 	'edgelist': filt_edges,
	# 	'edge_color': edge_colors,
	# 	'labels': nx.get_node_attributes(variant_graph, 'variantName')
	# }

	# nx.draw(variant_graph, **options)

	# plt.show()



# QUESTIONS:
# - where do you find the keys each attribute type has?
# - is there a list of ids for mutation types (Civic Attribute)