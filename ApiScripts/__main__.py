from civicpy import civic
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

# url for human-phenotype-ontology
HPO_OBO_URL = "http://purl.obolibrary.org/obo/hp.obo"

GRAPH_FILENAME = "civic_data.gml"

# keys to save in the networkx structure
VARIANT_KEYS = [
	"evidenceAccepted",
	"evidenceSubmitted",
	"associatedPhenotypes",
	"variantName",
	"variantTypes"

]

# 58 corresponds to the VHL gene 
GENE_ID =  58


# given two nodes, finds the intersection over union for a specific attribute
def score_iou(n1, n2, attr_name):
	set1 = set(n1[attr_name])
	set2 = set(n2[attr_name])

	intersect = len(set1.intersection(set2))
	union = len(set1.union(set2))
	score = 0

	if not union == 0:
		score = intersect/union	

	return score

if __name__ == '__main__':
	parser = argparse.ArgumentParser()


	parser.add_argument('-d', '--directory', help = '''Local driectory to save graph''', default="")
	parser.add_argument('-u', '--update', help = '''Update the local variants from civic''', action="store_true")

	parser.add_argument('-i', '--include_submitted', help = '''If set, includes unreviewed variants''', action="store_true")

	args = parser.parse_args()

	variant_graph = nx.Graph()
	if args.update:
		gene = civic.get_gene_by_id(GENE_ID)

		#information relating to HPO terms
		# hpo = obonet.read_obo(HPO_OBO_URL)
		# assert nx.is_directed_acyclic_graph(hpo)
		# print('Number of HPO terms: {}'.format(len(hpo)))


		print('Number of variants for VHL: {}'.format(len(gene.variants)))
		for variant in gene.variants:
			variant_graph.add_node(variant.id)
			variant_node = variant_graph.nodes[variant.id]

			variant_node["variantName"] = variant.name

			#in the following code, "NA" has to be added as a placeholder element,
			#or else the grah cannot be parsed into gml

			#finding how many evidence items exist for the variant.
			#TODO: actually verify that evidence: supports, is germline, and is case study
			variant_node["evidenceAccepted"] = ["NA"]
			variant_node["evidenceSubmitted"] = ["NA"]
			variant_node["associatedPhenotypes"] = ["NA"]
			for evidence in variant.evidence_items:
				phenotypes = [phenotype.hpo_id for phenotype in evidence.phenotypes]
				if evidence.status == "submitted":
					variant_node["evidenceSubmitted"].append(evidence.id)
					if args.include_submitted:
						variant_node["associatedPhenotypes"].extend(phenotypes)
				elif evidence.status == "accepted":
					variant_node["evidenceAccepted"].append(evidence.id)	
					variant_node["associatedPhenotypes"].extend(phenotypes)

			#finding the types of the variant and adding it to the node
			variant_node["variantTypes"] = ["NA"]
			for variant_type in variant.variant_types:
				variant_node["variantTypes"].append(variant_type.name)				
			



		nx.write_gml(variant_graph, os.path.join(args.directory, GRAPH_FILENAME))



	variant_graph = nx.read_gml(os.path.join(args.directory, GRAPH_FILENAME))

	print("# of variants loaded from GML: {}".format(len(variant_graph)))
	#filtering out variants with no accepted evidence statements, if enabled



	to_remove = []
	for (node, nodeData) in variant_graph.nodes.items():
		#clean out all "NA"s from data

		for key, val in nodeData.items():

			if val == "NA":
				variant_graph.nodes[node][key] = []
			elif isinstance(val, list):
				variant_graph.nodes[node][key] = [item for item in val if not item == "NA"]


		if len(nodeData["evidenceAccepted"]) == 0:
			to_remove.append(node)


	if not args.include_submitted:
		variant_graph.remove_nodes_from(to_remove)
		print("# of variants after unreviewed variants removed: {}".format(len(variant_graph)))


	# vartype_graph = variant_graph.copy()
	# phenotype_graph = variant_graph.copy()
	#calculating similarity scores
	for (node1, node1Data) in variant_graph.nodes.items():
		for (node2, node2Data) in variant_graph.nodes.items():
			phen_weight = score_iou(node1Data, node2Data, "associatedPhenotypes")
			vartype_weight = score_iou(node1Data, node2Data, "variantTypes")
			weight = phen_weight+vartype_weight
			variant_graph.add_edge(node1, node2, phenotype_weight = phen_weight, vartype_weight = vartype_weight, weight = weight)

	filt_edges = [(u,v,d) for (u,v,d) in variant_graph.edges(data=True) if d["weight"] !=0] 
	edge_widths = [ d["weight"]/10 for (u,v,d) in filt_edges ]	
	edge_colors = [[d["phenotype_weight"], d["vartype_weight"],0, 0.5] for (u,v,d) in filt_edges]



	#graph drawing
	options = {
		'node_color': 'blue',
		'node_size': 50,
		'line_color': 'grey',
		'font_size': 6,
		'width':edge_widths,
		'edgelist': filt_edges,
		'edge_color': edge_colors,
		'labels': nx.get_node_attributes(variant_graph, 'variantName')
	}

	nx.draw(variant_graph, **options)

	plt.show()




	#making a dictionary of variant types
	variant_type_dict = {}
		# for var_type in variant.types:
		# 	if not var_type.name in variant_type_dict:
		# 		variant_type_dict[var_type.name] = []
		# 	variant_type_dict[var_type.name].append(variant.name)
			# using the actual civic ids instead of mutation type names
			# if not var_type.id in variant_type_dict:
			# 	variant_type_dict[var_type.id] = []
			# variant_type_dict[var_type.id].append(variant)

# QUESTIONS:
# - where do you find the keys each attribute type has?
# - is there a list of ids for mutation types (Civic Attribute)