from .VariantGraph import VariantGraph
from . import constants
import math
import matplotlib.pyplot as plt
import networkx as nx
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
		VG.save_to_json_file("variant_nodes.json")

	else:
		VG.load_from_json_file("variant_nodes.json")

	VG.calculate_node_attributes()
	VG.calculate_similarities()
	type_dict = {}
	for node in VG.nodes(data=True):
		for vtype in node[1]["variantTypes"]:
			type_dict[vtype] =  type_dict.get(vtype, 0)+1

	print(type_dict)
	VG.save_to_json_file("calculated_variant_nodes.json")




# QUESTIONS:
# - where do you find the keys each attribute type has?
# - is there a list of ids for mutation types (Civic Attribute)