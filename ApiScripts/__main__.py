from civicpy import civic
import hgvs.parser
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
	"associatedPhenotypes"

]

# 58 corresponds to the VHL gene 
GENE_ID =  58

# should include submitted, along with accepted, evidence statements 
INCLUDE_SUBMITTED = False

if __name__ == '__main__':
	parser = argparse.ArgumentParser()


	parser.add_argument('-d', '--directory', help = '''Local driectory to save graph''', default="")
	parser.add_argument('-u', '--update', help = '''Update the local variants from civic''', action="store_true")

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

			#finding how many evidence items exist for the variant.
			#TODO: actually verify that evidence: supports, is germline, and is case study
			variant_node["evidenceAccepted"] = ["NA"]
			variant_node["evidenceSubmitted"] = ["NA"]
			variant_node["associatedPhenotypes"] = ["NA"]
			for evidence in variant.evidence_items:
				phenotypes = [phenotype.hpo_id for phenotype in evidence.phenotypes]
				if evidence.status == "submitted":
					variant_node["evidenceSubmitted"].append(evidence.id)
					if INCLUDE_SUBMITTED:
						variant_node["associatedPhenotypes"].extend(phenotypes)
				elif evidence.status == "accepted":
					variant_node["evidenceAccepted"].append(evidence.id)	
					variant_node["associatedPhenotypes"].extend(phenotypes)

						
			



		nx.write_gml(variant_graph, os.path.join(args.directory, GRAPH_FILENAME))

	else:

		variant_graph = nx.read_gml(os.path.join(args.directory, GRAPH_FILENAME))


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

		pprint.pprint(variant_graph)

# QUESTIONS:
# - where do you find the keys each attribute type has?
# - is there a list of ids for mutation types (Civic Attribute)