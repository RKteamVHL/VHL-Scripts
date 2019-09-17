from civicpy import civic
from . import constants
import networkx as nx
import json
import os 


		#NOTE: this code might be useful for proper pheotype checking,
		#since it stores hierarchical information relating to HPO terms
		# hpo = obonet.read_obo(HPO_OBO_URL)
		# assert nx.is_directed_acyclic_graph(hpo)
		# print('Number of HPO terms: {}'.format(len(hpo)))

class VariantGraph(nx.Graph):

	#TODO: cross-reference variants with other sources
	def add_nodes_from_civic_by_gene(self, geneId, ignore_submitted = False):
		# print('Number of variants for VHL: {}'.format(len(gene.variants)))

		for variant in civic.get_gene_by_id(geneId).variants:
			self.add_node(variant.id)
			variant_node = self.nodes[variant.id]

			variant_node["variantName"] = variant.name

			#finding how many evidence items exist for the variant.
			#TODO: actually verify that evidence: supports, is germline, and is case study

			# list to hold variants' accepted evidences by id
			variant_node["evidenceAccepted"] = []

			# list to hold variants' non-accepted evidences by id
			variant_node["evidenceSubmitted"] = []

			# list to hold variants' associated phenotypes,
			#	obtained from evidence items
			variant_node["associatedPhenotypes"] = []

			for evidence in variant.evidence_items:
				# get all phenotypes from evidences
				phenotypes = [phenotype.hpo_id for phenotype in evidence.phenotypes]

				# add evidence to appropriate list, and optionally filter out evidence items
				#	that haven't been accepted
				if evidence.status == "submitted":
					variant_node["evidenceSubmitted"].append(evidence.id)
					if not ignore_submitted:
						variant_node["associatedPhenotypes"].extend(phenotypes)


				elif evidence.status == "accepted":
					variant_node["evidenceAccepted"].append(evidence.id)	
					variant_node["associatedPhenotypes"].extend(phenotypes)

			#finding the types of the variant and adding it to the node
			variant_node["variantTypes"] = []
			for variant_type in variant.variant_types:
				variant_node["variantTypes"].append(variant_type.name)

	def save_to_json_file(self, directory=""):
		VG_json = nx.readwrite.json_graph.node_link_data(self)

		with open(os.path.join(directory, "processed_variant_"+constants.GRAPH_FILENAME+ ".json"), "w") as file:
			json.dump(VG_json, file)
