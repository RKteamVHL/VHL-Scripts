from civicpy import civic
from . import constants
from . import similarity_functions as sf
from . import variant_functions as vf
import networkx as nx
import json
import math
import os 


#NOTE: constants for variant analysis are declared here, but it may 
#make more sense for these to be static attributes in VarianGraph

# keys to save in the networkx structure
VARIANT_GRAPH_KEYS = [
	"evidenceAccepted",
	"evidenceSubmitted",
	"associatedPhenotypes",
	# civic(string): the civic name of the variant, 
	"variantName",
	# civic([string]): list of variant types, stored as sequence ontology names 
	"variantTypes", 
	# civic([string]): list of protein/chromosomal/mrna transcript ids for the variant
	"hgvsExpressions", 
	# civic(int): id used by civic to reference this variant
	"civicId" 

]

# list of dicts, where each dict references a similarity function, and stores
# keyword args for that function
SIMILARITY_METRICS = [
	{
		"function": sf.score_iou,
		"kwargs": {
			"attr_name": "associatedPhenotypes"
		}
	},
	{
		"function": sf.score_iou,
		"kwargs": {
			"attr_name": "variantTypes"
		}
	},
		{
		"function": sf.variant_score_domains,
		"kwargs": {
			"attr_name": "alpha"
		}
	},
		{
		"function": sf.variant_score_domains,
		"kwargs": {
			"attr_name": "beta"
		}
	}
]

# list of dicts, where each dict references a node metric function, and stores
# keyword args for that function
NODE_METRICS = [
	{
		"function": vf.affected_domains,
		"kwargs":{

		}	
	}
]

DEFAULT_DRAWING_OPTIONS = {
	'node_color': 'blue',
	'node_size': 50,
	'line_color': 'grey',
	'font_size': 6
}


class VariantGraph(nx.Graph):
	"""Represents a graph, with variant-specific added functionality

	Attributes:
		All attributes in the nx.Graph class 
	"""

	#TODO: cross-reference variants with other sources
	#TODO: potentially move these attribute-setting fuctions to node_functions.py
	def add_nodes_from_civic_by_gene(self, geneId, ignore_submitted = False):
		"""Fetches all variants for a given gene and adds them to self nodes

		Arguments:
			geneId (int): Civic id of the gene 
			ignore_submitted (bool): Determines whether variants with only
				submitted, not approved, evidence statements will be included
		"""
		to_remove = []
		for variant in civic.get_gene_by_id(geneId).variants:
			self.add_node(variant.id)
			variant_node = self.nodes[variant.id]

			variant_node["civicId"] = variant.id
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

			has_accepted_evidence = len(variant_node["evidenceAccepted"])>0

			#finding the types of the variant and adding it to the node
			variant_node["variantTypes"] = []
			for variant_type in variant.variant_types:
				variant_node["variantTypes"].append(variant_type.name)

			#finding the HGVS expressions
			variant_node["hgvsExpressions"] = []
			for hgvs in variant.hgvs_expressions:
				variant_node["hgvsExpressions"].append(hgvs)

			if ignore_submitted and not has_accepted_evidence:
				to_remove.append(variant_node["civicId"])

		self.remove_nodes_from(to_remove)


	def save_to_json_file(self, filename):
		"""Saves the graph to a file
		Arguments:
			filename (string): name of the saved file
		"""
		VG_json = nx.readwrite.json_graph.node_link_data(self)

		with open(filename, "w") as file:
			json.dump(VG_json, file)

		print("# of variants saved to file: {}".format(len(self)))

	def load_from_json_file(self, filename):
		"""Loads the graph from a file
		Arguments:
			filename (string): name of the saved file
		"""
		with open(filename, "r") as file:
			VG_json = json.load(file)
			VG = nx.readwrite.json_graph.node_link_graph(VG_json)
			self.add_nodes_from(VG.nodes(data=True))	

		print("# of variants loaded from file: {}".format(len(self)))

	#TODO: it may make sense to have the metric functions and arguments
	# inputted as arguments to this method
	def calculate_node_attributes(self):
		"""Calculates node metrics for each node, based on NODE_METRICS
		"""
		VG_nodes = list(self.nodes(data=True))
		for i in range(0,len(VG_nodes)):
			for metric in NODE_METRICS:
				(name, score) = metric["function"](VG_nodes[i][1], **metric["kwargs"])
				nId = VG_nodes[i][0]
				self.nodes[nId][name] = score

	#TODO: it may make sense to have the metric functions and arguments
	# inputted as arguments to this method
	def calculate_similarities(self):
		"""Calculates similarity metrics between nodes from SIMILARITY_METRICS
		"""

		VG_nodes = list(self.nodes(data=True))
		#assuming undirected graphs with no self connections
		for i in range(0,len(VG_nodes)):
			for j in range(i+1, len(VG_nodes)):

				similarities = {}
				for metric in SIMILARITY_METRICS:
					(name, score) = metric["function"](VG_nodes[i][1], VG_nodes[j][1], **metric["kwargs"])
					similarities[name] = score

				has_value = [ not math.isclose(v, 0, rel_tol=1e-5) for v in similarities.values()]
				if any(has_value):
					self.add_edge(VG_nodes[i][0], VG_nodes[j][0], **similarities)

	#TODO: this function may not be needed
	def draw_graph(self, **options):
		draw_options = DEFAULT_DRAWING_OPTIONS
		for key, val in options.items():
				draw_options[key] = val 




		nx.draw(self, **draw_options)

		plt.show()

		# filt_edges = [(u,v,d) for (u,v,d) in variant_graph.edges(data=True) if d["weight"] !=0] 
		# edge_widths = [ d["weight"]/10 for (u,v,d) in filt_edges ]	
		# edge_colors = [[d["phenotype_weight"], d["vartype_weight"],0, 0.5] for (u,v,d) in filt_edges]



