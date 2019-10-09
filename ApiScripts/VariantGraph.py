from civicpy import civic
from . import similarity_functions as sf
from . import variant_functions as vf
from sklearn.cluster import SpectralClustering
from snf import compute
from snf import metrics

import copy
import networkx as nx
import json
import math
import os 

# url for human-phenotype-ontology
HPO_OBO_URL = "http://purl.obolibrary.org/obo/hp.obo"


# NOTE: constants for variant analysis are declared here, but it may 
# make more sense for these to be static attributes in VariantGraph

# variant node attributes BEFORE any variant functions are done on a given variant
VARIANT_TEMPLATE = {
	"civic": {
		# list(int): the civic evidence ids that support this variant and are accepted
		"evidenceAccepted": [],

		# list(int): the civic evidence ids that support this variant and are not yet accepted
		"evidenceSubmitted": [],

		# list(string): holds variants' associated phenotypes obtained from evidence items
		"associatedPhenotypes": [],

		# string: the civic name of the variant, 
		"variantName": None,

		# list(string): list of variant types, stored as sequence ontology names 
		"variantTypes": [], 

		# list(string): list of protein/chromosomal/mrna transcript ids for the variant
		"hgvsExpressions": [], 

		# int: id used by civic to reference this variant
		"id": None,

		# string: the variant cdna change
		"cdnaChange": None,

		# string: the variant aa change
		"proteinChange": None

	},

	"gnomad": {
		# string: dbSNP id for the variant
		"rsID": None,

		# list(string): list of variant types, stored as sequence ontology names 
		"variantTypes": [],

		# string: the variant cdna change
		"cdnaChange": None,

		# string: the variant aa change
		"proteinChange": None,

		# int: the variant allele count
		"alleleCount": None,

		# float: the variant allele frequency
		"alleleFrequency": None

	}

	"students2019": {
		# int: pmid of the article where variant is affirmed 
		"pmid": None,

		# list(string): holds variants' associated phenotypes obtained from article 
		"associatedPhenotypes": [],

		# list(string): list of variant types, stored as sequence ontology names 
		"variantTypes": [],

		# string: the variant cdna change
		"cdnaChange": None,

		# string: the variant aa change
		"proteinChange": None,

	}

}

# list of dicts, where each dict references a similarity function, and stores
# keyword args for that function
SIMILARITY_METRICS = {
	"score_iou_associatedPhenotypes": {
		"function": sf.score_iou,
		"kwargs": {
			"attr_name": "associatedPhenotypes"
		}
	},
	"score_iou_variantTypes": {
		"function": sf.score_iou,
		"kwargs": {
			"attr_name": "variantTypes"
		}
	},
	"score_domain_alpha": {
		"function": sf.variant_score_domains,
		"kwargs": {
			"domain": "alpha"
		}
	},
	"score_domain_beta": {
		"function": sf.variant_score_domains,
		"kwargs": {
			"domain": "beta"
		}
	}
}

# list of dicts, where each dict references a node metric function, and stores
# keyword args for that function
NODE_METRICS = {
	"affected_domains": {
		"function": vf.affected_domains,
		"kwargs":{

		}	
	}
}

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
			node_id = len(self.nodes())

			# initializing the node with empty civic template. It's important here
			# to only include the civic dict, since it makes 
			# merging node data easier later
			node_template = {}
			node_template["civic"] = copy.deepcopy(VARIANT_TEMPLATE["civic"])

			self.add_node(node_id, attr_dict=node_template)
			variant_node = self.nodes[node_id]

			variant_node["civicId"] = variant.id
			variant_node["variantName"] = variant.name

			#finding how many evidence items exist for the variant.
			#TODO: actually verify that evidence: supports, is germline, and is case study

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
			self.add_edges_from(VG.edges(data=True))

		print("# of variants loaded from file: {}".format(len(self)))

	#NOTE: it may make sense to have the metric functions and arguments
	# inputted as arguments to this method
	def calculate_node_attributes(self):
		"""Calculates node metrics for each node, based on NODE_METRICS
		"""
		VG_nodes = list(self.nodes(data=True))
		for i in range(0,len(VG_nodes)):
			for name, metric in NODE_METRICS.items():
				score = metric["function"](VG_nodes[i][1], **metric["kwargs"])
				nId = VG_nodes[i][0]
				self.nodes[nId][name] = score

	#NOTE: it may make sense to have the metric functions and arguments
	# inputted as arguments to this method
	def calculate_similarities(self):
		"""Calculates similarity metrics between nodes from SIMILARITY_METRICS
		"""

		VG_nodes = list(self.nodes(data=True))
		#assuming undirected graphs with no self connections
		for i in range(0,len(VG_nodes)):
			for j in range(i, len(VG_nodes)):

				similarities = {}
				for name, metric in SIMILARITY_METRICS.items():
					score = metric["function"](VG_nodes[i][1], VG_nodes[j][1], **metric["kwargs"])
					similarities[name] = score

				has_value = [ not math.isclose(v, 0, rel_tol=1e-5) for v in similarities.values()]
				if any(has_value):
					self.add_edge(VG_nodes[i][0], VG_nodes[j][0], **similarities)

	def calculate_snf(self):

		adjmats = self.get_adjacency_mats(dense=True)

		#running SNF
		fused_ndarray = compute.snf(adjmats)
		clust_count1, clust_count2 = compute.get_n_clusters(fused_ndarray)


		sc = SpectralClustering(clust_count1, affinity='precomputed', n_init=100, assign_labels='discretize')
		sc.fit(fused_ndarray)

		fused_labels = sc.labels_

		silhouette = metrics.silhouette_score(fused_ndarray, fused_labels)

		print("Cluster estimates: {}, {}".format(clust_count1, clust_count2))
		print("Labels: ", fused_labels)


		# Merging the fused edge weights back into the original graph
		fused_graph = nx.convert_matrix.from_numpy_array(fused_ndarray)

		self.add_edges_from(fused_graph.edges(data=True))

		for i in range(0, len(self.nodes())):
			node = self.nodes[i]
			node['spectral_label'] = fused_labels[i].item()




	# NOTE: this code probably isn't needed- remove down the line			
	# def remove_isolates(self):
	# 	filtered_graph = self.copy()
	# 	filtered_graph.remove_nodes_from(list(nx.isolates(filtered_graph)))

	# 	relabelled_graph = nx.convert_node_labels_to_integers(filtered_graph, label_attribute="prev_id")
	# 	self.clear()
	# 	self.add_nodes_from(relabelled_graph.nodes(data=True))
	# 	self.add_edges_from(relabelled_graph.edges(data=True))
		


	def get_adjacency_mats(self, dense=False):
		mats = []
		nodes = self.nodes()
		for metric in SIMILARITY_METRICS.keys():
			mat = nx.adjacency_matrix(self, weight=metric)
			if dense:
				mats.append(mat.todense())
			else:
				mats.append(mat)

		return mats




