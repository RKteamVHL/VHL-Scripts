from .. import similarity_functions as sf
from .. import variant_functions as vf
from ..fetching.Fetcher import FETCHING_FACTORY
from sklearn.cluster import SpectralClustering
from snf import compute
from snf import metrics
from openpyxl import load_workbook

import copy
import networkx as nx
import csv
import json
import math
import os 
import re

# NOTE: constants for variant analysis are declared here, but it may 
# make more sense for these to be static attributes in VariantGraph

# variant node attributes BEFORE any variant functions are done on a given variant
# All attributes the the 'all' template also apply to individual sources
VARIANT_TEMPLATE = {
	"Civic": {
		# list(int): the civic evidence ids that support this variant and are accepted
		"evidenceAccepted": [],

		# list(int): the civic evidence ids that support this variant and are not yet accepted
		"evidenceSubmitted": [],

		# list(string): holds variants' phenotypes obtained from submitted evidence items
		"phenotypesSubmitted": [],

		# list(string): holds variants' phenotypes obtained from accepted evidence items
		"phenotypesAccepted": [],

		# string: the civic name of the variant, 
		"variantName": None,		

		# list(string): list of protein/chromosomal/mrna transcript ids for the variant
		"hgvsExpressions": [], 

		# int: id used by civic to reference this variant
		"civicId": None

	},

	"Gnomad": {
		# string: dbSNP id for the variant
		"rsID": None,

		# int: the variant allele count
		"alleleCount": None,

		# float: the variant allele frequency
		"alleleFrequency": None

	},

	"KimStudents2019": {
		# int: pmid of the article where variant is affirmed 
		"pmid": None	
	},

	"ClinVar": {
		# string: dbSNP id for the variant
		"RS# (dbSNP)": None
	},

	# attributes here should be common to ALL data sources
	"all":	{
		# list(string): holds variants' associated phenotypes obtained from article 
		"associatedPhenotypes": [],

		# list(string): list of variant types, stored as sequence ontology names 
		"variantTypes": [],

		# string: the variant cdna change
		"cdnaChange": None,

		# # string: the variant aa change
		# "proteinChange": None,

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

	def add_nodes_from_db(self, db):
		fetcher = FETCHING_FACTORY[db]()
		fetcher.process()
		rows = fetcher.rows

		for row in rows:
			node_id = len(self.nodes())

			# sets the next node id based on cdna change
			#this is used to automatically merge nodes across databases
			cdna_id = row['cdnaChange']

			valid_id = cdna_id if cdna_id else node_id
			self.add_node(valid_id, **{db: copy.deepcopy(VARIANT_TEMPLATE[db])})
			variant_node = self.nodes[valid_id][db]
			for key in VARIANT_TEMPLATE[db].keys():
				variant_node[key] = row[key]
			for key in VARIANT_TEMPLATE['all'].keys():
				variant_node[key] = row[key]

	def merge_nodes(self):
		"""Reconciles the attributes of separate sources
		"""

		for node_i, node_d in self.nodes(data=True):
			all_node = copy.deepcopy(VARIANT_TEMPLATE['all'])			

			# TODO: the merging technique for every field here
			# should be scrutinized and modified if needed
			for db, db_var in node_d.items():
				all_node['associatedPhenotypes'].extend(db_var['associatedPhenotypes'])
				all_node['variantTypes'].extend(db_var['variantTypes'])

			# trick for finding unique values
			all_node['associatedPhenotypes'] = list(set(all_node['associatedPhenotypes']))
			all_node['variantTypes'] = list(set(all_node['variantTypes']))
			self.nodes[node_i]['all'] = all_node


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
				score = metric['function'](VG_nodes[i][1], **metric['kwargs'])
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
					score = metric['function'](VG_nodes[i][1], VG_nodes[j][1], **metric['kwargs'])
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

		print(len(self))
		self.add_edges_from(fused_graph.edges(data=True))



		for i in range(0, len(self)):
			node = list(self.nodes(data=True))[i][1]
			node['spectral_label'] = fused_labels[i].item()

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




