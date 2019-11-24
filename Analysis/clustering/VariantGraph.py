from .. import similarity_functions as sf
from .. import variant_functions as vf
from ..fetching.Fetcher import FETCHING_FACTORY
from sklearn.cluster import SpectralClustering
from snf import compute
from snf import metrics

import copy
import networkx as nx
import numpy as np
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
		"rsid": None,

		# int: the variant allele count
		"alleleCount": None,

		# float: the variant allele frequency
		"alleleFrequency": None

	},

	"KimStudents2019": {
		# int: pmid of the article where variant is affirmed 
		"PMID": None	
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

		# string: the variant cdna change. this may not be needed, since
		#cdna change is encoded into node variant ids
		"cdnaChange": None,

		# # string: the variant aa change
		# "proteinChange": None,

	}
}

# list of dicts, where each dict references a similarity function, and stores
# keyword args for that function
SIMILARITY_METRICS = {
	"score_hpo": {
		"function": sf.variant_hpo_distance,
		"kwargs":{

		}
	},

	"score_so": {
		"function": sf.variant_so_distance,
		"kwargs":{

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

SCORE_THRESHOLD = 1e-2


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
			if isinstance(node_i, str):
				all_node['cdnaChange'] = node_i 
			self.nodes[node_i]['all'] = all_node


	def save_to_json_file(self, filename, nodes_only=False):
		"""Saves the graph to a file
		Arguments:
			filename (string): name of the saved file
		"""
		G = self
		if nodes_only:
			G = nx.Graph()
			G.add_nodes_from(self.nodes(data=True))

		VG_json = nx.readwrite.json_graph.node_link_data(G)

		with open(filename, "w") as file:
			json.dump(VG_json, file)

		print(f"# of variants saved to {filename}: {len(self)}")

	def load_from_json_file(self, filename, nodes_only=False):
		"""Loads the graph from a file
		Arguments:
			filename (string): name of the saved file
		"""
		with open(filename, "r") as file:
			VG_json = json.load(file)
			VG = nx.readwrite.json_graph.node_link_graph(VG_json)
			self.add_nodes_from(VG.nodes(data=True))	
			if not nodes_only:
				self.add_edges_from(VG.edges(data=True))

		print(f"# of variants loaded from {filename}: {len(self)}")

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

		This method also now calculates SNF metrics
		"""

		VG_nodes = list(self.nodes(data=True))
		#assuming undirected graphs with no self connections
		for i in range(0,len(VG_nodes)):
			for j in range(i+1, len(VG_nodes)):

				similarities = {}
				for name, metric in SIMILARITY_METRICS.items():
					score = metric['function'](VG_nodes[i][1], VG_nodes[j][1], **metric['kwargs'])
					similarities[name] = score

				has_value = [ not math.isclose(v, 0, abs_tol=SCORE_THRESHOLD) for v in similarities.values()]
				if any(has_value):
					self.add_edge(VG_nodes[i][0], VG_nodes[j][0], **similarities)

		adjmats = self.get_adjacency_mats(types=SIMILARITY_METRICS.keys())

		node_list = list(self.nodes())

		#running SNF
		fused_ndarray = compute.snf(adjmats)
		
		it = np.nditer(fused_ndarray, flags=['multi_index'])
		while not it.finished:
			# if the similarity is not 0
			if not math.isclose(it[0].item(), 0, abs_tol=SCORE_THRESHOLD):
				n1 = node_list[it.multi_index[0]]
				n2 = node_list[it.multi_index[1]]
				self.add_edge(n1, n2, fused_similarity=it[0].item())
			it.iternext()

	def cluster_by(self, affinity_type, num_clusters):

		adj_mat = self.get_adjacency_mats(types=[affinity_type])


		sc = SpectralClustering(num_clusters, affinity='precomputed', n_init=100, assign_labels='discretize')		
		sc.fit(adj_mat[0])

		fused_labels = sc.labels_

		silhouette = metrics.silhouette_score(adj_mat[0], fused_labels)

		print("Labels: ", fused_labels)

		node_list = list(self.nodes())
		it = np.nditer(fused_labels, flags=['multi_index'])
		while not it.finished:
			feature_label = f'{affinity_type}_label'
			n1 = node_list[it.multi_index[0]]
			self.nodes[n1][feature_label] = it[0].item()
			it.iternext()




	def get_adjacency_mats(self, types):
		assert not isinstance(types, str)
		mats = []
		nodes = self.nodes()
		for metric in types:
			mat = nx.to_numpy_array(self, weight=metric)
			mats.append(mat)

		return mats




