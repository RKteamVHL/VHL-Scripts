import networkx as nx
import numpy as np
import logging
import json
import math
import csv

from sklearn.cluster import SpectralClustering
from snf import compute
from snf import metrics

DELIMITER = "\t"
SCORE_THRESHOLD = 1e-2


class ClusteringGraph(nx.Graph):
	"""Represents a graph, with clustering, saving, and loading capabilities

	Attributes:
		All attributes in the nx.Graph class
	"""
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.name = "ClusteringGraph"
		self.logger = logging.getLogger(self.name)

	def save_to_json_file(self, filename, nodes_only=False):
		"""Saves the graph to a file
		Arguments:
			filename (string): name of the saved file
			nodes_only (bool): set to true if only nodes should be saved
		"""
		g = self
		if nodes_only:
			g = nx.Graph(**self.graph)
			g.add_nodes_from(self.nodes(data=True))

		g_json = nx.readwrite.json_graph.node_link_data(g)

		with open(filename, "w") as file:
			json.dump(g_json, file)

		print(f"# of nodes saved to {filename}: {len(self)}")

	def save_to_text_file(self, filename, headers):
		"""Saves the graph to a file
		Arguments:
			filename (string): name of the saved file
			headers (list(string)): list of string headers for text file
		"""

		with open(filename, "w", encoding="utf8", newline='\n') as file:
			writer = csv.DictWriter(file, delimiter=DELIMITER, fieldnames=headers)
			writer.writeheader()

			for node, data in self.nodes(data=True):
				row_dict = {}
				for k, v in data['all'].items():
					if isinstance(v, list):
						row_dict[k] = ";".join(v)
					else:
						row_dict[k] = v

				writer.writerow(row_dict)

		print(f"# of nodes saved to {filename}: {len(self)}")

	def load_from_json_file(self, filename, nodes_only=False):
		"""Loads the graph from a file
		Arguments:
			filename (string): name of the saved file
			nodes_only (bool): set to true if only nodes should be loaded
		"""
		with open(filename, "r") as file:
			g_json = json.load(file)
			g = nx.readwrite.json_graph.node_link_graph(g_json)
			self.add_nodes_from(g.nodes(data=True))
			if not nodes_only:
				self.add_edges_from(g.edges(data=True))

		print(f"# of nodes loaded from {filename}: {len(self)}")

	def calculate_node_attributes(self, node_metrics):
		"""Calculates metrics for each node in this graph

		:param node_metrics: a dictionary, where each key is the name of a metric, and its corresponding value is
		another dictionary with the following structure:

			node_metrics = {
				"metric1": {
					"function": function1,
					"kwargs":{
						"arg1": 1
					}
				},
				"metric2": {
					"function": function2,
					"kwargs":{
						"arg1": 10
					}
				}
			}

		The input to each of the functions is the node's data dictionary
		:return:
		"""
		g_nodes = list(self.nodes(data=True))
		for i in range(0, len(g_nodes)):
			for name, metric in node_metrics.items():
				try:
					score = metric['function'](g_nodes[i][1], **metric['kwargs'])
					n_id = g_nodes[i][0]
					self.nodes[n_id]['all'][name] = score
				except ValueError as e:
					self.logger.warning(repr(e))

	def calculate_similarities(self, similarity_metrics):
		"""Calculates similarity metrics for every pair of nodes in this graph

		:param similarity_metrics: a dictionary, where each key is the name of a metric, and its corresponding value is
		another dictionary with the following structure:

			similarity_metrics = {
				"metric1": {
					"function": function1,
					"kwargs":{
						"arg1": 1
					}
				},
				"metric2": {
					"function": function2,
					"kwargs":{
						"arg1": 10
					}
				}
			}

		The input to each of the functions is both of the nodes' data dictionary
		:return:
		"""

		g_nodes = list(self.nodes(data=True))
		# assuming undirected graphs with no self connections
		for i in range(0, len(g_nodes)):
			for j in range(i + 1, len(g_nodes)):

				similarities = {}
				for name, metric in similarity_metrics.items():
					score = metric['function'](g_nodes[i][1], g_nodes[j][1], **metric['kwargs'])
					similarities[name] = score

				has_value = [not math.isclose(v, 0, abs_tol=SCORE_THRESHOLD) for v in similarities.values()]
				if any(has_value):
					self.add_edge(g_nodes[i][0], g_nodes[j][0], **similarities)

		adjmats = self.get_adjacency_mats(types=similarity_metrics.keys())

		node_list = list(self.nodes())

		# running SNF
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

		# get the desired affinity matrix
		adj_mat = self.get_adjacency_mats(types=[affinity_type])

		# cluster based on given adjmat and # clusters
		sc = SpectralClustering(num_clusters, affinity='precomputed', n_init=1000)
		sc.fit(adj_mat[0])

		fused_labels = sc.labels_

		silhouette = metrics.silhouette_score(adj_mat[0], fused_labels)

		print("Labels: ", fused_labels)

		feature_label = f'{affinity_type}_label'

		# assign each node its appropriate label
		node_list = list(self.nodes())
		it = np.nditer(fused_labels, flags=['multi_index'])
		while not it.finished:
			label_i = it[0].item()
			n1 = node_list[it.multi_index[0]]
			self.nodes[n1]['all'][feature_label] = label_i

			it.iternext()

		self.sort_by_attribute(feature_label)

	def get_adjacency_mats(self, types):
		assert not isinstance(types, str)
		mats = []
		for metric in types:
			mat = nx.to_numpy_array(self, weight=metric)
			mats.append(mat)

		return mats

	def sort_by_attribute(self, feature_label):

		sorted_nodes = sorted(list(self.nodes(data=True)), key=lambda n: n[1]['all'][feature_label])
		# for n in sorted_nodes:
		# 	print(n[1][feature_label])
		g = nx.Graph()

		g.add_nodes_from(sorted_nodes)
		g.add_edges_from(self.edges(data=True))

		self.clear()
		self.add_nodes_from(g.nodes(data=True))
		self.add_edges_from(g.edges(data=True))
