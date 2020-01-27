import networkx as nx
import numpy as np
import matplotlib.pyplot as plt



class GraphPlotter:
	def __init__(self, *args, **kwargs):
		self.graph = None
		if "graph" in kwargs:
			self.set_graph(kwargs["graph"])

	def set_graph(self, graph):
		assert(nx.Graph in type(graph).mro())
		self.graph = graph

	def count_nodes_by(self, field):

		data_sums = {}

		for node, n_data in self.graph.nodes(data=True):
			all_node_data =  n_data['all']
			if field in all_node_data:
				if isinstance(all_node_data[field], list):
					for v in all_node_data[field]:
						data_sums[v] = data_sums.get(v, 0) + 1
				elif isinstance(all_node_data[field], str):
					data_sums[all_node_data[field]] = data_sums.get(all_node_data[field], 0) + 1

		_values = []
		_labels = []
		for k, v in data_sums.items():
			# filter values
			if v > 3:
				_values.append(v)
				_labels.append(k)


		i_sorted = np.argsort(_values, axis=0)

		sorted_values = np.take_along_axis(np.array(_values), i_sorted, 0)
		sorted_labels = np.take_along_axis(np.array(_labels), i_sorted, 0)


		# set width of bar
		barWidth = 1



		# Set position of bar on X axis
		r1 = np.arange(len(sorted_values))
		# r2 = [x + barWidth for x in r1]
		# r3 = [x + barWidth for x in r2]

		# Make the plot
		plt.bar(r1, sorted_values, color='#7f6d5f', width=barWidth, edgecolor='white', label='var1')
		# plt.bar(r2, bars2, color='#557f2d', width=barWidth, edgecolor='white', label='var2')
		# plt.bar(r3, bars3, color='#2d7f5e', width=barWidth, edgecolor='white', label='var3')

		# Add xticks on the middle of the group bars
		plt.xlabel(field, fontweight='bold')
		plt.xticks(r1, sorted_labels, rotation='vertical')

		plt.ylabel("# of variants", fontweight='bold')

		# Create legend & Show graphic
		plt.legend()
		plt.show()

# TODO: Implement all previous visualizations
	# 1) phenotype vs # variants
	# 2) phenotype vs phenotype correlation matrix
	# 3) phenotype + mutation type vs # of variants
	# 4) phenotype vs mutation type correlation matrix
	# 5) variant location histogram
	# 6) phenotype + exon location vs # of variants
	# 7) phenotype vs exon location correlation matrix
	# 8) phenotype + domain vs # of variants

	# TODO: Implement novel visualizations
	# 1) phenotype + cluster vs # of variants
	# 2) general phenotype + cluster vs # of variants
	# 3) variantType + cluster vs # of variants
	# 4)