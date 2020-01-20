import networkx as nx

class GraphPlotter:
	def __init__(self):
		self.graph = None

	def set_graph(self, graph):
		assert(nx.Graph in type(graph).mro())
		self.graph = graph

	def summarize_nodes_by(self, field):
		pass
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