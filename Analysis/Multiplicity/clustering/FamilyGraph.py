from . import ClusteringGraph
import logging

FAMILY_TEMPLATE = {
	# attributes here should be common to ALL data sources
	# list(string): holds variants' associated phenotypes obtained from article
	"KimStudents2019": {
	},
	"all": {
		"patient_id": [],
	}
}


class FamilyGraph(ClusteringGraph.ClusteringGraph):
	"""Represents a graph, with family-specific added functionality

	Attributes:
		All attributes in the nx.Graph class 
	"""
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.name = "FamilyGraph"
		self.logger = logging.getLogger(self.name)
