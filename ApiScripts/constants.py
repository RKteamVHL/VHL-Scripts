from . import GraphSimilarity

# 58 corresponds to the VHL gene 
VHL_GENE_ID =  58


# url for human-phenotype-ontology
HPO_OBO_URL = "http://purl.obolibrary.org/obo/hp.obo"

# default graph cache 
GRAPH_FILENAME = "graph"

# keys to save in the networkx structure
VARIANT_GRAPH_KEYS = [
	"evidenceAccepted",
	"evidenceSubmitted",
	"associatedPhenotypes",
	"variantName",
	"variantTypes"

]

SIMILARITY_METRICS = [
	{
		"function": GraphSimilarity.score_iou,
		"kwargs": {
			"attr_name": "associatedPhenotypes"
		}
	},
	{
		"function": GraphSimilarity.score_iou,
		"kwargs": {
			"attr_name": "variantTypes"
		}
	}
]

