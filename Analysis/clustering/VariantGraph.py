import logging

from . import ClusteringGraph

VARIANT_TEMPLATE = {
    # attributes here should be common to ALL data sources
    # list(string): holds variants' associated phenotypes obtained from article
    "KimStudents2019": {
        "PMID": "",
        "cDNA_Position": "",
        "Multiple Mutants in Case": "",
        "Mutation Event c.DNA.": "",
        "Predicted Consequence Protein Change": "",
        "variant_name": "",
        "Mutation Type": "",
        "Kindred Case": "",
        "Familial/Non-familial": "",
        "Phenotype": "",
        "Reference": "",
        "Age": "",
        "Notes": "",
        "Evidence Statements": "",
    },
    "all": {
        "cdna_start": "",
        "cdna_change": "",
        "protein_change": [],
        "variant_name": "",
        "variant_type": "",
        "variant_transition_transversion": "",
        "patient_id": []
    }
}


class VariantGraph(ClusteringGraph.ClusteringGraph):
    """Represents a graph, with variant-specific added functionality

    Attributes:
        All attributes in the nx.Graph class
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "VariantGraph"
        self.logger = logging.getLogger(self.name)
