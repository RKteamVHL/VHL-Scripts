import logging

from . import ClusteringGraph

PATIENT_TEMPLATE = {
    "KimStudents2019": {
        # from raw data
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
        "Notes": "",
        "Evidence Statements": "",
        # processed from raw data
        "Age": "",
        "Gender": "",
    },
    "all": {
        "age": None,
        "gender": "NA",
        "associated_phenotypes": [],
        "family_id": [],
        "variant_id": [],
        "evidence_statement": "",
        "pmid": []
    }

}


class PatientGraph(ClusteringGraph.ClusteringGraph):
    """Represents a graph, with patient-specific added functionality

    Attributes:
        All attributes in the nx.Graph class
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = "PatientGraph"
        self.logger = logging.getLogger(self.name)
