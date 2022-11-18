import networkx as nx
import os
import obonet
import urllib
from . import config
DISEASE_ENTITY_TO_HPO = {
        'asymptomatic': 'asymptomatic',
        'adrenalpheochromocytoma': 'neuroendocrineneoplasm',
        'cerebellarhemangioblastoma': 'hemangioblastoma',
        'clearcellrenalcellcarcinoma': 'renalcellcarcinoma',
        'hemangioblastoma': 'hemangioblastoma',
        'hemangioma': 'hemangioma',
        'neoplasmoftheliver': 'neoplasmoftheliver',
        'neuroendocrineneoplasm': 'neuroendocrineneoplasm',
        'pancreaticcysts': 'pancreaticcysts',
        'pancreaticendocrinetumor': 'pancreaticendocrinetumor',
        'pheochromocytoma': 'neuroendocrineneoplasm',
        'renalcellcarcinoma': 'renalcellcarcinoma',
        'renalcyst': 'renalcyst',
        'renalneoplasm': 'renalneoplasm',
        'retinalcapillaryhemangioma': 'hemangioblastoma',
        'spinalhemangioblastoma': 'hemangioblastoma'
}

GENERAL_HPO_TERMS = [
    'neuroendocrine neoplasm',
    'renal cell carcinoma',
    'hemangioblastoma',
    'retinal capillary hemangioma',
    'pancreatic endocrine tumor',
    'abnormality of the kidney',
    'abnormality of the pancreas',
    # 'abnormality of the epididymis',
    # 'abnormality of the ovary',
    # 'endolymphatic sac tumor'
]

HPO_ABBREVIATIONS = {
    'hemangioblastoma': 'CHB',
    'renal cell carcinoma': 'RCC',
    'retinal capillary hemangioma': 'RA',
    'neuroendocrine neoplasm': 'PPGL',
    'pancreatic endocrine tumor': 'PNET',
    'endolymphatic sac tumor': 'ELST',
    'abnormality of the pancreas': 'PCT',
    'abnormality of the kidney': 'RCT',
    'abnormality of the epididymis': 'ECT',
    'abnormality of the ovary': 'OCT'
}
_abbrv = {v: k for k, v in HPO_ABBREVIATIONS.items()}
HPO_ABBREVIATIONS.update(_abbrv)

## Phenotype and Sequency Ontology Utilities
SO_NAME = 'SequenceOntology'
SO_FILENAME = 'so.obo'
SO_HREF = 'https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so.obo'
SO_DIR = os.path.join(config.DIRS['lib'], 'so.obo.txt')

HPO_NAME = 'HumanPhenotypeOntology'
HPO_FILENAME = 'hp.obo'
HPO_HREF = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'
HPO_DIR = os.path.join(config.DIRS['lib'], 'hp.obo.txt')

if not config.USE_CACHE:
    h_request = urllib.request.Request(SO_HREF, method="GET")
    with urllib.request.urlopen(h_request) as response:
        with open(SO_DIR, 'wb') as file:
            file.write(response.read())
    h_request = urllib.request.Request(HPO_HREF, method="GET")
    with urllib.request.urlopen(h_request) as response:
        with open(HPO_DIR, 'wb') as file:
            file.write(response.read())

# gather HPO and SO into a network- node keys: id
_g = nx.compose(obonet.read_obo(SO_HREF), obonet.read_obo(HPO_HREF))

# add the id to the node attributes
for n, d in _g.nodes(data=True):
    d['id'] = n.casefold().strip()
    d['name'] = d['name'].casefold().strip()

# make a copy of merged network- node keys: term
_h = nx.relabel_nodes(_g, {n: d['name'] for n, d in _g.nodes(data=True)})
# combine together OBO network so that node keys: id or term
_f = nx.compose(_g, _h)
# casefold and strip the keys, so the result is an agglomeration of nodes:
# SO(term) + SO(id) + HPO(term) + HPO(id), all lowercase and stripped


OBONET = nx.DiGraph(nx.relabel_nodes(_f, {n: n.casefold().strip() for n in _f.nodes()}))
OBONET_UD = OBONET.to_undirected()


def get_valid_obo(term_or_id, obo_type='name'):
    # this function will take either a term (e.g., HP:0010797) or name (e.g., Hemangioblastoma) and return
    # a matching HPO term or name to verify it exists
    tid = term_or_id.strip().casefold()
    if tid in OBONET:
        node_data = OBONET.nodes[tid][obo_type]
        return node_data
    else:
        raise ValueError(f"Could not find an OBO node for {tid}")


