import networkx as nx
import os
import obonet
import urllib
import logging
from . import config

# sources: https://www.ncbi.nlm.nih.gov/protein/4507891, https://www.uniprot.org/uniprot/P40337,
# https://www.ebi.ac.uk/interpro/protein/UniProt/P40337/
VHL_FUNCTIONAL_REGIONS = {
    "tumour_suppression": set(range(64, 204)),
    "⍺-Domain": set(range(156, 205)),  # PF17211
    "β-Domain": set(range(63, 144)),  # PF01847
    "GXEEX8": set(range(14, 54)),
    "HIF1_alpha_binding": {67, 69, 75, 77, 78, 79, 88, 91, 98, 99, 105, 106, 107, 108, 109, 110, 111, 112, 115, 117},
    "ElonginB_ElonginC_binding": {79, 153, 159, 161, 162, 163, 165, 166, 174, 177, 178, 184}
}
VHL_FUNCTIONAL_REGIONS['Outside of ⍺-Domain and β-Domain'] = set(range(1, 214)) - (
            VHL_FUNCTIONAL_REGIONS['⍺-Domain'] | VHL_FUNCTIONAL_REGIONS['β-Domain'])

GENERAL_HPO_TERMS = [
    'neuroendocrine neoplasm',
    'renal cell carcinoma',
    'hemangioblastoma',
    'cerebellar hemangioblastoma',
    'spinal hemangioblastoma',
    'retinal capillary hemangioma',
    'pancreatic endocrine tumor',
    'abnormality of the kidney',
    'abnormality of the pancreas',
    # 'abnormality of the epididymis',
    # 'abnormality of the ovary',
    # 'endolymphatic sac tumor'
]

GENERAL_SO_TERMS = [
    'deletion',
    'exon_loss_variant',
    'missense_variant',
    'stop_gained',
    'utr_variant',
    'inframe_indel',
    'delins',
    'frameshift_variant',
    'splice_site_variant',
    'start_lost',
    'synonymous_variant',
    'intron_variant',
    'stop_lost'
]

HPO_ABBREVIATIONS = {
    'hemangioblastoma': 'HB',
    'cerebellar hemangioblastoma': 'CHB',
    'spinal hemangioblastoma': 'SHB',
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

SO_TERM_TYPES = {
    'truncating': ['stop_gained', 'deletion', 'exon_loss_variant', 'start_lost', 'frameshift_variant'],
    'non-truncating': ['missense_variant', 'inframe_indel'],
}

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
    d['name_spaceless'] = d['name'].replace(' ', '')
# make a copy of merged network- node keys: term
_h = nx.relabel_nodes(_g, {n: d['name'] for n, d in _g.nodes(data=True)})
_h2 = nx.relabel_nodes(_g, {n: d['name_spaceless'] for n, d in _g.nodes(data=True)})
# combine together OBO network so that node keys: id or term
_f1 = nx.compose(_g, _h)
_f = nx.compose(_f1, _h2)
# casefold and strip the keys, so the result is an agglomeration of nodes:
# SO(term) + SO(id) + HPO(term) + HPO(id) + HPO(no spaces), all lowercase and stripped

OBONET = nx.DiGraph(nx.relabel_nodes(_f, {n: n.casefold().strip() for n in _f.nodes()}))

def get_valid_obo(term_or_id, return_as=None):
    """
    Take either a term (e.g., HP:0010797) or name (e.g., Hemangioblastoma) and return
    a matching HPO term or name to verify it exists
    @param term_or_id: The term to match
    @param return_as: can be 'id', 'name', or 'name_spaceless'. returns the attribute, rather than node
    @return:
    """
    if isinstance(term_or_id, str):
        tid = term_or_id.strip().casefold()
        if tid in OBONET:
            to_return = OBONET.nodes[tid]
            if return_as == 'id':
                to_return = to_return['id']
            elif return_as == 'name':
                to_return = to_return['name']
            elif return_as == 'name_spaceless':
                to_return = to_return['name_spaceless']
            return to_return
        else:
            logging.warning(f"Could not find an OBO node for {tid}")
    else:
        logging.warning(f"{term_or_id} is not a string")


# TODO: this will break if 'return_as' is ever added as a command-line argument. error-checking may also be redundant
GENERAL_HPO_NODES = [get_valid_obo(term, return_as=config.OBO_RETURN_TYPE) for term in GENERAL_HPO_TERMS]
GENERAL_SO_NODES = [get_valid_obo(term, return_as=config.OBO_RETURN_TYPE) for term in GENERAL_SO_TERMS]


def generalized_obo(term, general_nodes, return_as=None, use_abbreviation=False):
    general_obo = None

    valid_obo = get_valid_obo(term, return_as=return_as)
    # we do not need to log if valid_obo is None, since get_valid_obo already does that
    if valid_obo is not None:
        for successor in nx.bfs_successors(OBONET, valid_obo):
            try:
                i = general_nodes.index(successor[0])
                general_obo = general_nodes[i]
                break
            except ValueError as e:
                pass

    if general_obo is None:
        logging.warning(f"Could not find a generalized term for {valid_obo}")

    if use_abbreviation:
        general_obo = HPO_ABBREVIATIONS[general_obo]

    return general_obo


def generalized_so_terms(so_type, **kwargs):
    """Given a node, find its general so_type
    """
    return generalized_obo(so_type, GENERAL_SO_NODES, **kwargs)


def generalized_vhl_phenotype(phenotype, **kwargs):
    """Given a node, find its general disease type
    """
    return generalized_obo(phenotype, GENERAL_HPO_NODES, **kwargs)
