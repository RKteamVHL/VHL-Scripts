import os
import re

import math
import networkx as nx
import obonet
from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.Seq import Seq

### The functions here are related to variant-based analysis for each individual node.
# Variant-Variant similarity functions are in similarity_functions.py


## Transcript CDS and Protein information
# The following code loads the VHL201 transcript, which contains sequences for:
# utr5-exon1-exon2-exon3-utr3
# introns are not included here
# TODO: make this less hardcoded, and function based on fasta descriptions

# The particular fasta file was obtained from: 
# https://useast.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000134086;r=3:10141008-10152220;t=ENST00000256474
VHL201_FASTA = os.path.join('Analysis', 'files', 'Homo_sapiens_VHL_201_sequence.fa')

vhl_seqs = [str(seqreq.seq) for seqreq in SeqIO.parse(VHL201_FASTA, "fasta")]

utr3 = Seq(vhl_seqs[5])
utr5 = Seq(vhl_seqs[6])
exon1 = Seq(vhl_seqs[0].replace(vhl_seqs[6], ''))
exon2 = Seq(vhl_seqs[1])
exon3 = Seq(vhl_seqs[2].replace(vhl_seqs[5], ''))

VHL_CDS = exon1 + exon2 + exon3

VHL_PROTEIN = VHL_CDS.translate()

## Protein Domain Information
#	resources used for finding these values:
#	https://www.ensembl.org/Homo_sapiens/Gene/Splice?db=core;g=ENSG00000134086;r=3:10141008-10152220

# for ENSG00000134086, ensembl links to Pfam ids: 
# PF01847 (VHL beta domain); http://pfam.xfam.org/family/PF01847
# PF17211 (short C-terminal alpha helical domain); http://pfam.xfam.org/family/PF17211
# summary: http://pfam.xfam.org/protein/VHL_HUMAN

# On these pfam pages, the aa ranges from Uniprot are found on the Structures tab.
# the most conservative ranges were used for the dictionary below

# Q: Gene3D has a different set of coordinates for these- which is more
# trustworthy, Gene3D or Pfam?

# Q: Although both human ref builds (GRCh38.p12, GRCh37.p13) have the same
# transcript ID for the current 213aa transcript (ENST00000256474.2),
# they have  different domains for the alpha/beta regions.
# Which one should be used? Below the most recent GRCH38.p12 are used.

CURRENT_VHL_TRANSCRIPT = {
    'ensembl': 'ENST00000256474.2',
    'ncbi': 'NM_000551.3'
}
CURRENT_VHL_PROTEIN = 'NP_000542.1'
CURRENT_VHL_GENE = 'ENSG00000134086'

PFAM_VHL_DOMAINS = {
    "ENST00000256474.2": {
        # 63-143, inclusive
        "beta": range(63, 144),
        # 149-204, inclusive
        "alpha": range(149, 205)
    },

    # Q: I couldn't find these values on pfam, only on ensembl. Did
    # ensembl do the mapping, or is it somewhere on pfam?
    "ENST00000345392.2": {
        # 63-114, inclusive
        "beta": range(63, 115),
        # 114-163, inclusive
        "alpha": range(114, 164)
    }
}

GENE3D_VHL_DOMAINS = {
    "ENST00000256474.2": {
        # 51-152, inclusive, 102bp
        "beta": range(51, 153),
        # 153-204, inclusive, 52 bp
        "alpha": range(153, 205)
    },

    "ENST00000345392.2": {
        # 51-117, inclusive
        "beta": range(51, 118),
        # 118-163, inclusive
        "alpha": range(118, 164)
    }
}
ALPHA_LEN = len(GENE3D_VHL_DOMAINS[CURRENT_VHL_TRANSCRIPT['ensembl']]['alpha'])
BETA_LEN = len(GENE3D_VHL_DOMAINS[CURRENT_VHL_TRANSCRIPT['ensembl']]['beta'])
CDS_LEN = len(VHL_PROTEIN)
## Extraction of Variant Amino Acid / Nucleotide Changes

# this regex doesn't get detailed variant info
CDNA_REGEX = re.compile('(?P<id>.*?)?(?P<gene>\(.+\))?[\:]?(?P<cdna>c\.[a-zA-Z0-9\+\-\_\>\?\*]+)')

DNA_REGEX = re.compile('{}{}{}{}{}{}{}{}{}{}'.format(
    r'(?:(?P<id>.*):*)?c\.',  # transcript reference
    r'(?P<start>[\d?]+)',  # start nt
    r'(?P<startNonCDS>[\+\-][\d?]+)?',  # nonCDS information
    r'(_(?P<end>[\d?]+))?',  # end nt
    r'(?P<stopNonCDS>[\+\-][\d?]+)?',  # nonCDS information
    r'(?P<ref>[ACTG])?',  # reference nt, for missense
    r'(?P<varType1>(del)?(ins)?(dup)?(>)?)',  # type of variation
    r'(?P<alt1>[ACTG]+)?',  # alternate nt
    r'(?P<varType2>(del)?(ins)?(dup)?(>)?)?',  # type of variation
    r'(?P<alt2>[ACTG]+)?',  # alternate nt
))
# these regular expressions are targeted for ensembl transcript references
# TODO: This regex was written from scratch and may be wrong. 
# It also doesn't currently support multiple variants at once.
# Find a library to replace the above regex. 
# (The hgvs python package currently doesn't work for windows)

# compiled CDS examples to test DNA_REGEX:
# ENST00000256474.2:c.1-?_642+?del
# ENST00000256474.2:c.227_229delTCT
# ENST00000256474.2:c.180delG
# ENST00000256474.2:c.445delG
# ENST00000256474.2:c.228_229insC
# ENST00000256474.2:c.478_479delGA
# ENST00000256474.2:c.481C>G
# ENST00000256474.2:c.349dupT
# ENST00000256474.2:c.516_517dupGTCAAGCCT
# ENST00000256474.2:c.216delC
# ENST00000256474.2:c.263_265delGGCinsTT 

# Q: from clinvar, what is the following variant?
# NM_000551.3(VHL):c.*2854G>T


## CDNA SNP Analysis
RING_TYPE = {
    'A': "purine",
    'C': "pyrimidine",
    'T': "pyrimidine",
    'G': "purine"
}


def TT_FUNCTION(ref, alt):
    isTransition = RING_TYPE[ref] == RING_TYPE[alt]
    isTransversion = not (RING_TYPE[ref] == RING_TYPE[alt])

    assert (isTransition != isTransversion)

    return "transition" if isTransition else "transversion"


AA_1TO3 = protein_letters_1to3
AA_1TO3['*'] = "Ter"
AA_1TO3['del'] = "del"
AA_1TO3['fs'] = "fs"
AA_1TO3['FS'] = AA_1TO3['fs']


def get_valid_cdna(cdna_str, check_version=False):
    return_cdna = None

    match = CDNA_REGEX.match(cdna_str)
    if match is not None:
        var = match.groupdict()

        cdna = var.get('cdna', None)
        # if there is some cdna change
        if cdna is not None and cdna != '':
            return_cdna = cdna

        # revert cdna to invalid if transcript version is wrong
        if check_version:
            t_id = var.get('id', '')
            # if the transcript ref is blank or not current
            if t_id == '' or (not t_id in CURRENT_VHL_TRANSCRIPT.values()):
                return_cdna = None

    return return_cdna


VHL201_FASTA = os.path.join('Analysis', 'files', 'Homo_sapiens_VHL_201_sequence.fa')

## Phenotype and Sequency Ontology Utilities
SO_NAME = 'SequenceOntology'
SO_FILENAME = 'so.obo'
# SO_HREF = 'https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so.obo'
SO_HREF = os.path.join('Analysis', 'files', 'so.obo.txt')

HPO_NAME = 'HumanPhenotypeOntology'
HPO_FILENAME = 'hp.obo'
# HPO_HREF = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'
HPO_HREF = os.path.join('Analysis', 'files', 'hp.obo.txt')

# gather HPO and SO into an netork- node keys: id
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
    tid = term_or_id.strip().casefold()
    if tid in OBONET:
        return OBONET.node[tid][obo_type]
    else:
        raise ValueError(f"Could not find an OBO node for {tid}")


## Scoring functions
# Q: hpo has a confusing hierarchy for neurendocrine neoplasms relating to the pancreas
# (namely, Pancreatic endocrine tumor vs neoplasms and cysts of the pancreas )


# GENERAL_HPO_TERMS = [
# 	'neuroendocrine neoplasm', 					# Pheochromocytoma + Paraganglioma + Pancreatic endocrine tumor + Pancreatic islet cell adenoma
# 	'renal neoplasm', 							# Renal cell carcinoma + Clear cell renal cell carcinoma
# 	'neoplasm of the central nervous system', 	# Cerebellar hemangioblastoma + Spinal hemangioblastoma
# 	'vascular neoplasm', 						# Retinal capillary hemangioma
# 	'neoplasm of the inner ear', 				# Endolymphatic sac tumor
# 	'neoplasm of the pancreas',					# Neoplasm of the pancreas
# 	'abnormal pancreas morphology',				# Pancreatic cysts
# 	'abnormal renal morphology',				# Renal cysts + Multiple Renal cysts
# 	'abnormality of the epididymis'				# Epididymal cyst
# ]

GENERAL_HPO_TERMS = [
    'neuroendocrine neoplasm',
    'renal cell carcinoma',
    'hemangioblastoma',
    'retinal capillary hemangioma',
    'pancreatic endocrine tumor',
    # 'endolymphatic sac tumor',
    'abnormality of the kidney',
    'abnormality of the pancreas',
    # 'abnormality of the epididymis',
    # 'Abnormality of the ovary'
]

# GENERAL_SO_TERMS = [
# 	'deletion',
# 	'exon_loss_variant',
# 	'missense_variant',
# 	'stop_gained',
# 	'utr_variant',
# 	'inframe_indel',
# 	'delins',
# 	'frameshift_variant',
# 	'splice_site_variant'
# ]
SO_TERM_TYPES = {
    # group a)
    'frameshift_variant': 'severe_LOF',
    'stop_gained': 'severe_LOF',
    'deletion': 'severe_LOF',
    'exon_loss_variant': 'severe_LOF',
    'start_lost': 'severe_LOF',
    'splice_site_variant': 'severe_LOF',

    # group b)
    'missense_variant': 'partial_LOF',
    'inframe_indel': 'partial_LOF',

    # group c)
    'synonymous_variant': 'minimal_LOF',
    'intron_variant': 'minimal_LOF',

    # group d)
    'utr_variant': 'misc_LOF',
    'stop_lost': 'misc_LOF',
    'delins': 'misc_LOF',

}

# TODO: this has been coded for phenotype entry, not Node
GENERAL_HPO_NODES = [get_valid_obo(term) for term in GENERAL_HPO_TERMS]


def generalized_vhl_phenotype(phenoype):
    '''Given a node, find its general disease type
    '''
    general_pheno = None

    valid_hpo = get_valid_obo(phenoype)
    for successor in nx.bfs_successors(OBONET, valid_hpo):
        try:
            i = GENERAL_HPO_NODES.index(successor[0])
            general_pheno = GENERAL_HPO_NODES[i]
            break
        except ValueError as e:
            pass

    if general_pheno is None:
        raise ValueError(f"Could not find a generalized term for {valid_hpo}")

    return general_pheno


GENERAL_SO_NODES = [get_valid_obo(term) for term in SO_TERM_TYPES.keys()]


def generalized_so_terms(so_type):
    '''Given a node, find its general so_type
    '''

    general_so = None

    valid_so = get_valid_obo(so_type)
    for successor in nx.bfs_successors(OBONET, valid_so):
        try:
            i = GENERAL_SO_NODES.index(successor[0])
            general_so = GENERAL_SO_NODES[i]
            break
        except ValueError as e:
            pass

    if general_so is None:
        raise ValueError(f"Could not find a generalized term for {valid_so}")

    return general_so


# TODO: code these

def protein_change(node):
    pass


def nucleotide_change(node):
    pass


def transition_translation(node):
    pass


def protein_substituion_score(node):
    pass


# TODO: figure out if this is needed; remove if not
# def aa_change_from_cds(cds):
# 	aa_change = None

# 	match = DNA_REGEX.match(cds)
# 		if match is not None:
# 			var = match.groupdict()

# 			#if the variant is not utr or intronic
# 			if (var['startNonCDS'] is None and var['stopNonCDS'] is None):
# 				reading_frame_i = math.floor(int(var['start'])/3)
# 				codon_i = int(var['start'])%3

# 				aa1= VHL_CDS[reading_frame_i:reading_frame_i+3].translate()
# 				aa1_i = reading_frame_i + 1
# 				# if ins, dup, or del of length != 3, then FS 

# 				start_aa_i = 
# 				stop_aa_i = math.floor(int(var['end'])/3) if var['end'] is not None else start_aa + 1

# 				start_aa = VHL_PROTEIN[start_aa_i]


# TODO: this changed to accept cdna, not a Node
def affected_domains(hgvs):
    '''Finds the affected VHL domains for a variant
    '''
    # TODO: make this less nested

    domains_affected = []

    # hgvs = node['all'].get('cdnaChange', None)
    if hgvs is not None:

        # there was a typo in Civic for VARIANT L89R(c.266T>G)
        # ENST00000256474.2:.266T>G
        match = DNA_REGEX.match(hgvs)
        if match is not None:
            var = match.groupdict()
            # if the variant is not utr or intronic
            if var['startNonCDS'] is None and var['stopNonCDS'] is None:
                try:
                    start_aa = math.floor(int(var['start']) / 3) + 1
                    stop_aa = math.floor(int(var['end']) / 3) + 1 if var['end'] is not None else start_aa + 1

                    affected_aa = range(start_aa, stop_aa)

                    for domain in GENE3D_VHL_DOMAINS[CURRENT_VHL_TRANSCRIPT['ensembl']]:
                        # if theres overlap in aa's between the variant and each domain
                        if len(list(set(GENE3D_VHL_DOMAINS[CURRENT_VHL_TRANSCRIPT['ensembl']][domain]) & set(
                                affected_aa))) > 0:
                            domains_affected.append(domain)

                    if len(domains_affected) == 0:
                        domains_affected = ["cds"]
                # couldnt cast start or end to integer
                except ValueError as e:
                    pass

    return domains_affected


def get_cdna_start(hgvs):
    match = DNA_REGEX.match(hgvs)
    cdna_start = 0
    if match is not None:
        var = match.groupdict()
        # if the variant is not utr or intronic
        if var['startNonCDS'] is None and var['stopNonCDS'] is None:
            cdna_start = int(var['start'])

    return cdna_start
