from Bio import SeqIO
from Bio.Seq import Seq
import math
import os
import re

#
# The following code loads the VHL201 transcript, which contains sequences for:
# utr5-exon1-exon2-exon3-utr3
# introns are not included here
#TODO: make this less hardcoded, and function based on fasta descriptions
#The particular fasta file was obtained from: 
# https://useast.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000134086;r=3:10141008-10152220;t=ENST00000256474
VHL201_FASTA = os.path.join('ApiScripts', 'files', 'Homo_sapiens_VHL_201_sequence.fa')

vhl_seqs = [str(seqreq.seq) for seqreq in SeqIO.parse(VHL201_FASTA, "fasta")]

utr3 = Seq(vhl_seqs[5])
utr5 = Seq(vhl_seqs[6])
exon1 = Seq(vhl_seqs[0].replace(vhl_seqs[6], ''))
exon2 = Seq(vhl_seqs[1])
exon3 = Seq(vhl_seqs[2].replace(vhl_seqs[5], ''))

VHL_CDS = exon1 + exon2 + exon3

VHL_PROTEIN = VHL_CDS.translate()


#resources used for finding these values:
#	https://www.ensembl.org/Homo_sapiens/Gene/Splice?db=core;g=ENSG00000134086;r=3:10141008-10152220

# for ENSG00000134086, ensembl links to Pfam ids: 
# PF01847 (VHL beta domain); http://pfam.xfam.org/family/PF01847
# PF17211 (short C-terminal alpha helical domain); http://pfam.xfam.org/family/PF17211
# summary: http://pfam.xfam.org/protein/VHL_HUMAN

# On these pfam pages, the aa ranges from Uniprot are found on the Structures tab.
# the most conservative ranges were used for the dictionary below

#Q: Gene3D has a different set of coordinates for these- which is more
#trustworthy, Gene3D or Pfam?

#Q: Although both human ref builds (GRCh38.p12, GRCh37.p13) have the same 
#transcript ID for the current 213aa transcript (ENST00000256474.2),
#they have  different domains for the alpha/beta regions.
#Which one should be used? Below the most recent GRCH38.p12 are used.

CURRENT_VHL_TRANSCRIPT = 'ENST00000256474.2'

PFAM_VHL_DOMAINS = {
	"ENST00000256474.2": {
		# 63-143, inclusive
		"beta": range(63, 144),
		# 149-204, inclusive
		"alpha": range(149, 205)
	},

	#Q: I couldn't find these values on pfam, only on ensembl. Did 
	#ensembl do the mapping, or is it somewhere on pfam?
	"ENST00000345392.2":{
		# 63-114, inclusive
		"beta": range(63, 115),
		# 114-163, inclusive
		"alpha": range(114, 164)
	}
}

GENE3D_VHL_DOMAINS = {
	"ENST00000256474.2": {
		# 51-152, inclusive
		"beta": range(51, 153),
		# 153-204, inclusive
		"alpha": range(153, 205)
	},


	"ENST00000345392.2":{
		# 51-117, inclusive
		"beta": range(51, 118),
		# 118-163, inclusive
		"alpha": range(118, 164)
	}
}


DNA_REGEX = re.compile('{}{}{}{}{}{}{}{}{}{}'.format(
	r'(?P<id>.*):c\.', 							#transcript reference
	r'(?P<start>[\d?]+)',						#start nt
	r'(?P<startNonCDS>[\+\-][\d?]+)?',			#nonCDS information
	r'(_(?P<end>[\d?]+))?',						#end nt
	r'(?P<stopNonCDS>[\+\-][\d?]+)?',			#nonCDS information
	r'(?P<ref>[ACTG])?', 						#reference nt, for missense
	r'(?P<varType1>(del)?(ins)?(dup)?(>)?)', 	#type of variation
	r'(?P<alt1>[ACTG]+)?', 						#alternate nt
	r'(?P<varType2>(del)?(ins)?(dup)?(>)?)?', 	#type of variation
	r'(?P<alt2>[ACTG]+)?', 						#alternate nt
))
#these regular expressions are targeted for ensembl transcript references
#TODO: This regex was written from scratch and may be wrong. 
# It also doesn't currently support multiple variants at once.
# Find a library to replace the above regex. 
# (The hgvs python package currently doesn't work for windows)


#NOTE: these were the old regexes before the one above was written
_DNA_REGEX =	{

	#TODO: try make this a 'catch-all' if possible
	'default': re.compile(r'([ACTG])>([ACTG])'),

	#eg. ENST00000256474.2:c.481C>G
	'missense_variant': re.compile(r'(.*):c\.([0-9]+)([ACTG])>([ACTG])'),
	'stop_gained': re.compile(r'(.*):c\.([0-9]+)([ACTG])>([ACTG])'),

	#eg. ENST00000256474.2:c.227_229delTCT, ENST00000256474.2:c.180delG
	'inframe_deletion': re.compile(r'(.*):c\.([0-9]+)(?:_([0-9]+))?del([ACTG]+)'),

	#eg. ENST00000256474.2:c.445delG, ENST00000256474.2:c.228_229insC, ENST00000256474.2:c.478_479delGA
	'frameshift_truncation': re.compile(r'(.*):c\.([0-9]+)(?:_([0-9]+))?(?:del)?(?:ins)?([ACTG]+)'),

	#eg. ENST00000256474.2:c.375_376insC, ENST00000256474.2:c.349dupT
	'inframe_insertion': re.compile(r'(.*):c\.([0-9]+)(?:_([0-9]*))?(?:dup)?(?:ins)?([ACTG]+)')
}

#compiled CDS examples to test DNA_REGEX:
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


def affected_domains(node):
	domains_affected = []

	#for simple missense
	for expression in node['hgvsExpressions']:
		if expression.find(CURRENT_VHL_TRANSCRIPT)>-1:

			#there was a typo in Civic for VARIANT L89R(c.266T>G)
			# ENST00000256474.2:.266T>G
			match = DNA_REGEX.match(expression)
			if match is not None:
				var = match.groupdict()
				#if the transcript isn't the current one
				if var['id'] in GENE3D_VHL_DOMAINS:		

					#if the variant is not utr or intronic
					if (var['startNonCDS'] is None and var['stopNonCDS'] is None):
						start_aa = math.floor(int(var['start'])/3)
						stop_aa = math.floor(int(var['end'])/3) if var['end'] is not None else start_aa + 1

						affected_aa = range(start_aa, stop_aa)

						for domain in GENE3D_VHL_DOMAINS[var['id']]:
							#if theres overlap in aa's
							if len(list(set(GENE3D_VHL_DOMAINS[var['id']][domain]) & set(affected_aa))) > 0:
								domains_affected.append(domain)

	return ("affected_domains", domains_affected)
