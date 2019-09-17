import inspect
import math

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
#translation ID for the correct 213aa transcript (ENST00000256474.2),
#they have  different domains for the alpha/beta regions.
#Which one should be used? Below the most recent GRCH38.p12 are used.

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


##Functions here should *only* take in n1 data, n2 data, and keyword args
# and must return a tuple of (fn_name, fn_score)

# given two nodes, finds the intersection over union for a specific attribute
# the attributes here must be of a class type that is iterable (i.e., a list)
def score_iou(n1, n2, attr_name):
	set1 = set(n1[attr_name])
	set2 = set(n2[attr_name])

	intersect = len(set1.intersection(set2))
	union = len(set1.union(set2))
	score = 0

	if not union == 0:
		score = intersect/union	

	#find better way of creating weight name
	return ("score_iou_"+attr_name, score)

# calculates similarity between variants based on whether they're in the same
# protein domains
def variant_score_domains(n1, n2):
	pass


# calculates similarity between snvs based on how far (in aa) they are,
#according to a gaussian distribution
def variant_aa_distance(n1, n2, sigma=5):
	#TODO: get these from the nodes
	aa1 = 5
	aa2 = 10

	score = math.exp(-0.5 * (((aa1-aa2)/sigma )**2) )
	return ("aa_distance", score)


# calculates similarity between snvs based on how far (in bp) they are,
#according to a gaussian distribution
def variant_nt_distance(n1, n2, sigma=200):
	#TODO: get these from the nodes
	nt1 = 5
	nt2 = 10

	score = math.exp(-0.5 * (((nt1-nt2)/sigma )**2) )
	return ("nt_distance", score)
