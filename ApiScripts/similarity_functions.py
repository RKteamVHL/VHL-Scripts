import math


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

# calculates similarity between variants based on whether they affect the
# same protein domains
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
