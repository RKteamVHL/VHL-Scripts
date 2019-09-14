import inspect

##Functions here should *only* take in n1 data, n2 data, and keyword args
# and must return a tuple of (fn_name, fn_score)

# given two nodes, finds the intersection over union for a specific attribute
def score_iou(n1, n2, attr_name):
	set1 = set(n1[attr_name])
	set2 = set(n2[attr_name])

	intersect = len(set1.intersection(set2))
	union = len(set1.union(set2))
	score = 0

	if not union == 0:
		score = intersect/union	

	#find better way of creating weight name
	return ("score_iou"+attr_name, score)