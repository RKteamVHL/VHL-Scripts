from .VariantGraph import VariantGraph
import matplotlib.pyplot as plt
import argparse
from snf import compute

COLORMAP= {
	0:(0,0,0,0.5), 
	1:(0,0,1,0.5),
	2:(0,1,0,0.5),
	3:(0,1,1,0.5),
	4:(0,1,0,0.5),
	5:(0,1,1,0.5),
	6:(1,0,0,0.5),
	7:(1,0,1,0.5),

	
}

FILE_OUTPUT = "output/"

PRE_CACHE = "cache/variant_nodes.json"
POST_CACHE = "cache/variant_nodes_with_similarities.json"



if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('similarity_type', help = '''The type of similarity to base clustering on''')
	parser.add_argument('-d', '--directory', help = '''Local directory to save graph''', default="")
	parser.add_argument('-f', '--filter', help='''Removes nodes that don't have values for the filter key''')
	parser.add_argument('-c', '--cluster', help='''Forces spectral clustering to use this many clusters''')
	parser.add_argument('-rc', '--pre_cache', help = '''Use the local graph raw node cache''', action="store_true")
	parser.add_argument('-oc', '--post_cache', help = '''Use the local graph calulated node cache''', action="store_true")
	parser.add_argument('-is', '--ignore_submitted', help = '''If set, ignores unreviewed variants (civic)''', action="store_true")

	args = parser.parse_args()

	VG = VariantGraph()

	#fetch/process all relevant data from all sources
	if not args.pre_cache:

		VG.add_nodes_from_db('Civic')
		VG.add_nodes_from_db('KimStudents2019')
		VG.add_nodes_from_db('Gnomad')
		VG.add_nodes_from_db('ClinVar')
		VG.merge_nodes()
		VG.save_to_json_file(PRE_CACHE)

	else:
		VG.load_from_json_file(PRE_CACHE)


	# NOTE: these lines for making new graph PG are temporary,
	# and should only be used for filtering the nodes
	PG = VariantGraph()

	# PG.add_nodes_from([(n, d) for n, d in VG.nodes(data=True) if d['all']['associatedPhenotypes']]) # needs pheno
	# PG.add_nodes_from([(n, d) for n, d in VG.nodes(data=True) if d['all']['variantTypes']]) #needs var type
	if args.filter:
		PG.add_nodes_from([(n, d) for n, d in VG.nodes(data=True) if d['all'][args.filter]]) # needs filter
	else:
		PG = VG # use everything

	PG.add_edges_from([(n1, n2, d) for n1, n2, d in VG.edges(data=True) if n1 in PG and n2 in PG])


	# PG.save_to_json_file("filtered_variant_nodes.json")
	PG.save_to_text_file("variant_names.tsv")
	if not args.post_cache:
		PG.calculate_node_attributes()
		PG.calculate_similarities()
		PG.save_to_json_file(POST_CACHE)
	else:
		PG.load_from_json_file(POST_CACHE)

	adj_mat = PG.get_adjacency_mats(types=[args.similarity_type])

	clust_count1, clust_count2 = compute.get_n_clusters(adj_mat[0])

	if args.cluster:
		clust_count1 = int(args.cluster)
		print(f"Forcing {clust_count1} clusters")
	else:
		print(f"Cluster estimates: {clust_count1}, {clust_count2}")



	PG.cluster_by(args.similarity_type, num_clusters=clust_count1)
	PG.aggregate_cluster(args.similarity_type)
	PG.save_to_json_file(f"{FILE_OUTPUT}{args.similarity_type}_variant_nodes.json", nodes_only=True)
	# Get newly sorted adjmat
	sorted_adj_mat = PG.get_adjacency_mats(types=[args.similarity_type])
	plt.matshow(sorted_adj_mat[0])
	plt.savefig(f'{FILE_OUTPUT}{args.similarity_type}_heatmap')
	plt.show()




	# plotting the network isn't too informative
	# nx.draw(PG, 
	# 	node_size=100, 
	# 	with_labels=True,		
	# 	pos=nx.spring_layout(PG),
	# 	width=0.001,
	# 	node_color=[COLORMAP[d[f'{args.similarity_type}_label']] for (u,d) in PG.nodes(data=True)]
	# )
	# plt.show()


# QUESTIONS:
# - where do you find the keys each attribute type has?
# - is there a list of ids for mutation types (Civic Attribute)