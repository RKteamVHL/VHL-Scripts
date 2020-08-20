import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

COUNT_THRESHOLD = 3

class GraphPlotter:
	def __init__(self, *args, **kwargs):
		self.graph = None
		self.unique_value_dict = {}
		if "graph" in kwargs:
			self.set_graph(kwargs["graph"])

	def set_graph(self, graph):
		assert(nx.Graph in type(graph).mro())
		self.graph = graph
		self.find_unique_values()

	def find_unique_values(self):
		for n, d in self.graph.nodes(data=True):
			# checks the merged data only; d is expected to have 'all' key
			if "all" in d:
				for field in d['all']:
					self.unique_value_dict[field] = self.unique_value_dict.get(field, set())

					to_add = d['all'][field]
					if not isinstance(d['all'][field], list):
						to_add = [d['all'][field]]

					for v in to_add:
						self.unique_value_dict[field].add(v)

	def count_nodes_by(self, field):
		data_sums = {}

		for node, n_data in self.graph.nodes(data=True):
			all_data =  n_data['all']
			if field in all_data:
				to_add = all_data[field]
				if not isinstance(all_data[field], list):
					to_add = [all_data[field]]

				for v in to_add:
					data_sums[v] = data_sums.get(v, 0) + 1

		_values = []
		_labels = []
		for k, v in data_sums.items():
			# filter values
			if v > COUNT_THRESHOLD:
				_values.append(v)
				_labels.append(k)


		i_sorted = np.argsort(_values, axis=0)

		sorted_values = np.take_along_axis(np.array(_values), i_sorted, 0)
		sorted_labels = np.take_along_axis(np.array(_labels), i_sorted, 0)


		# set width of bar
		barWidth = 1

		# Set position of bar on X axis
		r1 = np.arange(len(sorted_values))


		# Make the plot
		plt.bar(r1, sorted_values, width=barWidth, edgecolor='white', label='var1')


		# Add xticks on the middle of the group bars
		plt.xlabel(field, fontweight='bold')
		plt.xticks(r1, sorted_labels, rotation='vertical')

		plt.ylabel("# of variants", fontweight='bold')

		# Create legend & Show graphic
		plt.legend()
		plt.show()

	def grouped_bar(self, inner_field, outer_field, normalized=True):
		inner_labels = np.array(list(self.unique_value_dict[inner_field]))
		outer_labels = np.array(list(self.unique_value_dict[outer_field]))

		in_out = np.zeros(shape=(len(inner_labels), len(outer_labels)))


		for node, data in self.graph.nodes(data=True):
			all_data = data['all']
			valid = (inner_field in all_data) and (outer_field in all_data)

			if valid:
				outer_to_add = all_data[outer_field]
				inner_to_add = all_data[inner_field]
				if not isinstance(all_data[outer_field], list):
					outer_to_add = [all_data[outer_field]]

				if not isinstance(all_data[inner_field], list):
					inner_to_add = [all_data[inner_field]]

				for outer in outer_to_add:
					for inner in inner_to_add:
						in_out[np.argwhere(inner_labels == inner).item()][np.argwhere(outer_labels == outer).item()] += 1

		row_sorted_ind = np.argsort(np.sum(in_out, axis=1))[::-1]
		col_sorted_ind = np.argsort(np.sum(in_out, axis=0))[::-1]

		in_out_sorted = in_out[row_sorted_ind, :]
		in_out_sorted = in_out_sorted[:, col_sorted_ind]


		# filtering the correlation matrix
		row_sorted_filtered_ind  = np.any(in_out_sorted > COUNT_THRESHOLD, axis=1)
		col_sorted_filtered_ind = np.any(in_out_sorted > COUNT_THRESHOLD, axis=0)

		in_out_sorted_filtered = in_out_sorted[row_sorted_filtered_ind, :]
		in_out_sorted_filtered = in_out_sorted_filtered[:, col_sorted_filtered_ind]


		# set width of bar
		barWidth = 1/(in_out_sorted_filtered.shape[0]+1)

		bars =[]
		for row in in_out_sorted_filtered:
			bars.append(row)

		r = [None for i in range(in_out_sorted_filtered.shape[0])]
		r[0] = np.arange(in_out_sorted_filtered.shape[1])
		for i in range(1, in_out_sorted_filtered.shape[0]):
			r[i] = [x + barWidth for x in r[i-1]]

		for i in range(0, in_out_sorted_filtered.shape[0]):
			plt.bar(r[i], bars[i], width=barWidth, edgecolor='white', label=inner_labels[row_sorted_ind][i])


		# Add xticks on the middle of the group bars
		# plt.xlabel('group', fontweight='bold')
		plt.xticks([x for x in range(len(bars[0]))], outer_labels[col_sorted_ind], rotation='vertical')

		plt.title(f'# of variants vs {outer_field}[{inner_field}]')
		# Create legend & Show graphic
		plt.legend()
		plt.show()

	def correlation_matrix(self, x_field, y_field, normalized=True):
		x_labels = np.array(list(self.unique_value_dict[x_field]))
		y_labels = np.array(list(self.unique_value_dict[y_field]))
		xy = np.zeros(shape=(len(y_labels), len(x_labels)))

		for node, data in self.graph.nodes(data=True):
			all_data = data['all']
			if x_field in all_data:
				if y_field in all_data:
					x_to_add = all_data[x_field]
					y_to_add = all_data[y_field]
					if not isinstance(all_data[x_field], list):
						x_to_add = [all_data[x_field]]

					if not isinstance(all_data[y_field], list):
						y_to_add = [all_data[y_field]]


					for x_term in x_to_add:
						for y_term in y_to_add:
							xy[np.argwhere(y_labels == y_term).item()][np.argwhere(x_labels == x_term).item()] += 1

		row_sorted_ind = np.argsort(np.sum(xy, axis=1))[::-1]
		col_sorted_ind = np.argsort(np.sum(xy, axis=0))[::-1]

		xy_sorted = xy[row_sorted_ind, :]
		xy_sorted = xy_sorted[:, col_sorted_ind]


		# filtering the correlation matrix
		row_sorted_filtered_ind  = np.any(xy_sorted > COUNT_THRESHOLD, axis=1)
		col_sorted_filtered_ind = np.any(xy_sorted > COUNT_THRESHOLD, axis=0)

		xy_sorted_filtered = xy_sorted[row_sorted_filtered_ind, :]
		xy_sorted_filtered = xy_sorted_filtered[:, col_sorted_filtered_ind]

		if normalized:
			#by row sums
			xy_sorted_filtered = xy_sorted_filtered / np.abs(xy_sorted_filtered).max(axis=0)
			# by col sums
			# xy_sorted_filtered = xy_sorted_filtered / np.abs(xy_sorted_filtered).max(axis=1)
			#by row and col l1 norms
			# xy_sorted_filtered = xy_sorted_filtered / np.linalg.norm(xy_sorted_filtered, axis=0)
			# xy_sorted_filtered = xy_sorted_filtered / np.linalg.norm(xy_sorted_filtered, axis=1).transpose()
			#by total l1 norm
			# xy_sorted_filtered = xy_sorted_filtered / np.linalg.norm(xy_sorted_filtered)




		f = plt.figure(figsize=(19, 15))
		plt.matshow(xy_sorted_filtered, fignum=f.number)
		plt.xticks(range(xy_sorted_filtered.shape[1]), x_labels[col_sorted_ind], fontsize=6, rotation='vertical')
		plt.yticks(range(xy_sorted_filtered.shape[0]), y_labels[row_sorted_ind], fontsize=6, rotation=45)
		cb = plt.colorbar()
		# cb.ax.tick_params(labelsize=14)
		plt.show()
		pass

# TODO: Implement all previous visualizations
	# 1) phenotype vs # variants - done
	# 2) phenotype vs phenotype correlation matrix - done
	# 3) phenotype + mutation type vs # of variants - done
	# 4) phenotype vs mutation type correlation matrix - done
	# 5) variant location histogram
	# 6) phenotype + exon location vs # of variants
	# 7) phenotype vs exon location correlation matrix
	# 8) phenotype + domain vs # of variants - done

# TODO: Implement novel visualizations
	# 1) phenotype + cluster vs # of variants
	# 2) general phenotype + cluster vs # of variants
	# 3) variantType + cluster vs # of variants
	# 4)