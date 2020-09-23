from .KimStudentsFeatures import *
import pandas as pd


class FeatureTable:
	def __init__(self, *args, name="FeatureTable", main_feature=None, row_classes=[], **kwargs):
		self.name = name
		self.main_feature = main_feature
		self.row_feature_classes = row_classes
		self.row_feature_columns = [{} for feat in self.row_feature_classes]

		self.dataframe = pd.DataFrame()

	def add_row(self, row):
		if self.main_feature is not None:
			self.main_feature.add_row(row)

	def to_csv(self, filename):
		self.dataframe.to_csv(filename, index=False)

	def cluster(self, out_dir=None, num_clusters=2):
		pass

	def normalize(self, norm_types=["minmax"]):
		combined = pd.DataFrame()
		for norm_type in norm_types:
			df = self.dataframe.select_dtypes(include=[np.number])
			if norm_type == "minmax":
				df = (df - df.min()) / (df.max() - df.min())
			elif norm_type == "zscore":
				df = (df - df.mean()) / df.std()
			elif norm_type == "ratio":
				df = df / df.sum()
			elif norm_type == "cdf":
				df = df / df.sum()
				# CDF
				df = df.cumsum()
			else:
				return
			df.columns = [f'{norm_type}_{colname}' for colname in df.columns]
			combined = pd.concat([combined, df], sort=False, axis=1)
		self.dataframe = pd.concat([self.dataframe, combined], sort=False, axis=1)

	def make_dataframe(self, include_count=True):
		for category, row_list in self.main_feature.rows.items():
		# for category, row_list in self.main_feature.rows.items():
			row_dict = {self.main_feature.name: category, "row_count": len(row_list)}

			row_features = [feat_class() for feat_class in self.row_feature_classes]
			for i in range(len(row_features)):
				for row in row_list:
					row_features[i].add_row(row)

			for i in range(len(row_features)):
				_d = row_features[i].to_dict()
				feat_dict = {f'{row_features[i].name}_{k}': v for k, v in _d.items()}
				for key in feat_dict.keys():
					# use a dict here b/c they keep key order
					self.row_feature_columns[i][key] = None
				row_dict.update(feat_dict)

			self.dataframe = self.dataframe.append(row_dict, ignore_index=True)

		all_cols = [self.main_feature.name]

		if include_count:
			all_cols.append("row_count")
		for feat_cols in self.row_feature_columns:
			all_cols.extend(list(feat_cols.keys()))

		self.dataframe = self.dataframe.fillna(0)

		self.dataframe = self.dataframe[all_cols]

