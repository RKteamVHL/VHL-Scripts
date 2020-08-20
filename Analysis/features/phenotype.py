from .KimStudents import *
import pandas as pd


class Phenotype:
	def __init__(self):
		self.name = "Phenotype"
		self.main_feature = PhenotypeFeature()
		self.dataframe = pd.DataFrame()

	def add_row(self, row):
		self.main_feature.add_row(row)

	def to_csv(self, filename):
		self.dataframe.to_csv(filename)

	def normalize(self, norm_types=["minmax"]):
		combined = pd.DataFrame()
		for norm_type in norm_types:
			df = self.dataframe.select_dtypes(include=[np.number])
			if norm_type == "minmax":
				df = (df - df.min()) / (df.max() - df.min())
			elif norm_type == "zscore":
				df = (df - df.mean()) / df.std()
			else:
				return
			df.columns = [f'{norm_type}_{colname}' for colname in df.columns]
			combined = pd.concat([combined, df], sort=False, axis=1)
		self.dataframe = pd.concat([self.dataframe, combined], sort=False, axis=1)

	def make_dataframe(self):

		for category, row_list in self.main_feature.rows.items():
			row_dict = {"phenotype": category, "row_count": len(row_list)}
			row_features = [
				VariantTypeFeature(),
				SexFeature(),
				DomainFeature()
			]

			for row in row_list:
				for i in range(len(row_features)):
					row_features[i].add_row(row)

			for feat in row_features:
				feat_dict = {k: len(v) for k, v in feat.to_dict().items()}
				row_dict.update(feat_dict)

			self.dataframe = self.dataframe.append(row_dict, ignore_index=True)

		self.dataframe = self.dataframe.fillna(0)

