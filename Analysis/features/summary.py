from .KimStudents import *
from pandas import DataFrame


class Summary:
	def __init__(self):
		self.name = "Summary"
		self.features = [
			PhenotypeFeature(),
			VariantTypeFeature(),
			DeNovoFeature(),
			SexFeature(),
			DomainFeature(),
			MultipleMutationFeature(),
			EvaluatedAgeFeature(),
			LastKnownAgeFeature()
		]
		self.dataframe = DataFrame()

	def add_row(self, row):
		for feature in self.features:
			feature.add_row(row)

	def to_csv(self, filename):
		self.dataframe.to_csv(filename)

	def normalize(self, norm_types=None):
		pass

	def make_dataframe(self):
		row_obj = {}
		for feature in self.features:
			feat_dict = {k: len(v) for k, v in feature.to_dict().items()}
			row_obj.update(feat_dict)

		self.dataframe = self.dataframe.append(row_obj, ignore_index=True)
		self.dataframe = self.dataframe.fillna(0)
