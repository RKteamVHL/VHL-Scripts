from .features import *
import csv

class SummaryCalculator:
	def __init__(self):
		self.features = [
			PhenotypeFeature(),
			VariantTypeFeature(),
			DeNovoFeature(),
			SexFeature(),
			DomainFeature(),
			MultipleMutationFeature(),
			EvaluatedAgeFeature(),
			LastKnownAgeFeature(),

		]
		# # histogram features
		# # 10 years
		# max_age = 1200
		# age_interval = 120
		# self.month_ranges = {}
		# for i in range(0, max_age/age_interval):
		# 	age_bin = range(i*age_interval, (i+1)*age_interval)
		# 	self.month_ranges[age_bin] = 0



	def add_row(self, row):
		for feature in self.features:
			feature.add_row(row)

	def to_csv(self, filename):
		row_obj = {}
		for feature in self.features:
			feat_dict = feature.to_dict()
			row_obj.update(feat_dict)

		with open(filename, 'w', newline='', encoding='utf-8') as file:
			fields = list(row_obj.keys())
			writer = csv.DictWriter(file, fieldnames=fields, delimiter=",")
			writer.writeheader()
			writer.writerows([row_obj])
