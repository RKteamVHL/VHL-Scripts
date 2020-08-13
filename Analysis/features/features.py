import math
import re
from ..constants import *
from .. import variant_functions as vf
import logging
import numpy as np

EVALUATED_AGE_REGEX = re.compile("E((?P<Y>[0-9]+)Y)?((?P<M>[0-9]+)M)?")
LASTKNOWN_AGE_REGEX = re.compile("lk((?P<Y>[0-9]+)Y)?((?P<M>[0-9]+)M)?")

class Feature:
	def __init__(self, *args, name="Feature", column=None, **kwargs):
		self.name = name
		self.column = column
		self.logger = logging.getLogger(self.name)

	def filter(self, row):
		if row is not None:
			return True

	def update_by_column(self, row, column):
		category = row[column].casefold()
		self.update(row, category)

	def add_category(self, category):
		self.rows[category] = self.rows.get(category, [])

	def update(self, row, category):
		self.add_category(category)
		self.rows[category].append(row)

class NominalFeature(Feature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, *args, **kwargs)
		self.rows = {}

	def to_dict(self):
		keys = np.array(list(self.rows.keys()))
		counts = [len(v) for v in self.rows.values()]
		c_order = np.argsort(counts)
		k_order = keys[c_order]

		return {f'{self.name}_{k}': len(self.rows[k]) for k in k_order}
		# return {f'{self.name}_{k}': len(v) for k, v in self.rows.items()}


class OrdinalFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, *args, **kwargs)


class IntervalFeature(Feature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, *args, **kwargs)
		self.rows = []
		self.values = []
		self.histogram = NominalFeature(name=f'{self.name}_hist')

	def update(self, row, value):
		self.values.append(value)
		self.rows.append(row)

	def make_hist(self, n_bins=10):
		values = np.array(self.values)
		bins = np.histogram_bin_edges(values, bins=n_bins)
		values_bins = np.digitize(values, bins)

		for i in range(len(bins)-1):
			category = f'{int(bins[i])}-{int(bins[i+1])}'
			self.histogram.add_category(category)

		for i in range(len(self.rows)):
			v_bin = values_bins[i]
			if v_bin == len(bins):
				v_bin -= 1
			bin_lower = int(bins[v_bin-1])
			bin_upper = int(bins[v_bin])

			self.histogram.update(self.rows[i], f'{bin_lower}-{bin_upper}')

	def to_dict(self):
		self.make_hist()

		out_dict = self.histogram.to_dict()
		# out_dict = {}
		# out_dict[f"{self.name}_min"] = sorted_values[0]
		# out_dict[f"{self.name}_max"] = sorted_values[-1]
		#
		# for i in range(len(bins)-1):
		# 	out_dict[f"{self.name}_{int(bins[i])}-{int(bins[i+1])}"] = hist[i]

		return out_dict


class RatioFeature(IntervalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, *args, **kwargs)


class PhenotypeFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="phenotypes", *args, **kwargs)

	def add_row(self, row):
		hpo_list = re.split('[;,]', row['Phenotype'])
		for term in hpo_list:
			try:
				if term.casefold() is 'none':
					super().update(row, var_hpo)
				if term.casefold() not in NULL_TERMS:
					var_hpo = vf.generalized_vhl_phenotype(term.strip())
					super().update(row, var_hpo)
			except ValueError as e:
				self.logger.warning(repr(e))


class VariantTypeFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="variant_type", *args, **kwargs)

	def add_row(self, row):
		so_list = re.split('[;,]', row['Mutation Type'])
		for term in so_list:
			try:
				if term.casefold() not in NULL_TERMS:
					var_obo = vf.get_valid_obo(term.strip())
					super().update(row, var_obo)
			except ValueError as e:
				self.logger.warning(repr(e))

class DeNovoFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="de_novo", *args, **kwargs)

	def add_row(self, row):
		super().update_by_column(row, "Confirmed De Novo")

class SexFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="sex", *args, **kwargs)

	def add_row(self, row):
		super().update_by_column(row, "Sex")


class DomainFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="domain", *args, **kwargs)

	def add_row(self, row):
		cdna = vf.get_valid_cdna(row["Mutation Event c.DNA."])
		for domain in vf.affected_domains(cdna):
			super().update(row, domain)

class MultipleMutationFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="multiple", *args, **kwargs)

	def add_row(self, row):
		super().update_by_column(row, "Multiple Mutants in Case")



class EvaluatedAgeFeature(RatioFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="evaluated_age", *args, **kwargs)

	def add_row(self, row):
		match = EVALUATED_AGE_REGEX.search(row['Age'])
		if match is not None:
			var = match.groupdict()
			months = 0
			if var['Y'] is not None:
				months += int(var['Y']) * 12
			if var['M'] is not None:
				months += int(var['M'])

			if not(var['Y'] is None and var['M'] is None):
				super().update(row, months)


class LastKnownAgeFeature(RatioFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="last_known_age", *args, **kwargs)

	def add_row(self, row):

		lk_match = LASTKNOWN_AGE_REGEX.search(row['Age'])
		e_match = EVALUATED_AGE_REGEX.search(row['Age'])
		if lk_match is not None:
			var = lk_match.groupdict()
			months = 0
			if var['Y'] is not None:
				months += int(var['Y']) * 12
			if var['M'] is not None:
				months += int(var['M'])

			if not(var['Y'] is None and var['M'] is None):
				super().update(row, months)


		elif e_match is not None:
			var = e_match.groupdict()
			months = 0
			if var['Y'] is not None:
				months += int(var['Y']) * 12
			if var['M'] is not None:
				months += int(var['M'])
			if not(var['Y'] is None and var['M'] is None):
				super().update(row, months)

		else:
			pass
