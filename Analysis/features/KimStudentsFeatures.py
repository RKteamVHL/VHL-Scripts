from .Feature import *
from ..constants import *
from .. import variant_functions as vf

EVALUATED_AGE_REGEX = re.compile("E((?P<Y>[0-9]+)Y)?((?P<M>[0-9]+)M)?")
LASTKNOWN_AGE_REGEX = re.compile("lk((?P<Y>[0-9]+)Y)?((?P<M>[0-9]+)M)?")


class SummaryFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="summary", *args, **kwargs)

	def add_row(self, row):
		super().update(row, "totals")


class MutationEventFeature(OrdinalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="mutation_event", *args, **kwargs)

	def add_row(self, row):
		mutlist = re.split('[;]', row['Mutation Event c.DNA.'])

		try:
			if len(mutlist) > 1:
				assert(row['Multiple Mutants in Case'].strip().casefold() == "yes")
		except ValueError as e:
			self.logger.warning(repr(e))

		for term in mutlist:
			term = term.casefold().strip()
			if term not in NULL_TERMS:
				super().update(row, term)

	def sort(self):
		# sorts rows by cdna position

		keys = np.array(list(self.rows.keys()))
		cdnas = [vf.get_cdna_start(k) for k in self.rows.keys()]
		c_order = np.argsort(cdnas)
		k_order = keys[c_order]
		old_rows = self.rows
		self.rows = {k: old_rows[k] for k in k_order}

	def to_dict(self):
		self.sort()
		super().to_dict()


class ResolutionFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="resolution", *args, **kwargs)

	def add_row(self, row):
		super().update_by_column(row, "Resolution")


class PhenotypeFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="phenotypes", *args, **kwargs)

	def add_row(self, row):
		hpo_list = re.split('[;,]', row['Phenotype'])
		row_hpos = {}
		for term in hpo_list:
			term = term.casefold().strip()
			try:

				if term not in NULL_TERMS:
					var_hpo = vf.generalized_vhl_phenotype(term)
					row_hpos[var_hpo] = row_hpos.get(var_hpo, 0) + 1
				elif term.casefold() is 'none':
					row_hpos[term] = row_hpos.get(term, 0) + 1

			except ValueError as e:
				self.logger.warning(repr(e))

		for term in row_hpos.keys():
			super().update(row, term)

class IsolatedPhenotypeFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="isolated_phenotypes", *args, **kwargs)

	def add_row(self, row):
		hpo_list = re.split('[;,]', row['Phenotype'])
		for term in hpo_list:
			term = term.casefold().strip()
			try:
				if term not in NULL_TERMS:
					var_hpo = vf.generalized_vhl_phenotype(term)
					# if this is the only phenotype in the column
					if len(hpo_list) == 1:
						super().update(row, var_hpo)
			except ValueError as e:
				self.logger.warning(repr(e))


class VariantTypeFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="variant_type", *args, **kwargs)

	def add_row(self, row):
		so_list = re.split('[;,]', row['Mutation Type'])
		for term in so_list:
			term = term.casefold().strip()
			try:
				if term not in NULL_TERMS:
					var_obo = vf.generalized_so_terms(term)
					super().update(row, var_obo)
			except ValueError as e:
				self.logger.warning(repr(e))


class VariantTypeUngroupedFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="variant_type_ungrouped", *args, **kwargs)

	def add_row(self, row):
		super().update_by_column(row, 'Mutation Type')


class LossOfFunctionFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="lof", *args, **kwargs)

	def add_row(self, row):
		so_list = re.split('[;,]', row['Mutation Type'])
		for term in so_list:
			term = term
			try:
				if term not in NULL_TERMS:
					var_obo = vf.generalized_so_terms(term)
					lof_type = vf.SO_TERM_TYPES[var_obo]
					super().update(row, lof_type)
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

class DomainWeightedFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="domain", *args, **kwargs)

	def add_row(self, row):
		cdna = vf.get_valid_cdna(row["Mutation Event c.DNA."])
		for domain in vf.affected_domains(cdna):
			super().update(row, domain)

	def to_dict(self):
		return_dict = {}
		for key in self.rows.keys():
			if key == "alpha":
				return_dict[key] = len(self.rows[key]) * (vf.CDS_LEN/vf.ALPHA_LEN)
			elif key == "beta":
				return_dict[key] = len(self.rows[key]) * (vf.CDS_LEN / vf.BETA_LEN)
			elif key == "cds":
				return_dict[key] = len(self.rows[key]) * (vf.CDS_LEN / (vf.ALPHA_LEN + vf.BETA_LEN))
		return return_dict


class MultipleMutationFeature(NominalFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="multiple", *args, **kwargs)

	def add_row(self, row):
		super().update_by_column(row, "Multiple Mutants in Case")


class EvaluatedAgeFeature(RatioFeature):
	def __init__(self, *args, **kwargs):
		super().__init__(self, name="evaluated_age", *args, **kwargs)

	def make_hist(self, bins=range(0, 1200, 120)):
		super().make_hist(bins=bins)


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

	def make_hist(self, bins=range(0, 1200, 120)):
		super().make_hist(bins=bins)

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
