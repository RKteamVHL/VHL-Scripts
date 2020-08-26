from .KimStudentsFeatures import *
from .FeatureTable import FeatureTable
import pandas as pd


class EvaluatedAgeTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = EvaluatedAgeFeature()
		row_classes = [
			IsolatedPhenotypeFeature,
			PhenotypeFeature,
			VariantTypeFeature,
			SexFeature,
			DomainFeature
		]
		super().__init__(self, name="EvaluatedAge", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)

	# self.main_feature = LastKnownAgeFeature()

	def make_dataframe(self):
		self.main_feature.make_hist()
		self.main_feature = self.main_feature.histogram
		super().make_dataframe()

class EvaluatedAgeUnbinnedTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = EvaluatedAgeFeature()
		row_classes = [
			IsolatedPhenotypeFeature,
			PhenotypeFeature,
			VariantTypeFeature,
			SexFeature,
			DomainFeature
		]
		super().__init__(self, name="EvaluatedAge", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)

	# self.main_feature = LastKnownAgeFeature()

	def make_dataframe(self):
		self.main_feature.make_hist()
		self.main_feature = self.main_feature.histogram
		super().make_dataframe()


class LastKnownAgeTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = LastKnownAgeFeature()
		row_classes = [
			IsolatedPhenotypeFeature,
			PhenotypeFeature,
			VariantTypeFeature,
			SexFeature,
			DomainFeature
		]
		super().__init__(self, name="LastKnownAge", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)

	# self.main_feature = LastKnownAgeFeature()

	def make_dataframe(self):
		self.main_feature.make_hist()
		self.main_feature = self.main_feature.histogram
		super().make_dataframe()


class PhenotypeTable(FeatureTable):

	def __init__(self, *args, **kwargs):
		main_feat = PhenotypeFeature()
		row_classes = [
			VariantTypeFeature,
			SexFeature,
			DomainFeature
		]
		super().__init__(self, name="Phenotype", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)


class IsolatedPhenotypeTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = IsolatedPhenotypeFeature()
		row_classes = [
			VariantTypeFeature,
			SexFeature,
			DomainFeature
		]
		super().__init__(self, name="IsolatedPhenotype", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)


class VariantTypeTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = VariantTypeFeature()
		row_classes = [
				PhenotypeFeature,
				LastKnownAgeFeature,
				SexFeature,
				DomainFeature
		]
		super().__init__(self, name="VariantType", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)


class SummaryTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = ResolutionFeature()
		row_classes = [
			IsolatedPhenotypeFeature,
			PhenotypeFeature,
			VariantTypeFeature,
			DeNovoFeature,
			SexFeature,
			DomainFeature,
			MultipleMutationFeature,
			EvaluatedAgeFeature,
			LastKnownAgeFeature
		]
		super().__init__(self, name="Summary", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)

	def normalize(self, norm_types=["minmax"]):
		pass
