from .KimStudentsFeatures import *
from .FeatureTable import FeatureTable
from sklearn.cluster import SpectralClustering
import pandas as pd


class EvaluatedAgeBinnedTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = EvaluatedAgeFeature()
		row_classes = [
			IsolatedPhenotypeFeature,
			PhenotypeFeature,
			VariantTypeFeature,
			SexFeature,
			DomainFeature
		]
		super().__init__(self, name="EvaluatedAgeBinned", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)

	def make_dataframe(self):
		self.main_feature.make_hist()
		self.main_feature = self.main_feature.histogram
		super().make_dataframe()


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

	def make_dataframe(self):
		self.main_feature.sort()
		super().make_dataframe()

class LastKnownAgeBinnedTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = LastKnownAgeFeature()
		row_classes = [
			IsolatedPhenotypeFeature,
			PhenotypeFeature,
			VariantTypeFeature,
			SexFeature,
			DomainFeature
		]
		super().__init__(self, name="LastKnownAgeBinned", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)

	def make_dataframe(self):
		self.main_feature.make_hist()
		self.main_feature = self.main_feature.histogram
		super().make_dataframe()


class PhenotypeTable(FeatureTable):

	def __init__(self, *args, **kwargs):
		main_feat = PhenotypeFeature()
		row_classes = [
			LossOfFunctionFeature,
			VariantTypeFeature,
			SexFeature,
			DomainWeightedFeature
		]
		super().__init__(self, name="Phenotype", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)


class MutationEventTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = MutationEventFeature()
		row_classes = [
			LossOfFunctionFeature,
			VariantTypeFeature,
			SexFeature,
			DomainWeightedFeature
		]
		super().__init__(self, name="MutationEvent", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)

	def make_dataframe(self):
		self.main_feature.sort()
		super().make_dataframe()

class IsolatedPhenotypeTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = IsolatedPhenotypeFeature()
		row_classes = [
			VariantTypeFeature,
			SexFeature,
			DomainWeightedFeature
		]
		super().__init__(self, name="IsolatedPhenotype", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)


class VariantTypeTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = VariantTypeUngroupedFeature()
		row_classes = [
				PhenotypeFeature,
				SexFeature,
				DomainWeightedFeature
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
			DomainWeightedFeature,
			MultipleMutationFeature,
			EvaluatedAgeFeature,
			LastKnownAgeFeature
		]
		super().__init__(self, name="Summary", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)

	def normalize(self, norm_types=["minmax"]):
		pass


class PhenotypeCoocurrenceTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = PhenotypeFeature()
		row_classes = [
			PhenotypeFeature,
		]
		super().__init__(self, name="PhenotypeCoocurrence", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)

	def make_dataframe(self):
		super().make_dataframe()
		feat_name = self.main_feature.name
		columns = [f'{feat_name}_{phen}' for phen in self.dataframe[feat_name]]
		# columns.insert(0, feat_name)
		self.dataframe = self.dataframe[columns]

	def cluster(self, num_clusters=4):
		# cluster based on given adjmat and # clusters
		sc = SpectralClustering(num_clusters, affinity='precomputed', n_init=1000)

		sc.fit(self.dataframe.to_numpy())

		fused_labels = sc.labels_

	def normalize(self, norm_types=["minmax"]):
		pass

class LossOfFunctionTable(FeatureTable):
	def __init__(self, *args, **kwargs):
		main_feat = LossOfFunctionFeature()
		row_classes = [
			PhenotypeFeature,
		]
		super().__init__(self, name="VariantLOF", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)
