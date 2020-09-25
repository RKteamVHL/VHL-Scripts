import pandas as pd
import snf
from matplotlib import pyplot
from sklearn.cluster import SpectralClustering

from .FeatureTable import FeatureTable
from .KimStudentsFeatures import *


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
        super().__init__(self, name="EvaluatedAgeBinned", main_feature=main_feat, row_classes=row_classes, *args,
                         **kwargs)

    def make_dataframe(self, *args, **kwargs):
        self.main_feature.make_hist()
        self.main_feature = self.main_feature.histogram
        super().make_dataframe(*args, **kwargs)


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

    def make_dataframe(self, *args, **kwargs):
        self.main_feature.sort()
        super().make_dataframe(*args, **kwargs)


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
        super().__init__(self, name="LastKnownAgeBinned", main_feature=main_feat, row_classes=row_classes, *args,
                         **kwargs)

    def make_dataframe(self, *args, **kwargs):
        self.main_feature.make_hist()
        self.main_feature = self.main_feature.histogram
        super().make_dataframe(*args, **kwargs)


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

    def make_dataframe(self, *args, **kwargs):
        self.main_feature.sort()
        super().make_dataframe(*args, **kwargs)


class IsolatedPhenotypeTable(FeatureTable):
    def __init__(self, *args, **kwargs):
        main_feat = IsolatedPhenotypeFeature()
        row_classes = [
            VariantTypeFeature,
            SexFeature,
            DomainWeightedFeature
        ]
        super().__init__(self, name="IsolatedPhenotype", main_feature=main_feat, row_classes=row_classes, *args,
                         **kwargs)


class VariantTypeTable(FeatureTable):
    def __init__(self, *args, **kwargs):
        main_feat = VariantTypeFeature()
        row_classes = [
            PhenotypeFeature,
            # SexFeature,
            # DomainWeightedFeature
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
        super().__init__(self, name="PhenotypeCoocurrence", main_feature=main_feat, row_classes=row_classes, *args,
                         **kwargs)

    def make_dataframe(self, *args, **kwargs):
        super().make_dataframe(*args, **kwargs)
        feat_name = self.main_feature.name
        columns = [f'{feat_name}_{phen}' for phen in self.dataframe[feat_name]]
        # columns.insert(0, feat_name)
        self.dataframe = self.dataframe[columns]

    def cluster(self, out_dir=None, num_clusters=5):

        df_mat = self.dataframe.to_numpy()
        norm_df_mat = np.zeros_like(df_mat)
        for i in range(df_mat.shape[0]):
            for j in range(df_mat.shape[1]):
                norm_df_mat[i, j] = df_mat[i, j] / np.max([np.max(df_mat[i]), np.max(df_mat[j])])

        col_labels = np.array(self.dataframe.columns)

        if num_clusters is None:
            best, second = snf.get_n_clusters(norm_df_mat)
            num_clusters = best
        # cluster based on given adjmat and # clusters
        sc = SpectralClustering(num_clusters, affinity='precomputed', n_init=1000)

        sc.fit(norm_df_mat)

        fused_labels = sc.labels_

        l_ids = np.argsort(fused_labels)
        df_mat_sorted = norm_df_mat[l_ids, :]
        df_mat_sorted = df_mat_sorted[:, l_ids]
        labels_out = np.array([col_labels[l_ids], fused_labels[l_ids]])

        # adj_plot = pyplot.matshow(self.dataframe.to_numpy())
        if out_dir is not None:
            pd.DataFrame(labels_out).to_csv(os.path.join(out_dir, f'{self.name}_labels.csv'), index=False)

            mat_plot = pyplot.matshow(df_mat_sorted)
            pyplot.savefig(figure=mat_plot, fname=os.path.join(out_dir, f'{self.name}_adj.png'))

    def normalize(self, norm_types=["minmax"]):
        pass


class LossOfFunctionTable(FeatureTable):
    def __init__(self, *args, **kwargs):
        main_feat = LossOfFunctionFeature()
        row_classes = [
            PhenotypeFeature,
        ]
        super().__init__(self, name="VariantLOF", main_feature=main_feat, row_classes=row_classes, *args, **kwargs)
