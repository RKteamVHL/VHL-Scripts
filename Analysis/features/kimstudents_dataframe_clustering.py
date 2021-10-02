import snf
# from sklearn.cluster import spectral_clustering
from sklearn.cluster import SpectralClustering


from .kimstudents_dataframe_preprocessing import COMPUTED_COLUMNS


def dataframe_snf(df):


    col_names = ["generalized_phenotype", "generalized_mutant_type", "codon"]

    col_groups = [COMPUTED_COLUMNS[col_type] for col_type in col_names]
    all_cols = []
    for col_group in col_groups:
        for col in col_group:
            all_cols.append(col)
    # df = df.dropna(subset=COMPUTED_COLUMNS["codon"]).sort_values(by="codon_start")
    df_cols = df[all_cols]

    # codon is not filled as 0's because it is a ratio value, not nominal like pheno or muttype
    df_cols[COMPUTED_COLUMNS["generalized_phenotype"]] = df_cols[COMPUTED_COLUMNS["generalized_phenotype"]].fillna(0)
    df_cols[COMPUTED_COLUMNS["generalized_mutant_type"]] = df_cols[COMPUTED_COLUMNS["generalized_mutant_type"]].fillna(0)


    #converting to bools:
    df_cols = df_cols >= 1

    feat_metrics = []
    for col_group in col_groups:
        df_group = df_cols[col_group]

        # affinity = snf.make_affinity(df_group.to_numpy(), metric='sqeuclidean', K=len(df_cols.index), mu=0.5)
        affinity = snf.make_affinity(df_group.to_numpy(), normalize=False,  K=len(df_cols.index),
                                     metric=['jaccard', 'jaccard', 'sqeuclidean'], mu=0.7)
        feat_metrics.append(affinity)

        pass

    fused_net = feat_metrics[0]
    #uncomment the following to do snf
    # if len(feat_metrics) == 1:
    #     fused_net =feat_metrics[0]
    # else:
    #     fused_net = snf.snf(feat_metrics, K=50)
        # fused_net = snf.snf(feat_metrics, K=len(df_cols.index))

    best, second = snf.get_n_clusters(fused_net)
    sc1 = SpectralClustering(best, affinity='precomputed', n_init=100, assign_labels='discretize')
    labels1 = sc1.fit_predict(fused_net) + 1
    df["cluster_labels_best"] = labels1

    sc2 = SpectralClustering(second, affinity='precomputed', n_init=100, assign_labels='discretize')
    labels2 = sc2.fit_predict(fused_net) + 1
    df["cluster_labels_second"] = labels2

    df = df.sort_values(by="cluster_labels_best")


    #plotting affinity mats
    # plt.close('all')
    # sorted_i = np.argsort(labels1, axis=0)
    # for i in range(len(col_names)):
    #     sorted_fused = feat_metrics[i][sorted_i, :]
    #     sorted_fused = sorted_fused[:, sorted_i]
    #     ax = plt.imshow(sorted_fused, cmap='hot', interpolation='nearest')
    #     plt.savefig(os.path.join(directory, f'sorted_{col_names[i]}_affinities.pdf'))
    #     plt.close('all')

    return df
