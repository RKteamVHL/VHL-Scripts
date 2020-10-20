from .kimstudents_dataframe_preprocessing import COMPUTED_COLUMNS
from . kimstudents_dataframe_views import plot_clustered_stacked
import snf
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from .. import variant_functions as vf
# from sklearn.cluster import spectral_clustering
from sklearn.cluster import SpectralClustering

FIGURE_DIR = "figures"
DOMAIN_TICKS = [1, 62, 154, 192, 204, 213]

def dataframe_snf(df, analysis_type):
    fig_path = os.path.join(FIGURE_DIR, analysis_type)
    if not os.path.isdir(fig_path):
        os.makedirs(fig_path)

    col_groups = [
        COMPUTED_COLUMNS["generalized_phenotype"],
        # COMPUTED_COLUMNS["generalized_mutant_type"],
        # COMPUTED_COLUMNS["codon"]
    ]
    all_cols = []
    for col_group in col_groups:
        for col in col_group:
            all_cols.append(col)
    # df = df.dropna(subset=COMPUTED_COLUMNS["codon"]).sort_values(by="codon_start")
    df_cols = df[all_cols].fillna(0)



    feat_metrics = []
    for col_group in col_groups:
        df_group = df_cols[col_group]

        affinity = snf.make_affinity(df_group.to_numpy(), metric='sqeuclidean', K=len(df_cols.index), mu=0.5)
        # affinity = snf.make_affinity(df_group.to_numpy(), metric='cosine', K=20, mu=0.5)
        feat_metrics.append(affinity)

        pass
    fused_net = None
    if len(feat_metrics) == 1:
        fused_net =feat_metrics[0]
    else:
        fused_net = snf.snf(feat_metrics, K=20)
# fused_net = snf.snf(feat_metrics, K=len(df_cols.index))
    fused_df = pd.DataFrame(data=fused_net, index=df.index.to_list(), columns=df.index.to_list())



    n_clust = 4
    # try:
    #     best, second = snf.get_n_clusters(fused_net)
    #     n_clust = best
    # except np.linalg.LinAlgError as e:
    #     print(repr(e))


    best, second = snf.get_n_clusters(fused_net)
    n_clust = best

    sc = SpectralClustering(n_clust, affinity='precomputed', n_init=100, assign_labels='discretize')
    labels = sc.fit_predict(fused_net)
    df["cluster_labels"] = labels
    df = df.sort_values(by="cluster_labels")
    df.to_csv(os.path.join(fig_path, f"clustered_out.tsv"), sep='\t')


    plt.close('all')
    # fused_net = fused_net - np.diag(np.diag(fused_net))
    ax = plt.imshow(fused_net, cmap='hot', interpolation='nearest')
    plt.savefig(os.path.join(fig_path, f'affinities.pdf'))


    phens = df.set_index("cluster_labels")[COMPUTED_COLUMNS["generalized_phenotype"]].fillna(0)
    plt.close('all')

    pheno_sums = phens.sum().sort_values()

    pheno_counts = phens.groupby("cluster_labels").sum()

    pheno_means = phens.groupby("cluster_labels").mean()
    pheno_means = pheno_means[pheno_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1])
    pheno_stds = phens.groupby("cluster_labels").std()
    pheno_stds = pheno_stds[pheno_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1])
    # pheno_means.plot(kind="bar", yerr=pheno_stds, figsize=(12,8))
    pheno_counts.plot(kind="bar", figsize=(12,8))
    plt.savefig(os.path.join(fig_path, f'clustered_pheno_means.pdf'))
    plt.close('all')


    muttypes = df.set_index("cluster_labels")[COMPUTED_COLUMNS["generalized_mutant_type"]].fillna(0)

    muttype_sums = muttypes.sum().sort_values()

    muttype_counts = muttypes.groupby("cluster_labels").sum()

    muttype_means = muttypes.groupby("cluster_labels").mean()
    muttype_means = muttype_means[muttype_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1])
    muttype_stds = muttypes.groupby("cluster_labels").std()
    muttype_stds = muttype_stds[muttype_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1])
    # muttype_means.plot(kind="bar", yerr=muttype_stds, figsize=(12,8))
    muttype_counts.plot(kind="bar", figsize=(12,8))
    plt.savefig(os.path.join(fig_path, f'clustered_muttype_means.pdf'))
    plt.close('all')

    ix = list(df['cluster_labels'].unique())
    alpha_cols = ["⍺-Domain", "ElonginB_ElonginC_binding"]
    beta_cols = ["β-Domain", "HIF1_alpha_binding"]
    cds_cols = ["Outside of ⍺-Domain and β-Domain", "GXEEX8"]
    clust_alpha = pd.DataFrame(columns=alpha_cols, index=ix)
    clust_beta = pd.DataFrame(columns=beta_cols, index=ix)
    clust_cds = pd.DataFrame(columns=cds_cols, index=ix)


    for cluster in ix:
        df_clust = df[df['cluster_labels'] == cluster]

        clust_alpha.loc[cluster, "⍺-Domain"] = len(df_clust[(df_clust['region.⍺-Domain'] >= 1) & (df_clust['region.ElonginB_ElonginC_binding'].isna())])

        clust_alpha.loc[cluster, "ElonginB_ElonginC_binding"] = len(df_clust[(df_clust['region.⍺-Domain'] >= 1) & (df_clust['region.ElonginB_ElonginC_binding'] >= 1)])

        clust_beta.loc[cluster, "β-Domain"] = len(df_clust[(df_clust['region.β-Domain'] >= 1) & (df_clust['region.HIF1_alpha_binding'].isna())])

        clust_beta.loc[cluster, "HIF1_alpha_binding"] = len(df_clust[(df_clust['region.β-Domain'] >= 1) & (df_clust['region.HIF1_alpha_binding'] >= 1)])

        clust_cds.loc[cluster, "Outside of ⍺-Domain and β-Domain"] = len(df_clust[(df_clust['region.Outside of ⍺-Domain and β-Domain'] >= 1) & (df_clust['region.GXEEX8'].isna())])

        clust_cds.loc[cluster, "GXEEX8"] = len(df_clust[(df_clust['region.Outside of ⍺-Domain and β-Domain'] >= 1) & (df_clust['region.GXEEX8'] >= 1)])

    all_phen_dfs = [clust_alpha, clust_beta, clust_cds]
    for i in range(len(all_phen_dfs)):
        all_phen_dfs[i] = all_phen_dfs[i].rename_axis("Cluster", axis=0).rename_axis("Region", axis=1)


    plot_clustered_stacked(all_phen_dfs, ["⍺-Domain", 'β-Domain', 'Outside of ⍺-Domain and β-Domain'])

    plt.savefig(os.path.join(fig_path, f'clustered_regions.pdf'))
    plt.close('all')


    codon = df.set_index("codon_start")
    codon = pd.get_dummies(codon[codon.index.notnull()]["cluster_labels"]).sort_index()
    codon = codon.groupby(codon.index).sum()

    axs = codon.plot(kind="bar", xticks=[], subplots=True, figsize=(12, 10), title=["" for v in codon.columns ])
    # for ax in axs:
    #     ax.set_xticks(DOMAIN_TICKS)
    plt.savefig(os.path.join(fig_path, f'clustered_codon_start.pdf'))
    plt.close('all')

    return df


