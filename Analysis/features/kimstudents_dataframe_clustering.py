from .kimstudents_dataframe_preprocessing import COMPUTED_COLUMNS
from . kimstudents_dataframe_views import plot_clustered_stacked
import snf
import os
import pandas as pd
import matplotlib.pyplot as plt
from . import constants
import seaborn as sns
import numpy as np

from .. import variant_functions as vf
# from sklearn.cluster import spectral_clustering
from sklearn.cluster import SpectralClustering


DOMAIN_TICKS = [1, 62, 154, 192, 204, 213]

def dataframe_snf(df, analysis_type):
    fig_path = os.path.join(constants.FIGURE_DIR, analysis_type)
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

    #converting to bools:
    df_cols = df_cols >= 1

    feat_metrics = []
    for col_group in col_groups:
        df_group = df_cols[col_group]

        # affinity = snf.make_affinity(df_group.to_numpy(), metric='sqeuclidean', K=len(df_cols.index), mu=0.5)
        affinity = snf.make_affinity(df_group.to_numpy(), normalize=False,  K=len(df_cols.index), metric='jaccard', mu=0.7)
        feat_metrics.append(affinity)

        pass

    fused_net = None
    if len(feat_metrics) == 1:
        fused_net =feat_metrics[0]
    else:
        fused_net = snf.snf(feat_metrics, K=50)
        # fused_net = snf.snf(feat_metrics, K=len(df_cols.index))




    best, second = snf.get_n_clusters(fused_net)
    n_clust = best

    sc = SpectralClustering(n_clust, affinity='precomputed', n_init=100, assign_labels='discretize')
    labels = sc.fit_predict(fused_net)
    df["cluster_labels"] = labels
    df = df.sort_values(by="cluster_labels")
    df.to_csv(os.path.join(fig_path, f"clustered_out.tsv"), sep='\t')


    plt.close('all')
    # fused_net = fused_net - np.diag(np.diag(fused_net))

    sorted_i = np.argsort(labels, axis=0)
    sorted_fused = fused_net[sorted_i, :]
    sorted_fused = sorted_fused[:, sorted_i]
    ax = plt.imshow(sorted_fused, cmap='hot', interpolation='nearest')
    plt.savefig(os.path.join(fig_path, f'sorted_affinities.pdf'))
    plt.close('all')

    properties = ["generalized_phenotype", "grouped_mutation_type"]
    for prop in properties:
        prop_df = df.set_index("cluster_labels")[COMPUTED_COLUMNS[prop]].fillna(0)
        plot_cluster_by_df(prop_df, figure_path=fig_path, property_name=prop)

    prop = "age"
    prop_df = df.set_index("cluster_labels")[COMPUTED_COLUMNS[prop]]
    plot_cluster_by_df(prop_df, figure_path=fig_path, property_name=prop, use_mean=True)

    # ix = list(df['cluster_labels'].unique())
    # alpha_cols = ["⍺-Domain", "ElonginB_ElonginC_binding"]
    # beta_cols = ["β-Domain", "HIF1_alpha_binding"]
    # cds_cols = ["Outside of ⍺-Domain and β-Domain", "GXEEX8"]
    # clust_alpha = pd.DataFrame(columns=alpha_cols, index=ix)
    # clust_beta = pd.DataFrame(columns=beta_cols, index=ix)
    # clust_cds = pd.DataFrame(columns=cds_cols, index=ix)
    #
    #
    # for cluster in ix:
    #     df_clust = df[df['cluster_labels'] == cluster]
    #
    #     clust_alpha.loc[cluster, "⍺-Domain"] = len(df_clust[(df_clust['region.⍺-Domain'] >= 1) & (df_clust['region.ElonginB_ElonginC_binding'].isna())])
    #
    #     clust_alpha.loc[cluster, "ElonginB_ElonginC_binding"] = len(df_clust[(df_clust['region.⍺-Domain'] >= 1) & (df_clust['region.ElonginB_ElonginC_binding'] >= 1)])
    #
    #     clust_beta.loc[cluster, "β-Domain"] = len(df_clust[(df_clust['region.β-Domain'] >= 1) & (df_clust['region.HIF1_alpha_binding'].isna())])
    #
    #     clust_beta.loc[cluster, "HIF1_alpha_binding"] = len(df_clust[(df_clust['region.β-Domain'] >= 1) & (df_clust['region.HIF1_alpha_binding'] >= 1)])
    #
    #     clust_cds.loc[cluster, "Outside of ⍺-Domain and β-Domain"] = len(df_clust[(df_clust['region.Outside of ⍺-Domain and β-Domain'] >= 1) & (df_clust['region.GXEEX8'].isna())])
    #
    #     clust_cds.loc[cluster, "GXEEX8"] = len(df_clust[(df_clust['region.Outside of ⍺-Domain and β-Domain'] >= 1) & (df_clust['region.GXEEX8'] >= 1)])
    #
    # all_phen_dfs = [clust_alpha, clust_beta, clust_cds]
    # for i in range(len(all_phen_dfs)):
    #     all_phen_dfs[i] = all_phen_dfs[i].rename_axis("Cluster", axis=0).rename_axis("Region", axis=1)
    # plot_clustered_stacked(all_phen_dfs, ["⍺-Domain", 'β-Domain', 'Outside of ⍺-Domain and β-Domain'])
    #
    # plt.savefig(os.path.join(fig_path, f'clustered_regions.pdf'))
    # plt.close('all')


    codon = df.set_index("codon_start")
    # codon = codon[codon['generalized_mutant_type.missense_variant'] >= 1]
    codon = codon.dropna(subset=['generalized_mutant_type.missense_variant'])
    codon = pd.get_dummies(codon[codon.index.notnull()]["cluster_labels"]).sort_index()
    codon = codon.groupby(codon.index).sum()

    axs = codon.plot(kind="bar", xticks=[], subplots=True, figsize=(12, 10), title=["" for v in codon.columns ])
    # for ax in axs:
    #     ax.set_xticks(DOMAIN_TICKS)
    plt.savefig(os.path.join(fig_path, f'clustered_codon_start.pdf'))
    plt.close('all')

    return df


def plot_cluster_by_df(df, figure_path, property_name, use_mean=False):


    df_sums = df.sum().sort_values()

    if use_mean:
        df_means = df.groupby("cluster_labels").mean()
        df_means = df_means[df_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1] if "." in x else x)
        df_stds = df.groupby("cluster_labels").std()
        df_stds = df_stds[df_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1] if "." in x else x)
        df_means.plot(kind="bar", yerr=df_stds, figsize=(12, 8))
        plt.savefig(os.path.join(figure_path, f'clustered_{property_name}_means.pdf'))
        plt.close('all')

    else:
        df_counts = df.groupby("cluster_labels").sum()
        df_counts = df_counts[df_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1] if "." in x else x)
        ax = df_counts.plot(kind="bar", figsize=(12, 8))

        plt.savefig(os.path.join(figure_path, f'clustered_{property_name}_counts.pdf'))
        plt.close('all')

        df_sums = df_counts.sum(axis=1)

        df_ratio = df_counts.divide( df_sums, axis='index')
        ax = df_ratio.plot(kind="bar", figsize=(12, 8))
        autolabel(ax, ax.patches)

        plt.savefig(os.path.join(figure_path, f'clustered_{property_name}_ratios.pdf'))
        plt.close('all')


def autolabel(ax, rects, xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    xpos = xpos.lower()  # normalize the case of the parameter
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                f'{height:.2f}', ha=ha[xpos], va='bottom')