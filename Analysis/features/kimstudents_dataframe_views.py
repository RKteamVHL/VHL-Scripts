import os
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler

from .kimstudents_dataframe_preprocessing import COMPUTED_COLUMNS
from .. import variant_functions as vf
from .kimstudents_dataframe_clustering import *

STATS_DIR = "statistics"
if not os.path.isdir(STATS_DIR):
    os.makedirs(STATS_DIR)

FIGURE_DIR = "figures"
DATA_DIR = "data"


TESTS_DIR = "tests"

DOMAIN_TICKS = [1, 62, 154, 192, 204, 213]

COLOR20 = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
           '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']

DOMAIN_COLORS = ['#fa0707', '#f56c6c', '#0000f7', '#6565fc', '#808080', '#c9c9c9']

mpl.rcParams['axes.prop_cycle'] = cycler(color=COLOR20)
# TAKE_TOP = -1


def regions_alpha_beta(df):
    ab_cols = ['region.⍺-Domain', 'region.β-Domain', 'region.Outside of ⍺-Domain and β-Domain']
    df = df[[*COMPUTED_COLUMNS["generalized_phenotype"], *ab_cols]]
    pheno_domains = pd.DataFrame(columns=ab_cols, index=COMPUTED_COLUMNS["generalized_phenotype"])
    for col in COMPUTED_COLUMNS["generalized_phenotype"]:
        phen_agg = df[df[col] >= 1].sum()
        pheno_domains.loc[col] = phen_agg[ab_cols]

    pheno_domains = pheno_domains.rename_axis("Phenotype", axis=0).rename_axis("Region", axis=1).rename(
        lambda x: x.split(".")[1])
    pheno_domains = pheno_domains.sort_index()

    ax = pheno_domains.plot(kind="bar", legend=True, figsize=(10, 6), color=DOMAIN_COLORS[0::2])
    ax.figure.subplots_adjust(left=0.1, bottom=0.4)
    plt.ylabel("Number of Occurences")

    return pheno_domains


def regions_elongin_hifa(df):
    ehif_cols = ['region.ElonginB_ElonginC_binding', 'region.HIF1_alpha_binding', 'region.GXEEX8']
    df = df[[*COMPUTED_COLUMNS["generalized_phenotype"], *ehif_cols]]
    pheno_domains = pd.DataFrame(columns=ehif_cols, index=COMPUTED_COLUMNS["generalized_phenotype"])
    for col in COMPUTED_COLUMNS["generalized_phenotype"]:
        phen_agg = df[df[col] >= 1].sum()
        pheno_domains.loc[col] = phen_agg[ehif_cols]

    pheno_domains = pheno_domains.rename_axis("Phenotype", axis=0).rename_axis("Region", axis=1).rename(
        lambda x: x.split(".")[1])
    pheno_domains = pheno_domains.sort_index()

    ax = pheno_domains.plot(kind="bar", legend=True, figsize=(10, 6), color=DOMAIN_COLORS[1::2])
    ax.figure.subplots_adjust(left=0.1, bottom=0.4)
    plt.ylabel("Number of Occurences")
    return pheno_domains


def missense_regions_alpha_beta(df):
    df = df.dropna(subset=['generalized_mutant_type.missense_variant'])

    return regions_alpha_beta(df)


def missense_regions(df):
    df = df.dropna(subset=['generalized_mutant_type.missense_variant'])
    return regions(df)


def regions(df):
    df = df[[*COMPUTED_COLUMNS["generalized_phenotype"], *COMPUTED_COLUMNS["region"]]]

    alpha_cols = ["ElonginB_ElonginC_binding", "⍺-Domain"]
    beta_cols = ["HIF1_alpha_binding", "β-Domain"]
    cds_cols = ["GXEEX8", "Outside of ⍺-Domain and β-Domain"]

    pheno_alpha = pd.DataFrame(columns=alpha_cols, index=COMPUTED_COLUMNS["generalized_phenotype"])
    pheno_beta = pd.DataFrame(columns=beta_cols, index=COMPUTED_COLUMNS["generalized_phenotype"])
    pheno_cds = pd.DataFrame(columns=cds_cols, index=COMPUTED_COLUMNS["generalized_phenotype"])

    for row in COMPUTED_COLUMNS["generalized_phenotype"]:
        df_phen = df[df[row] >= 1]

        pheno_alpha.loc[row, "⍺-Domain"] = len(
            df_phen[(df_phen['region.⍺-Domain'] >= 1) & (df_phen['region.ElonginB_ElonginC_binding'].isna())])

        pheno_alpha.loc[row, "ElonginB_ElonginC_binding"] = len(
            df_phen[(df_phen['region.⍺-Domain'] >= 1) & (df_phen['region.ElonginB_ElonginC_binding'] >= 1)])

        pheno_beta.loc[row, "β-Domain"] = len(
            df_phen[(df_phen['region.β-Domain'] >= 1) & (df_phen['region.HIF1_alpha_binding'].isna())])

        pheno_beta.loc[row, "HIF1_alpha_binding"] = len(
            df_phen[(df_phen['region.β-Domain'] >= 1) & (df_phen['region.HIF1_alpha_binding'] >= 1)])

        pheno_cds.loc[row, "Outside of ⍺-Domain and β-Domain"] = len(
            df_phen[(df_phen['region.Outside of ⍺-Domain and β-Domain'] >= 1) & (df_phen['region.GXEEX8'].isna())])

        pheno_cds.loc[row, "GXEEX8"] = len(
            df_phen[(df_phen['region.Outside of ⍺-Domain and β-Domain'] >= 1) & (df_phen['region.GXEEX8'] >= 1)])

    all_phen_dfs = [pheno_alpha, pheno_beta, pheno_cds]
    for i in range(len(all_phen_dfs)):
        all_phen_dfs[i] = all_phen_dfs[i].rename_axis("Phenotype", axis=0).rename_axis("Region", axis=1).rename(
            lambda x: x.split(".")[1])

    combined_df = pheno_alpha.join(pheno_beta).join(pheno_cds)

    order = combined_df.sum().sort_values(ascending=False).index
    combined_df = combined_df[order]

    plot_clustered_stacked(all_phen_dfs, ["⍺-Domain", 'β-Domain', 'Outside of ⍺-Domain and β-Domain'])
    return combined_df


def domains_adjusted(df):
    df = df[[*COMPUTED_COLUMNS["generalized_phenotype"], *COMPUTED_COLUMNS["domain"]]]
    pheno_domains = pd.DataFrame(columns=COMPUTED_COLUMNS["domain"], index=COMPUTED_COLUMNS["generalized_phenotype"])

    for col in COMPUTED_COLUMNS["generalized_phenotype"]:
        phen_agg = df[df[col] >= 1].sum()
        pheno_domains.loc[col] = phen_agg[COMPUTED_COLUMNS["domain"]]

    pheno_domains = pheno_domains.rename_axis("Phenotype", axis=0).rename_axis("Domain", axis=1).rename(
        lambda x: x.split(".")[1])
    pheno_domains = pheno_domains.sort_index()

    pheno_domains["region.β-Domain"] = pheno_domains["region.β-Domain"] * (vf.CDS_LEN / vf.BETA_LEN)
    pheno_domains["region.⍺-Domain"] = pheno_domains["region.⍺-Domain"] * (vf.CDS_LEN / vf.ALPHA_LEN)
    pheno_domains["region.Outside of ⍺-Domain and β-Domain"] = pheno_domains[
                                                                   "region.Outside of ⍺-Domain and β-Domain"] * (
                                                                           vf.CDS_LEN / (
                                                                               1 - vf.ALPHA_LEN + vf.BETA_LEN))

    ax = pheno_domains.plot(kind="bar", legend=True, figsize=(10, 6))
    ax.figure.subplots_adjust(left=0.1, bottom=0.4)
    plt.ylabel("Number of Occurences, adjusted by domain lengths")

    return pheno_domains


def mutant_type_counts(df):
    df = df[[*COMPUTED_COLUMNS["generalized_phenotype"], *COMPUTED_COLUMNS["generalized_mutant_type"]]]
    pheno_muttypes = pd.DataFrame(columns=COMPUTED_COLUMNS["generalized_mutant_type"],
                                  index=COMPUTED_COLUMNS["generalized_phenotype"])

    for col in COMPUTED_COLUMNS["generalized_phenotype"]:
        phen_agg = df[df[col] >= 1].sum()
        pheno_muttypes.loc[col] = phen_agg[COMPUTED_COLUMNS["generalized_mutant_type"]]

    pheno_muttypes = pheno_muttypes.rename_axis("Phenotype", axis=0).rename_axis("Mutant Type", axis=1).rename(
        lambda x: x.split(".")[1])
    pheno_muttypes = pheno_muttypes.sort_index().rename(columns=lambda x: x.split(".")[1])

    ax = pheno_muttypes.plot(kind="bar", legend=True, figsize=(12, 6), fontsize=8)
    ax.legend(ncol=2)
    ax.figure.subplots_adjust(left=0.1, bottom=0.35)

    plt.ylabel("Number of Occurences")

    return pheno_muttypes


def mutant_type_ratios(df):
    df = df[[*COMPUTED_COLUMNS["generalized_phenotype"], *COMPUTED_COLUMNS["generalized_mutant_type"]]]
    pheno_muttypes = pd.DataFrame(columns=COMPUTED_COLUMNS["generalized_mutant_type"],
                                  index=COMPUTED_COLUMNS["generalized_phenotype"])

    for col in COMPUTED_COLUMNS["generalized_phenotype"]:
        phen_agg = df[df[col] >= 1].sum()
        pheno_muttypes.loc[col] = phen_agg[COMPUTED_COLUMNS["generalized_mutant_type"]]

    pheno_muttypes = pheno_muttypes.rename_axis("Phenotype", axis=0).rename_axis("Mutant Type", axis=1).rename(
        lambda x: x.split(".")[1])
    pheno_muttypes = pheno_muttypes.sort_index().rename(columns=lambda x: x.split(".")[1])

    pheno_muttypes = pheno_muttypes.loc[:, :].div(pheno_muttypes.sum(axis=1), axis=0)

    ax = pheno_muttypes.plot(kind="bar", legend=True, figsize=(12, 6), fontsize=8)
    ax.legend(ncol=2)
    ax.figure.subplots_adjust(left=0.1, bottom=0.35)

    plt.ylabel("Number of Occurences")

    return pheno_muttypes


def grouped_mutant_type_ratios(df):
    df = df[[*COMPUTED_COLUMNS["generalized_phenotype"], *COMPUTED_COLUMNS["grouped_mutation_type"]]]
    pheno_muttypes = pd.DataFrame(columns=COMPUTED_COLUMNS["grouped_mutation_type"],
                                  index=COMPUTED_COLUMNS["generalized_phenotype"])

    for col in COMPUTED_COLUMNS["generalized_phenotype"]:
        phen_agg = df[df[col] >= 1].sum()
        pheno_muttypes.loc[col] = phen_agg[COMPUTED_COLUMNS["grouped_mutation_type"]]

    pheno_muttypes = pheno_muttypes.rename_axis("Phenotype", axis=0).rename_axis("Mutant Type", axis=1).rename(
        lambda x: x.split(".")[1])
    pheno_muttypes = pheno_muttypes.sort_index().rename(columns=lambda x: x.split(".")[1])

    pheno_muttypes = pheno_muttypes.loc[:, :].div(pheno_muttypes.sum(axis=1), axis=0)

    ax = pheno_muttypes.plot(kind="bar", legend=True, figsize=(12, 6), fontsize=8)
    ax.legend(ncol=2)
    ax.figure.subplots_adjust(left=0.1, bottom=0.35)

    autolabel(ax, ax.patches)

    plt.ylabel("Ratios of Mutation Types")

    return pheno_muttypes


def grouped_mutant_type_counts(df):
    df = df[[*COMPUTED_COLUMNS["generalized_phenotype"], *COMPUTED_COLUMNS["grouped_mutation_type"]]]
    pheno_muttypes = pd.DataFrame(columns=COMPUTED_COLUMNS["grouped_mutation_type"],
                                  index=COMPUTED_COLUMNS["generalized_phenotype"])

    for col in COMPUTED_COLUMNS["generalized_phenotype"]:
        phen_agg = df[df[col] >= 1].sum()
        pheno_muttypes.loc[col] = phen_agg[COMPUTED_COLUMNS["grouped_mutation_type"]]

    pheno_muttypes = pheno_muttypes.rename_axis("Phenotype", axis=0).rename_axis("Mutant Type", axis=1).rename(
        lambda x: x.split(".")[1])
    pheno_muttypes = pheno_muttypes.sort_index().rename(columns=lambda x: x.split(".")[1])

    # pheno_muttypes =  pheno_muttypes.loc[:, :].div(pheno_muttypes.sum(axis=1), axis=0)

    ax = pheno_muttypes.plot(kind="bar", legend=True, figsize=(12, 6), fontsize=8)
    ax.legend(ncol=2)
    ax.figure.subplots_adjust(left=0.1, bottom=0.35)

    autolabel(ax, ax.patches)

    plt.ylabel("Ratios of Mutation Types")

    return pheno_muttypes


def codon_phenotype_subplots(df):
    df = df.set_index("codon_start")
    df = df[df.index.notnull()]

    phens = df[COMPUTED_COLUMNS["generalized_phenotype"]]

    codon = phens.groupby(phens.index).agg("sum")
    codon = codon.rename_axis("Codon Position").rename(columns=lambda x: x.split(".")[1])

    combined = codon.sum(axis=1)

    axs = codon.plot(kind="bar", xticks=[], subplots=True, figsize=(12, 10), title=["" for v in codon.columns])

    return codon


def codon_histogram(df):
    df = df.set_index("codon_start")
    df = df[df.index.notnull()]

    phens = df[COMPUTED_COLUMNS["generalized_phenotype"]]

    codon = phens.groupby(phens.index).agg("sum")
    codon = codon.rename_axis("Codon Position").rename(columns=lambda x: x.split(".")[1])

    combined = codon.sum(axis=1)
    split = [combined.sort_values().iloc[-2], combined.sort_values().iloc[-1]]

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [1, 3]})
    fig.subplots_adjust(hspace=0.05)
    fig.set_figwidth(8)
    # adjust space between axes

    ax1.bar(combined.index.to_numpy(), combined.to_numpy())
    ax2.bar(combined.index.to_numpy(), combined.to_numpy())

    ax1.set_ylim(split[1] - 5, split[1] + 5)  # outliers only
    ax2.set_ylim(0, split[0] + 10)  # most of the data

    # hide the spines between ax and ax2
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.xaxis.tick_top()
    ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()

    d = .3  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                  linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

    ax1.set_xticks(DOMAIN_TICKS)
    ax2.set_xticks(DOMAIN_TICKS)

    for p in ax1.patches:
        thresh = 50

        if p.get_height() >= thresh:
            # ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
            plt.annotate(str(int(p.get_x() + p.get_width() / 2.)), (p.get_x() + p.get_width(), p.get_height()),
                         ha='left',
                         va='center', xytext=(-2, 5), textcoords='offset points', rotation=30)
    for p in ax2.patches:
        thresh = 50

        if p.get_height() >= thresh:
            # ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
            plt.annotate(str(int(p.get_x() + p.get_width() / 2.)), (p.get_x() + p.get_width(), p.get_height()),
                         ha='left',
                         va='center', xytext=(-2, 5), textcoords='offset points', rotation=30)

    plt.ylabel('# of Mutations')
    plt.xlabel('Codon Position')

    return combined


def ratio_of_phenotypes(df, generalized=True):
    type_ = "phenotype"
    if generalized:
        type_ = "generalized_" + type_
    phens = df[COMPUTED_COLUMNS[type_]].fillna(0)
    phens_counts = phens.agg("sum").sort_values()

    phens_ratio = phens_counts.rename_axis("Phenotype").rename("Share of Observations").rename(
        lambda x: x.split(".")[1])

    fig = plt.figure(figsize=(8, 6))
    fig.subplots_adjust(bottom=0.4)
    ax = phens_ratio.plot.bar()
    ax.set_xlabel(phens_ratio.axes[0].name)
    ax.set_ylabel(phens_ratio.name)
    for p in ax.patches:
        # ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
        ax.annotate(str(int(p.get_height())), (p.get_x() + p.get_width() / 2., p.get_height()), ha='center',
                    va='center', xytext=(0, 10), textcoords='offset points')

    # plt.close(fig)
    #
    # phens_ratio = phens_ratio / phens_counts.sum()
    #
    # fig = plt.figure(figsize=(8, 6))
    # fig.subplots_adjust(bottom=0.4)
    # ax = phens_ratio.plot.bar()
    # ax.set_xlabel(phens_ratio.axes[0].name)
    # ax.set_ylabel(phens_ratio.name)
    # for p in ax.patches:
    #     # ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    #     ax.annotate(np.round(p.get_height(),decimals=2), (p.get_x() + p.get_width() / 2., p.get_height()), ha='center',
    #             va='center', xytext=(0, 10), textcoords='offset points')
    # plt.savefig(os.path.join(STATS_DIR, f'{type_}ratio.pdf'))
    # plt.close(fig)
    return phens_ratio


def penetrance(df):
    # this threshold is arbitrary right now
    # threshold = 10
    df = df.dropna(subset=COMPUTED_COLUMNS["age"])

    summed_phenos = df[COMPUTED_COLUMNS["generalized_phenotype"]].sum(axis=1)
    iso_ind = summed_phenos[summed_phenos == 1].index

    df = df.loc[iso_ind]
    # df["evaluated_age"] = df["evaluated_age"].apply(np.floor)
    df = df.set_index("evaluated_age").sort_index()
    iso_sorted = df[COMPUTED_COLUMNS["generalized_phenotype"]].fillna(0)
    iso_sorted.index = iso_sorted.index.astype(int)
    # iso_sorted = iso_sorted.groupby(iso_sorted.index).sum()
    # iso_sorted = iso_sorted.reindex(labels= list(range(iso_sorted.index.min(), iso_sorted.index.max()+1)), fill_value=0.0)
    # summed = iso_sorted.sum()
    # meets_thresh = summed[summed>threshold]
    # iso_sorted = iso_sorted[meets_thresh.index]

    pdf = iso_sorted / iso_sorted.sum()
    cdf = pdf.cumsum()
    cdf = cdf.rename(columns=lambda x: x.split(".")[1])

    fig = plt.figure(figsize=(8, 6))

    ax = plt.step(cdf.index.to_numpy(), cdf.to_numpy())

    plt.xlabel("Age (Years)")
    plt.ylabel("Cumulative Distribution")
    plt.legend(ax, cdf.columns)

    return iso_sorted


def _phenotype_correlation(df):
    pheno_pheno = pd.DataFrame(columns=COMPUTED_COLUMNS["generalized_phenotype"],
                               index=COMPUTED_COLUMNS["generalized_phenotype"])

    for col in COMPUTED_COLUMNS["generalized_phenotype"]:
        phen_agg = df[df[col] >= 1].sum()
        pheno_pheno.loc[col] = phen_agg

    pheno_pheno = pheno_pheno.reindex(pheno_pheno.sum().sort_values(ascending=False).index)
    pheno_pheno = pheno_pheno[pheno_pheno.index.to_list()]

    pheno_pheno = pheno_pheno.rename(index=lambda x: x.split(".")[1]).rename(columns=lambda x: x.split(".")[1])

    return pheno_pheno


def phenotype_correlation_counts(df):
    pheno_pheno = _phenotype_correlation(df)
    ax = sns.heatmap(pheno_pheno.astype(float))
    return pheno_pheno


def phenotype_correlation_ratio(df):
    pheno_pheno = _phenotype_correlation(df)
    # normalizing so that diagonals = 1; could probably do better normalization
    for col_i in pheno_pheno.columns:
        for col_j in pheno_pheno.columns:
            col_max = max(pheno_pheno[col_j].max(), pheno_pheno.loc[col_i].max())
            pheno_pheno.loc[col_i][col_j] = pheno_pheno.loc[col_i][col_j] / col_max

    ax = sns.heatmap(pheno_pheno.astype(float))

    return pheno_pheno


def phenotype_codon_heatmap(df):
    codon_df = df[df["codon_start"] >= 1]
    codon_df = codon_df.dropna(subset=['generalized_mutant_type.missense_variant'])
    codon_df = codon_df[[*COMPUTED_COLUMNS["generalized_phenotype"], *COMPUTED_COLUMNS["codon"]]]
    codon_df = codon_df.groupby("codon_start").sum()
    codon_df = codon_df.transpose()

    order = codon_df.sum().sort_values(ascending=False).index
    codon_df = codon_df[order]
    # codon_df = codon_df.iloc[:, 0:TAKE_TOP]

    ax = sns.heatmap(codon_df.astype(float))

    return codon_df


def phenotype_aachange_heatmap(df):
    missense_df = df[df["codon_start"] >= 1]
    missense_df = missense_df.dropna(subset=['generalized_mutant_type.missense_variant'])

    missense_df = missense_df[[*COMPUTED_COLUMNS["generalized_phenotype"], *COMPUTED_COLUMNS["aa_change"]]]
    aachange_df = pd.DataFrame(index=COMPUTED_COLUMNS["generalized_phenotype"], columns=COMPUTED_COLUMNS["aa_change"])

    for phen in aachange_df.index:
        aachange_df.loc[phen] = missense_df.dropna(subset=[phen]).sum()[COMPUTED_COLUMNS["aa_change"]]

    order = aachange_df.sum().sort_values(ascending=False).index
    aachange_df = aachange_df[order]
    # aachange_df = aachange_df.iloc[:, 0:TAKE_TOP]

    # aachange_df = aachange_df.loc[:, (aachange_df.sum() >= 5)]
    ax = sns.heatmap(aachange_df.astype(float))

    return aachange_df


def plot_clustered_stacked(dfall, labels=None, title="multiple stacked bar plot", H="/", **kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot.
    labels is a list of the names of the dataframe, used for the legend
    title is a string for the title of the plot
    H is the hatch used for identification of the different dataframe"""

    n_df = len(dfall)
    n_col = len(dfall[0].columns)
    n_ind = len(dfall[0].index)
    axes = plt.subplot(111)
    axes.figure.set_figwidth(8)
    for i in range(len(dfall)):  # for each data frame
        axe = dfall[i].plot(kind="bar",
                            linewidth=0,
                            stacked=True,
                            ax=axes,
                            legend=False,
                            grid=False,
                            color=DOMAIN_COLORS[i * n_col:],
                            **kwargs)

        # make bar plots

    h, l = axe.get_legend_handles_labels()  # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col):  # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i + n_col]):
            for rect in pa.patches:  # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(dfall[-1].index, rotation=0)
    axe.set_ylabel("Number of Observations")
    # axe.set_title(title)

    l1 = axe.legend(h, l, loc=[0.53, 0.55])
    # if labels is not None:
    #     l2 = plt.legend(n, labels, loc=[0.7, 0.6])
    axe.add_artist(l1)
    return axe


def autolabel(ax, rects, xpos='center', as_percentage=False):
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
        text = None
        if as_percentage:
            text = f'{height:.0f}' if float(int(height)) == height else f'{height:.1%}'
        else:
            text = f'{height:.0f}' if float(int(height)) == height else f'{height:.2f}'
        ax.text(rect.get_x() + rect.get_width() * offset[xpos], 1.01 * height,
                text, ha=ha[xpos], va='bottom')


def plot_cluster_property(df, figure_path, property_name, cluster_column="cluster_labels", use_mean=False,
                          save_csv=False, ratio_type='ratio_of_total'):
    df_sums = df.sum().sort_values()
    if not os.path.isdir(figure_path):
        os.makedirs(figure_path)

    if use_mean:
        df_means = df.groupby(cluster_column).mean()
        df_means = df_means[df_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1] if "." in x else x)
        df_stds = df.groupby(cluster_column).std()
        df_stds = df_stds[df_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1] if "." in x else x)
        df_means.plot(kind="bar", yerr=df_stds, figsize=(12, 8))
        # plt.savefig(os.path.join(figure_path, f'clustered_{property_name}_means.pdf'))
        plt.savefig(os.path.join(figure_path, f'clustered_{property_name}_means.eps'), format='eps')
        plt.close('all')

    else:
        df_counts = df.groupby(cluster_column).sum()
        df_counts = df_counts[df_sums.index.to_list()].rename(columns=lambda x: x.split(".")[1] if "." in x else x)
        ax = df_counts.plot(kind="bar", figsize=(12, 8))
        if save_csv:
            df_counts.to_csv(os.path.join(figure_path, f'clustered_{property_name}_counts.csv'))

        ax.set_xlabel("Cluster", rotation=0)
        ax.set_ylabel("Number of Observations")
        autolabel(ax, ax.patches)
        # plt.savefig(os.path.join(figure_path, f'clustered_{property_name}_counts.pdf'))
        plt.savefig(os.path.join(figure_path, f'clustered_{property_name}_counts.eps'), format='eps')
        plt.close('all')

        if ratio_type == "within":
            df_sums = df_counts.sum(axis=1)
            df_ratio = df_counts.divide(df_sums, axis='index')
        elif ratio_type == "across":
            df_sums = df_counts.sum()
            df_ratio = df_counts.divide(df_sums)
        elif ratio_type == "ratio_of_total":
            df_sums = df.groupby(cluster_column).count()
            df_sums = df_sums.rename(columns=lambda x: x.split(".")[1] if "." in x else x)
            df_ratio = df_counts.divide(df_sums, axis='index')
            df_ratio = df_ratio[df_counts.columns]
        ax = df_ratio.plot(kind="bar", figsize=(12, 8))
        if save_csv:
            df_ratio.to_csv(os.path.join(figure_path, f'clustered_{property_name}_ratios.csv'))
        ax.set_xlabel("Cluster", rotation=0)
        ax.set_ylabel("Ratio of Observations")
        autolabel(ax, ax.patches, as_percentage=True)

        # plt.savefig(os.path.join(figure_path, f'clustered_{property_name}_ratios.pdf'))
        plt.savefig(os.path.join(figure_path, f'clustered_{property_name}_ratios.eps'), format='eps')
        plt.close('all')


def create_descriptive_figures(dfs):
    for df_type, df_out in dfs.items():
        fns = [
            regions_alpha_beta,
            regions_elongin_hifa,
            regions,
            # missense_domains,
            mutant_type_counts,
            mutant_type_ratios,
            codon_phenotype_subplots,
            codon_histogram,
            ratio_of_phenotypes,
            phenotype_correlation_counts,
            phenotype_correlation_ratio,
            penetrance,
            grouped_mutant_type_ratios,
            grouped_mutant_type_counts,
            phenotype_codon_heatmap,
            phenotype_aachange_heatmap
        ]

        for fn in fns:
            stats_name = fn.__name__
            plt.close('all')
            dataframe = fn(df_out)
            stats_path = os.path.join(STATS_DIR, df_type)
            fig_path = os.path.join(stats_path, FIGURE_DIR)
            data_path = os.path.join(stats_path, DATA_DIR)
            if not os.path.isdir(stats_path):
                os.makedirs(stats_path)

            if not os.path.isdir(fig_path):
                os.makedirs(fig_path)

            if not os.path.isdir(data_path):
                os.makedirs(data_path)

            # plt.savefig(os.path.join(fig_path, f'{stats_name}.pdf'))
            plt.savefig(os.path.join(fig_path, f'{stats_name}.eps'), format='eps')
            dataframe.to_csv(os.path.join(data_path, f'{stats_name}.csv'))

def create_cluster_figures(dfs):
    for df_type, df_out in dfs.items():

        clustered = dataframe_snf(df_out).fillna(0)

        create_cluster_summaries(clustered, df_type)
        create_cluster_phenotype_summaries(clustered, df_type)

def create_cluster_summaries(df, analysis_type):
    base_path = os.path.join(STATS_DIR, analysis_type, "cluster")
    if not os.path.isdir(base_path):
        os.makedirs(base_path)

    for clust_type in ["cluster_labels_best", "cluster_labels_second"]:
        fig_path = os.path.join(base_path, clust_type)
        if not os.path.isdir(fig_path):
            os.makedirs(fig_path)

        properties = ["generalized_phenotype", "grouped_mutation_type", "domain", "aa_change"]
        for prop in properties:
            prop_df = df.set_index(clust_type)[COMPUTED_COLUMNS[prop]].fillna(0)
            plot_cluster_property(prop_df, cluster_column=clust_type, figure_path=fig_path, property_name=prop,
                                  save_csv=True)

        prop = "age"
        prop_df = df.set_index(clust_type)[COMPUTED_COLUMNS[prop]]
        plot_cluster_property(prop_df, cluster_column=clust_type, figure_path=fig_path, property_name=prop,
                              use_mean=True)
        df["codon_start"] = df["codon_start"].astype(int)
        codon = df.set_index("codon_start")


        codon = codon[codon['generalized_mutant_type.missense_variant'] != 0]
        # codon = codon.dropna(subset=['generalized_mutant_type.missense_variant'])
        codon = pd.get_dummies(codon[codon.index.notnull()][clust_type]).sort_index()
        codon = codon.groupby(codon.index).sum()
        codons = pd.DataFrame(index=set(range(0, 214)), columns=codon.columns)
        codons[:] = 0
        codons.loc[set(codon.index)] = codon
        axs = codons.plot(kind="bar", xticks=[], subplots=True, figsize=(12, 10), title=["" for v in codon.columns])
        # for ax in axs:
        #     ax.set_xticks(DOMAIN_TICKS)
        # plt.savefig(os.path.join(fig_path, f'clustered_codon_start.pdf'))
        plt.savefig(os.path.join(fig_path, f'clustered_codon_start.eps'), format='eps')

        df.to_csv(os.path.join(fig_path, f"clustered_out.tsv"), sep='\t')
        plt.close('all')


def create_cluster_phenotype_summaries(df, analysis_type):
    base_path = os.path.join(STATS_DIR, analysis_type, "cluster")
    for clust_type in ["cluster_labels_best", "cluster_labels_second"]:
        fig_path = os.path.join(base_path, clust_type)
        if not os.path.isdir(fig_path):
            os.makedirs(fig_path)

        df = df.set_index(clust_type)
        for phenotype in COMPUTED_COLUMNS["generalized_phenotype"]:
            fig_path_pheno = os.path.join(fig_path, phenotype)

            df_phen = df[df[phenotype] >= 1]
            properties = ["grouped_mutation_type", "domain", "aa_change"]
            for prop in properties:
                prop_df = df_phen[COMPUTED_COLUMNS[prop]].fillna(0)
                plot_cluster_property(prop_df, cluster_column=clust_type, figure_path=fig_path_pheno,
                                      property_name=prop, ratio_type="across")

            prop_df = df_phen[phenotype].fillna(0)
            df_counts = prop_df.groupby(clust_type).sum()
            df_ratio = df_counts.divide(df_counts.sum())
            ax = df_ratio.plot(kind="bar", figsize=(12, 8))
            autolabel(ax, ax.patches)
            plt.savefig(os.path.join(fig_path, f'clustered_phenotype_ratios.eps'), format='eps')
            # plt.savefig(os.path.join(fig_path, f'clustered_phenotype_ratios.pdf'))
            plt.close('all')
