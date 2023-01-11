import os
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler

from .. import variant_functions as vf
from .kimstudents_dataframe_clustering import *


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
    df = df.sort_index()

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


def _plot_codon(combined):
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
    ax2.set_xticklabels(ax2.get_xticks(), rotation=45)

    for p in ax1.patches:
        thresh = 17

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

def codon_histogram(df):
    df = df.set_index("codon_start")
    df = df[df.index.notnull()]
    df["default_score"] = 1

    codon_score = df["default_score"]

    codon = codon_score.groupby(codon_score.index).agg("sum")
    codon = codon.rename("Codon Count")

    _plot_codon(codon)
    return codon


def codon_blosum62_histogram(df):
    df = df.set_index("codon_start")
    df = df[df.index.notnull()]

    codon_blosum = df["blosum62_score"]

    codon = codon_blosum.groupby(codon_blosum.index).agg("sum")
    codon = codon.rename("BLOSUM Score")


    _plot_codon(codon)
    return codon

def codon_blosum90_histogram(df):
    df = df.set_index("codon_start")
    df = df[df.index.notnull()]

    codon_blosum = df["blosum90_score"]

    codon = codon_blosum.groupby(codon_blosum.index).agg("sum")
    codon = codon.rename("BLOSUM Score")

    _plot_codon(codon)
    return codon

def penetrance(df):

    df = df.dropna(subset=COMPUTED_COLUMNS["age"])

    summed_phenos = df[COMPUTED_COLUMNS["generalized_phenotype"]].sum(axis=1)
    iso_ind = summed_phenos[summed_phenos == 1].index

    df = df.loc[iso_ind]

    df = df.set_index("evaluated_age").sort_index()
    iso_df = df[COMPUTED_COLUMNS["generalized_phenotype"]].fillna(0).rename(columns=lambda x: x.split(".")[1])
    iso_df.index = iso_df.index.astype(int)

    sorted = iso_df.sum().sort_values(ascending=False)
    sorted_filtered = sorted[sorted != 0]
    iso_sorted = iso_df[sorted_filtered.index]
    pdf = iso_sorted / iso_sorted.sum()
    cdf = pdf.cumsum()

    fig = plt.figure(figsize=(8, 6))

    ax = plt.step(cdf.index.to_numpy(), cdf.to_numpy())

    plt.xlabel("Age (Years)")
    plt.ylabel("Cumulative Distribution")
    plt.legend(ax, [f"{ind} (N={int(val)})" for ind, val in sorted_filtered.iteritems()])

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

def create_descriptive_figures(directory, dfs):
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
            codon_blosum62_histogram,
            codon_blosum90_histogram,
            phenotype_correlation_counts,
            phenotype_correlation_ratio,
            penetrance,
            grouped_mutant_type_ratios,
            grouped_mutant_type_counts,
        ]

        for fn in fns:
            stats_name = fn.__name__
            plt.close('all')
            dataframe = fn(df_out)
            stats_path = os.path.join(directory, df_type)
            fig_path = os.path.join(stats_path, "figures")
            data_path = os.path.join(stats_path, "data")
            if not os.path.isdir(stats_path):
                os.makedirs(stats_path)

            if not os.path.isdir(fig_path):
                os.makedirs(fig_path)

            if not os.path.isdir(data_path):
                os.makedirs(data_path)

            # plt.savefig(os.path.join(fig_path, f'{stats_name}.pdf'))
            plt.savefig(os.path.join(fig_path, f'{stats_name}.eps'), format='eps')
            dataframe.to_csv(os.path.join(data_path, f'{stats_name}.csv'))