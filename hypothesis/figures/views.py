import os
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cycler import cycler
from ..annotations.Annotation import AnnotationHeader

from .. import variant_functions as vf



DOMAIN_TICKS = [1, 62, 154, 192, 204, 213]

COLOR20 = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
           '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']

DOMAIN_COLORS = ['#fa0707', '#f56c6c', '#0000f7', '#6565fc', '#808080', '#c9c9c9']

mpl.rcParams['axes.prop_cycle'] = cycler(color=COLOR20)
# TAKE_TOP = -1


def _vhl_disease_colnames(df: pd.DataFrame, generalized=True):
    header = AnnotationHeader.GENERALIZED_VHL_DISEASE_ENTITY if generalized else AnnotationHeader.DISEASE_ENTITY_CLEAN
    return list(df.columns[df.columns.str.startswith(str(header))])


def _vhl_aop_colnames(df: pd.DataFrame, generalized=True):
    header = AnnotationHeader.GENERALIZED_VHL_AGE_OF_PRESENTATION if generalized else AnnotationHeader.AGE_OF_PRESENTATION_CLEAN
    return list(df.columns[df.columns.str.startswith(str(header))])


def _vhl_muttype_colnames(df: pd.DataFrame, generalized=True):
    header = AnnotationHeader.GENERALIZED_MUTATION_TYPE if generalized else AnnotationHeader.MUTATION_TYPE_CLEAN
    return list(df.columns[df.columns.str.startswith(str(header))])


def _vhl_grouped_muttype_colnames(df: pd.DataFrame):
    header = AnnotationHeader.GROUPED_MUTATION_TYPE
    return list(df.columns[df.columns.str.startswith(str(header))])


def penetrance(df: pd.DataFrame):
    out_df = df.copy()

    aop_cols = _vhl_aop_colnames(df)

    out_df = out_df[aop_cols].dropna(how='all')

    min_age = 0
    max_age = (out_df.max()).max()

    aop_df = pd.DataFrame(index=np.arange(min_age, max_age+1), dtype=int)
    for col in out_df.columns:
        counts = out_df[col].value_counts().astype(int)
        counts.name = col

        aop_df = pd.concat([aop_df, counts], axis=1)

    aop_df = aop_df.fillna(0)
    pdf_aop_df = aop_df / aop_df.sum()
    cdf_aop_df = pdf_aop_df.cumsum()

    fig = plt.figure(figsize=(8, 6))

    ax = plt.step(cdf_aop_df.index.to_numpy(), cdf_aop_df.to_numpy())

    plt.xlabel("Age (Years)")
    plt.ylabel("Cumulative Distribution")
    plt.legend(ax, [f"{ind.split('.')[-1]} (N={int(val)})" for ind, val in aop_df.sum().items()])

    return aop_df


def regions(df, columns):
    df = df[[*_vhl_disease_colnames(df), *columns]]
    pheno_domains = pd.DataFrame(columns=columns, index=_vhl_disease_colnames(df))
    for col in _vhl_disease_colnames(df):
        phen_agg = df[df[col] >= 1].sum()
        pheno_domains.loc[col] = phen_agg[columns]

    pheno_domains = pheno_domains.rename_axis("Phenotype", axis=0).rename_axis("Region", axis=1)
    pheno_domains = pheno_domains.rename(lambda x: x.split(".")[-1]).rename(lambda x: x.split(".")[-1], axis=1)
    pheno_domains = pheno_domains.sort_index()

    ax = pheno_domains.plot(kind="bar", legend=True, figsize=(10, 6), color=DOMAIN_COLORS[0::2])
    ax.figure.subplots_adjust(left=0.1, bottom=0.4)
    plt.ylabel("Number of Occurrences")

    return pheno_domains


def regions_alpha_beta(df):
    ab_cols = [f'{str(AnnotationHeader.FUNCTIONAL_REGION)}.{x}'
               for x in ['⍺-Domain', 'β-Domain', 'Outside of ⍺-Domain and β-Domain']]
    return regions(df, ab_cols)


def regions_elongin_hifa(df):
    ehif_cols = [f'{str(AnnotationHeader.FUNCTIONAL_REGION)}.{x}'
               for x in ['ElonginB_ElonginC_binding', 'HIF1_alpha_binding', 'GXEEX8']]
    return regions(df, ehif_cols)


def missense_regions_alpha_beta(df):
    df = df.dropna(subset=['generalized_mutant_type.missense_variant'])

    return regions_alpha_beta(df)


def missense_regions(df):
    df = df.dropna(subset=['generalized_mutant_type.missense_variant'])
    return regions(df)


def mutant_type(df, as_ratio=False):
    mutant_colnames = _vhl_muttype_colnames(df)
    disease_colnames = _vhl_disease_colnames(df)
    df = df[[*disease_colnames, *mutant_colnames]]
    pheno_muttypes = pd.DataFrame(columns=mutant_colnames,
                                  index=disease_colnames)

    for col in disease_colnames:
        phen_agg = df[df[col] >= 1].sum()
        pheno_muttypes.loc[col] = phen_agg[mutant_colnames]

    pheno_muttypes = pheno_muttypes.rename_axis("Phenotype", axis=0).rename_axis("Mutant Type", axis=1).rename(
        lambda x: x.split(".")[-1])
    pheno_muttypes = pheno_muttypes.sort_index().rename(columns=lambda x: x.split(".")[-1])

    if as_ratio:
        pheno_muttypes = pheno_muttypes.loc[:, :].div(pheno_muttypes.sum(axis=1), axis=0)

    ax = pheno_muttypes.plot(kind="bar", legend=True, figsize=(12, 6), fontsize=8)
    ax.legend(ncol=2)
    ax.figure.subplots_adjust(left=0.1, bottom=0.35)

    plt.ylabel("Number of Occurences")

    return pheno_muttypes


def mutant_type_counts(df):
    return mutant_type(df, as_ratio=False)


def mutant_type_ratios(df):
    return mutant_type(df, as_ratio=True)


def grouped_mutant_type(df, as_ratio=False):
    dis_colnames = _vhl_disease_colnames(df)
    mut_colnames = _vhl_grouped_muttype_colnames(df)
    df = df[[*dis_colnames, *mut_colnames]]
    pheno_muttypes = pd.DataFrame(columns=mut_colnames,
                                  index=dis_colnames)

    for col in dis_colnames:
        phen_agg = df[df[col] >= 1].sum()
        pheno_muttypes.loc[col] = phen_agg[mut_colnames]

    pheno_muttypes = pheno_muttypes.rename_axis("Phenotype", axis=0).rename_axis("Mutant Type", axis=1).rename(
        lambda x: x.split(".")[-1])
    pheno_muttypes = pheno_muttypes.sort_index().rename(columns=lambda x: x.split(".")[-1])

    if as_ratio:
        pheno_muttypes =  pheno_muttypes.loc[:, :].div(pheno_muttypes.sum(axis=1), axis=0)

    ax = pheno_muttypes.plot(kind="bar", legend=True, figsize=(12, 6), fontsize=8)
    ax.legend(ncol=2)
    ax.figure.subplots_adjust(left=0.1, bottom=0.35)

    autolabel(ax, ax.patches)

    plt.ylabel("Ratios of Mutation Types")

    return pheno_muttypes

def grouped_mutant_type_ratios(df):
    return grouped_mutant_type(df, as_ratio=True)


def grouped_mutant_type_counts(df):
    return grouped_mutant_type(df, as_ratio=False)


def codon_phenotype_subplots(df):
    df = df.set_index(str(AnnotationHeader.PROTEIN_POSITION_CLEAN))
    df = df[df.index.notnull()]
    pheno_colnames = _vhl_disease_colnames(df)
    phens = df[pheno_colnames]

    codon = phens.groupby(phens.index).agg("sum")
    codon = codon.rename_axis("Codon Position").rename(columns=lambda x: x.split(".")[-1])

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
    df = df.set_index(str(AnnotationHeader.PROTEIN_POSITION_CLEAN))
    df = df[df.index.notnull()]
    df["default_score"] = 1

    codon_score = df["default_score"]

    codon = codon_score.groupby(codon_score.index).agg("sum")
    codon = codon.rename("Codon Count")

    _plot_codon(codon)
    return codon


def codon_blosum62_histogram(df):
    df = df.set_index(str(AnnotationHeader.PROTEIN_POSITION_CLEAN))
    df = df[df.index.notnull()]

    codon_blosum = df["blosum62_score"]

    codon = codon_blosum.groupby(codon_blosum.index).agg("sum")
    codon = codon.rename("BLOSUM Score")


    _plot_codon(codon)
    return codon

def codon_blosum90_histogram(df):
    df = df.set_index(str(AnnotationHeader.PROTEIN_POSITION_CLEAN))
    df = df[df.index.notnull()]

    codon_blosum = df["blosum90_score"]

    codon = codon_blosum.groupby(codon_blosum.index).agg("sum")
    codon = codon.rename("BLOSUM Score")

    _plot_codon(codon)
    return codon


def _phenotype_correlation(df):
    pheno_colnames = _vhl_disease_colnames(df)
    pheno_pheno = pd.DataFrame(columns=pheno_colnames,
                               index=pheno_colnames)

    for col in pheno_colnames:
        phen_agg = df[df[col] >= 1].sum()
        pheno_pheno.loc[col] = phen_agg

    pheno_pheno = pheno_pheno.reindex(pheno_pheno.sum().sort_values(ascending=False).index)
    pheno_pheno = pheno_pheno[pheno_pheno.index.to_list()]

    pheno_pheno = pheno_pheno.rename(index=lambda x: x.split(".")[-1]).rename(columns=lambda x: x.split(".")[-1])

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
            # missense_domains,
            mutant_type_counts,
            mutant_type_ratios,
            codon_phenotype_subplots,
            codon_histogram,
            # codon_blosum62_histogram,
            # codon_blosum90_histogram,
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