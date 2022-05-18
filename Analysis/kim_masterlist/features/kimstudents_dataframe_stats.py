import os
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import binom_test
from scipy.stats import chi2_contingency
from scipy.stats import ks_2samp
from .kimstudents_dataframe_views import TESTS_DIR, DATA_DIR

PVALUE = 0.05

HIGHEST_N = 10



def get_asterisks_for_pval(p_val, alpha):
    """Receives the p-value and returns asterisks string."""
    if p_val > alpha:
        p_text = "ns"  # above threshold => not significant
    elif p_val/alpha < 1e-4:
        p_text = '****'
    elif p_val/alpha < 1e-3:
        p_text = '***'
    elif p_val/alpha < 1e-2:
        p_text = '**'
    else:
        p_text = '*'

    return p_text


def chisq_and_posthoc_corrected(df, outfile):
    """Receives a dataframe and performs chi2 test and then post hoc.
    Prints the p-values and corrected p-values (after FDR correction)"""
    # start by running chi2 test on the matrix
    chi2, p, dof, ex = chi2_contingency(df, correction=True)
    print(f"Chi2 result of the contingency table: {chi2}, p-value: {p}")

    # post-hoc
    all_combinations = list(combinations(df.index, 2))  # gathering all combinations for post-hoc chi2
    p_vals = []
    report_df = pd.DataFrame(index=all_combinations, columns=["p_value", "corrected", "reject", "significance"])

    for comb in all_combinations:
        new_df = df[(df.index == comb[0]) | (df.index == comb[1])]
        chi2, p, dof, ex = chi2_contingency(new_df, correction=True)
        p_vals.append(p)
        # print(f"For {comb}: {p}")  # uncorrected

    # checking significance
    # correction for multiple testing
    # reject_list, corrected_p_vals = multipletests(p_vals, method='fdr_bh')[:2]
    for p_val, comb in zip(p_vals, all_combinations):
        report_df.loc[[comb], "p_value"] = p_val
        report_df.loc[[comb], "corrected"] = PVALUE/len(p_vals)
        report_df.loc[[comb], "reject"] = p_val < PVALUE/len(p_vals)
        report_df.loc[[comb], "significance"] = get_asterisks_for_pval(p_val, PVALUE/len(p_vals))
        # print(
        #     f"{comb}: p_value: {p_val:5f}; corrected: {corr_p_val:5f} ({get_asterisks_for_pval(p_val)}) reject: {reject}")
    report_df.to_csv(outfile)


def binomial_and_posthoc_corrected(df, outfile):
    """Receives a dataframe and performs binomial test and then post hoc.
    Prints the p-values and corrected p-values (after FDR correction)"""

    n = df.iloc[:, 0].sum()
    n_categories = 213

    adjusted_p_thresh = PVALUE/n_categories

    report_df = pd.DataFrame(index=df.index, columns=["p_value", "corrected", "reject", "significance"])
    p_vals = []

    for i in range(len(df.index)):
        p = binom_test(df.iloc[i, 0], n=n, p=1/n_categories, alternative='greater')
        p_vals.append(p)
        # print(f"For {comb}: {p}")  # uncorrected

    # checking significance
    # correction for multiple testing
    # reject_list, corrected_p_vals = multipletests(p_vals, method='bonferroni')[:2]
    for p_val, cat in zip(p_vals, df.index):
        report_df.loc[cat, "p_value"] = p_val
        report_df.loc[cat, "corrected"] = adjusted_p_thresh
        report_df.loc[cat, "reject"] = p_val < adjusted_p_thresh
        report_df.loc[cat, "significance"] = get_asterisks_for_pval(p_val,adjusted_p_thresh)

        # print(
        #     f"{cat}: p_value: {p_val:5f}; corrected: {corr_p_val:5f} ({get_asterisks_for_pval(p_val)}) reject: {reject}")
    report_df.to_csv(outfile)

def ks_test(df, outfile):
    df = df.loc[:, (df != 0).any(axis=0)]
    index_list = list(df.index)
    column_list = list(df.columns)

    p_values = pd.DataFrame(index=column_list, columns=column_list)
    d_stat = pd.DataFrame(index=column_list, columns=column_list)
    test_df = pd.DataFrame(index=column_list, columns=column_list)
    num_tests = 0

    for j1 in range(0, len(column_list)-1):
        j1_data = df.index[df[column_list[j1]] == 1].tolist()
        for j2 in range(j1+1, len(column_list)):
            j2_data = df.index[df[column_list[j2]] == 1].tolist()

            stat, pvalue = ks_2samp(j1_data, j2_data)
            num_tests = num_tests + 1
            p_values.loc[column_list[j1], column_list[
                j2]] = pvalue
            d_stat.loc[column_list[j1], column_list[
                j2]] = stat


    alpha = PVALUE/num_tests

    for j1 in range(0, len(column_list)-1):
        j1_data = df.index[df[column_list[j1]] == 1].tolist()
        for j2 in range(j1+1, len(column_list)):
            j2_data = df.index[df[column_list[j2]] == 1].tolist()
            pvalue = p_values.loc[column_list[j1], column_list[j2]]
            stat = d_stat.loc[column_list[j1], column_list[j2]]
            ksresult = ""
            if pvalue < (alpha):
                ksresult = "**\n"
            test_df.loc[column_list[j1], column_list[j2]] = f"alpha:{alpha}\n{ksresult}p-value: {pvalue}\nD({len(j1_data)}, {len(j2_data)}): {stat}"
    test_df.to_csv(outfile)

def run_stats(root_dir):
    data_dir = os.path.join(root_dir, DATA_DIR)
    test_dir = os.path.join(root_dir, TESTS_DIR)
    if not os.path.isdir(test_dir):
        os.makedirs(test_dir)


    codon_pheno_df = pd.read_csv(os.path.join(data_dir, "phenotype_codon_heatmap.csv"), index_col=0)
    codon_df = pd.read_csv(os.path.join(data_dir, "codon_histogram.csv"), index_col=0)
    codon_blosum62_df = pd.read_csv(os.path.join(data_dir, "codon_blosum62_histogram.csv"), index_col=0)
    codon_blosum90_df = pd.read_csv(os.path.join(data_dir, "codon_blosum90_histogram.csv"), index_col=0)
    aachange_df = pd.read_csv(os.path.join(data_dir, "phenotype_aachange_heatmap.csv"), index_col=0)
    penetrance_df = pd.read_csv(os.path.join(data_dir, "penetrance.csv"), index_col=0)
    regions_ab_df = pd.read_csv(os.path.join(data_dir, "regions_alpha_beta.csv"), index_col=0)
    regions_elohif_df = pd.read_csv(os.path.join(data_dir, "regions_elongin_hifa.csv"), index_col=0)
    muttype_df = pd.read_csv(os.path.join(data_dir, "grouped_mutant_type_counts.csv"), index_col=0)

    # regions_df = regions_df.drop(columns=["ElonginB_ElonginC_binding", "HIF1_alpha_binding", "GXEEX8"])
    regions_ab_df = regions_ab_df[["region.⍺-Domain", "region.β-Domain"]]
    regions_elohif_df = regions_elohif_df.drop(columns=["region.GXEEX8"])

    binomial_and_posthoc_corrected(codon_df, os.path.join(test_dir, "codon_binom.csv"))
    binomial_and_posthoc_corrected(codon_blosum62_df, os.path.join(test_dir, "codon_binom_62.csv"))
    binomial_and_posthoc_corrected(codon_blosum90_df, os.path.join(test_dir, "codon_binom_90.csv"))
    chisq_and_posthoc_corrected(muttype_df, os.path.join(test_dir, "muttype_chisquare.csv"))
    chisq_and_posthoc_corrected(regions_ab_df, os.path.join(test_dir, "regions_ab_chisquare.csv"))
    chisq_and_posthoc_corrected(regions_elohif_df, os.path.join(test_dir, "regions_elohif_chisquare.csv"))
    ks_test(penetrance_df, os.path.join(test_dir, "penetrance_ks.csv") )

