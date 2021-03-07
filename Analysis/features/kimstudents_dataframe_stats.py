from .kimstudents_dataframe_preprocessing import COMPUTED_COLUMNS
from .kimstudents_dataframe_views import TESTS_DIR, DATA_DIR
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import ks_2samp
import pandas as pd
import os

PVALUE = 0.05

HIGHEST_N = 10


def total_observations(df):
    return len(df.index)


def total_with_sex_data(df):
    df_sex = df.dropna(subset=["sex.f", "sex.m"], how="all")
    return len(df_sex.index)

def total_female(df):
    df_sex = df.dropna(subset=COMPUTED_COLUMNS['sex'], how="all")
    df_female = df_sex.dropna(subset=["sex.f"])
    return len(df_female.index)

def total_male(df):
    df_sex = df.dropna(subset=COMPUTED_COLUMNS['sex'], how="all")
    df_male = df_sex.dropna(subset=["sex.m"])
    return len(df_male.index)

def total_with_age_data(df):
    df_age = df.dropna(subset=COMPUTED_COLUMNS['age'], how="all")
    return len(df_age.index)

def total_with_lk_age(df):
    df_age = df.dropna(subset=COMPUTED_COLUMNS['age'], how="all")
    df_lk = df_age.dropna(subset=["last_known_age"], how="all")
    return len(df_lk.index)

def total_with_onset_age(df):
    df_age = df.dropna(subset=COMPUTED_COLUMNS['age'], how="all")
    df_onset = df_age.dropna(subset=["evaluated_age"], how="all")
    return len(df_onset.index)

def lk_age_mean(df):
    df_age = df.dropna(subset=COMPUTED_COLUMNS['age'], how="all")
    df_lk = df_age.dropna(subset=["last_known_age"], how="all")
    return df_lk["last_known_age"].mean()

def onset_age_mean(df):
    df_age = df.dropna(subset=COMPUTED_COLUMNS['age'], how="all")
    df_onset = df_age.dropna(subset=["evaluated_age"], how="all")
    return df_onset["evaluated_age"].mean()

def total_denovo(df):
    df_denovo = df[df["denovo.yes"] == 1]
    return len(df_denovo.index)

def total_with_phenotype(df):
    df = df.replace(0, np.NaN)
    df_pheno = df.dropna(subset=COMPUTED_COLUMNS["generalized_phenotype"], how='all')
    return len(df_pheno.index)

def total_with_muttype(df):
    df_muttype = df.dropna(subset=COMPUTED_COLUMNS["generalized_mutant_type"], how='all')
    return len(df_muttype.index)

def fisher_combinations(df, outfile):
    
    index_list = list(df.index)
    column_list = list(df.columns)
    data = df.to_numpy()[:, 0:HIGHEST_N]

    p_values = pd.DataFrame(index=column_list, columns=column_list)
    d_stat = pd.DataFrame(index=column_list, columns=column_list)
    test_df = pd.DataFrame(index=column_list, columns=column_list)

    num_tests = 0
    for i1 in range(0, data.shape[0] - 1):
        for i2 in range(i1 + 1, data.shape[0]):
            for j1 in range(0, data.shape[1] - 1):
                for j2 in range(j1 + 1, data.shape[1]):
                    num_tests = num_tests + 1

    alpha = PVALUE / num_tests

    with open(outfile, 'w', encoding="utf-8") as file:
        for i1 in range(0, data.shape[0]-1):
            for i2 in range(i1+1, data.shape[0]):
                for j1 in range(0, data.shape[1]-1):
                    for j2 in range(j1+1, data.shape[1]):
                        contingency = np.array([
                            [data[i1, j1], data[i1, j2]],
                            [data[i2, j1], data[i2, j2]],
                        ])
                        oddsratio, pvalue = fisher_exact(contingency)
                        if pvalue < alpha:
                            file.write(f"alpha: {alpha}\np-value: {pvalue}\nOdds ratio: {oddsratio}\n\t{index_list[i1]}\t{index_list[i2]}\n\t{column_list[j1]}\t{column_list[j2]}\n")

def chi_test(df, outfile):
    df = df.loc[:, (df != 0).any(axis=0)]
    index_list = list(df.index)
    column_list = list(df.columns)

    p_values = pd.DataFrame(index=column_list, columns=column_list)
    chi_stat = pd.DataFrame(index=column_list, columns=column_list)
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
            chi_stat.loc[column_list[j1], column_list[
                j2]] = stat


    alpha = PVALUE/num_tests

    for j1 in range(0, len(column_list)-1):
        j1_data = df.index[df[column_list[j1]] == 1].tolist()
        for j2 in range(j1+1, len(column_list)):
            j2_data = df.index[df[column_list[j2]] == 1].tolist()
            pvalue = p_values.loc[column_list[j1], column_list[j2]]
            stat = chi_stat.loc[column_list[j1], column_list[j2]]
            ksresult = ""
            if pvalue < (alpha):
                ksresult = "**\n"
            test_df.loc[column_list[j1], column_list[j2]] = f"alpha:{alpha}\n{ksresult}p-value: {pvalue}\nD({len(j1_data)}, {len(j2_data)}): {stat}"
            test_df.to_csv(outfile)

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

    codon_name = "phenotype_codon_heatmap.csv"
    aachange_name = "phenotype_aachange_heatmap.csv"
    penetrance_name = "penetrance.csv"
    regions_name = "regions.csv"
    muttype_name = "grouped_mutant_type_counts.csv"

    codon_df = pd.read_csv(os.path.join(data_dir, codon_name), index_col=0)
    aachange_df = pd.read_csv(os.path.join(data_dir, aachange_name), index_col=0)
    penetrance_df = pd.read_csv(os.path.join(data_dir, penetrance_name), index_col=0)
    regions_df = pd.read_csv(os.path.join(data_dir, regions_name), index_col=0)
    muttype_df = pd.read_csv(os.path.join(data_dir, muttype_name), index_col=0)

    regions_df = regions_df.drop(columns=["ElonginB_ElonginC_binding", "HIF1_alpha_binding", "GXEEX8"])

    fisher_combinations(codon_df, os.path.join(test_dir, "codon_fisher.txt"))
    fisher_combinations(aachange_df, os.path.join(test_dir, "aachange_fisher.txt"))
    fisher_combinations(regions_df, os.path.join(test_dir, "region_fisher.txt"))
    fisher_combinations(muttype_df, os.path.join(test_dir, "muttype_fisher.txt"))
    ks_test(penetrance_df, os.path.join(test_dir, "penetrance_ks.csv") )


STAT_NAMES_FUNCTIONS = {
    "Number of Observations": total_observations,
    "Total N with sex data": total_with_sex_data,
    "Total male": total_male,
    "Total female": total_female,
    "Total N with age data": total_with_age_data,
    "Total with last known age": total_with_lk_age,
    "Total with onset age": total_with_onset_age,
    "Mean last known age": lk_age_mean,
    "Mean onset age": onset_age_mean,
    "N confirmed de novo": total_denovo,
    "N with phenotype": total_with_phenotype,
    "N with muttype": total_with_muttype,
}
