import os
import numpy as np
import pandas as pd
from .kimstudents_dataframe_preprocessing import COMPUTED_COLUMNS
STATS_DIR = "statistics"
SUMMARY_DIR = "summary"
summary_path = os.path.join(STATS_DIR, SUMMARY_DIR)
if not os.path.isdir(summary_path):
    os.makedirs(summary_path)

SUPPLEMENTARY_HEADERS = {
    "Kindred Case": "Kindred Case",
    "Mutation Type": "Mutation Type",
    "Transcript Reference": "Transcript Reference",
    "HGVS_transcript": "HGVS Transcript",
    "HGVS_Predicted_Protein": "HGVS Predicted Protein",
    "Multiple Mutants in Case":  "Multiple Mutants in Case",
    "Confirmed De Novo": "Confirmed De Novo",
    "Phenotype": "Phenotype",
    "Age": "Age",
    "Sex": "Sex",
    "Resolution": "Resolution",
    "Reference": "Reference",
    "PMID": "PMID"
}





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


#this table will have all entries that have at least a phenotype and/or mutant type
def create_filtered_table(df):
    df.to_csv(os.path.join(summary_path, "filtered_out.csv"))

def create_refs_table(df):
    refs = df[["PMID", "Reference"]]
    refs = refs.groupby(["Reference", "PMID"])
    refs.first().to_csv(os.path.join(summary_path, "all_refs.csv"))

def create_vars_table(df):
    vars = df.copy()[["HGVS_transcript"]]
    vars.loc[:, "HGVS_transcript"] = vars["HGVS_transcript"].str.split(";")
    vars = vars.explode(column="HGVS_transcript")
    vars.loc[:, "HGVS_transcript"] = vars["HGVS_transcript"].str.replace("NM_000551.3:", "")
    vars = vars.replace('', np.nan)
    vars = vars.dropna()
    vars.to_csv(os.path.join(summary_path,"all_variants.csv"))


def create_type_summary_tables(dfs):
    cols = list([*COMPUTED_COLUMNS["generalized_phenotype"], *COMPUTED_COLUMNS["generalized_mutant_type"]])
    summary = pd.DataFrame(columns=cols, index=list(dfs.keys()))
    for df_type, df_out in dfs.items():
        summary.loc[df_type, cols] = df_out[cols].sum()

    summary.to_csv(os.path.join(summary_path, "summary_by_type.csv"))
def create_summary_table(dfs, prefix):
    summary = pd.DataFrame(columns=list(STAT_NAMES_FUNCTIONS.keys()), index=list(dfs.keys()))

    for df_type, df_out in dfs.items():
        for statname, statfunc in STAT_NAMES_FUNCTIONS.items():
            summary.loc[df_type, statname] = statfunc(df_out)

    summary.to_csv(os.path.join(summary_path, prefix+"summary.csv"))

def create_predrop_summary_table(out_df):
    create_summary_table(out_df, "predrop")

def create_postdrop_summary_table(out_df):
    create_summary_table(out_df, "postdrop")

def create_supplementary_table(df, prefix):
    df_trimmed = df[list(SUPPLEMENTARY_HEADERS.keys())]
    df_trimmed = df_trimmed.rename(columns=SUPPLEMENTARY_HEADERS)
    df_sorted = df_trimmed.sort_values(by=["Reference"])
    df_sorted = df_sorted.replace(r'^\s*$', np.nan, regex=True)
    df_sorted = df_sorted.dropna(how='all')
    df_sorted = df_sorted.replace(np.nan, "N/A", regex=True)
    supplementary_table = df_sorted
    supplementary_table.to_csv(os.path.join(summary_path,prefix+"supplementary_1.csv"), index=False)

def create_predropsupplementary_table(df):
    create_supplementary_table(df, "predrop")

def create_postdropsupplementary_table(df):
    create_supplementary_table(df, "postdrop")


