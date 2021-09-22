import argparse

from .fetching.KimStudents import KimStudents

from .features.kimstudents_dataframe_preprocessing import *
from .features.kimstudents_dataframe_views import *
from .features.kimstudents_dataframe_clustering import *
from .features.kimstudents_dataframe_decisiontree import *
from .features.kimstudents_dataframe_stats import STAT_NAMES_FUNCTIONS, run_stats

from .validation.core import get_umd_variants

import numpy as np
import pandas as pd

OUTPUT_DIR = "output"

if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def groupby_patient(df):
    patient_df = df[df["Resolution"].str.casefold() == "patient"]
    return patient_df

def groupby_variant(df):
    variant_df_all = df[df["Resolution"].isin(["patient", "family", "tumour", "variant"])]
    variant_df = variant_df_all.set_index(["Mutation Event c.DNA."])

    variant_df_phens = variant_df[COMPUTED_COLUMNS["generalized_phenotype"]].groupby(variant_df.index).agg("sum")
    variant_df_phens[variant_df_phens >= 1] = 1
    variant_df_rest = variant_df.drop(columns=COMPUTED_COLUMNS["generalized_phenotype"]).groupby(variant_df.index).first()
    variant_df = variant_df_phens.join(variant_df_rest)
    return variant_df

def groupby_kindred(df):
    kindred_df_all = df[df["Resolution"].isin(["patient", "family"])]
    kindred_df = kindred_df_all.set_index(["PMID", "Kindred Case", "Mutation Event c.DNA."])
    kindred_df_phens = kindred_df[COMPUTED_COLUMNS["generalized_phenotype"]].groupby(kindred_df.index).agg("sum")
    kindred_df_phens[kindred_df_phens >= 1] = 1
    kindred_df_rest = kindred_df.drop(columns=COMPUTED_COLUMNS["generalized_phenotype"]).groupby(kindred_df.index).first()
    kindred_df = kindred_df_phens.join(kindred_df_rest)
    return kindred_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--cached', help="Load data from local cache", action="store_true")
    parser.add_argument('-figs', '--createfigs', help="Create all figures", action="store_true")
    parser.add_argument('-cl', '--cluster', help="Only redo clustering, rather than create all figures", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
        datefmt="%H:%M:%S")


    # test all functions
    fetcher = KimStudents()

    if args.cached:
        fetcher.load_from_dsv((os.path.join(OUTPUT_DIR, "data*.csv")))
    else:
        fetcher.process()
        fetcher.save_raw_file(os.path.join(OUTPUT_DIR, "data.csv"))


    out_table = pd.DataFrame(fetcher.rows)
    out_table = out_table.pipe(kimstudents_preprocessing)
    out_table["Resolution"] = out_table["Resolution"].str.casefold()

    #dropping all rows that dont have a single mutation type and phenotype
    # see pmid 28388566 in sheet: row is recorded as CHB, even though there is no mutation for this patient
    out_table = out_table.dropna(subset=[*COMPUTED_COLUMNS["generalized_phenotype"],
                                         *COMPUTED_COLUMNS["generalized_mutant_type"]], how='all')

    out_table.to_csv("unfiltered_out.csv")
    supplementary_table = out_table.pipe(create_supplementary_tables)
    supplementary_table.to_csv("supplementary_1.csv", index=False)
    if not os.path.isfile("umd.csv"):
        umd_variant_df = get_umd_variants()
        umd_variant_df.to_csv("umd.csv")

    umd_variant_df = pd.read_csv('umd.csv', dtype={'PMID':str})

    umd_variant_df['cdna_in_students'] = umd_variant_df['Mutation Event c.DNA.'].isin(out_table['Mutation Event c.DNA.'])
    umd_variant_df['pmid_in_students'] = umd_variant_df['PMID'].isin(out_table['PMID'])
    umd_variant_df.to_csv('umd_out.csv')


    out_df = {
        "patient": out_table.pipe(groupby_patient),
        "kindred": out_table.pipe(groupby_kindred),
        "variant": out_table.pipe(groupby_variant)
    }
    summary = pd.DataFrame(columns=list(STAT_NAMES_FUNCTIONS.keys()), index=list(out_df.keys()))

    for df_type, df_out in out_df.items():
        df_out.sum(skipna=True, numeric_only=True).to_csv(f"pre_dropna_{df_type}.csv")
        df_out[df_out[COMPUTED_COLUMNS["generalized_phenotype"]] == 0] = np.NaN
        df = df_out.dropna(subset=COMPUTED_COLUMNS["generalized_phenotype"], how='all')
        df = df.dropna(subset=COMPUTED_COLUMNS["generalized_mutant_type"], how='all')

        for statname, statfunc in STAT_NAMES_FUNCTIONS.items():
            summary.loc[df_type, statname] = statfunc(df)

        df.sum(skipna=True, numeric_only=True).to_csv(f"post_dropna_{df_type}.csv")

        if args.createfigs:
            if not args.cluster:
                create_descriptive_figures(df, df_type)

            clustered = dataframe_snf(df, df_type).fillna(0)


            create_cluster_summaries(clustered, df_type)
            create_cluster_phenotype_summaries(clustered, df_type)

            #disable the decision tree for now
            # decision = dataframe_decisiontree(df, df_type)


    plt.close('all')
    summary.to_csv("summary.csv")

    run_stats(os.path.join(STATS_DIR, "patient"))
    run_stats(os.path.join(STATS_DIR, "kindred"))

    refs = out_table[["PMID", "Reference"]]
    refs = refs.groupby(["Reference", "PMID"])
    refs.first().to_csv("all_refs.csv")

    vars = out_table.copy()[["HGVS_transcript"]]
    vars.loc[:, "HGVS_transcript"] = vars["HGVS_transcript"].str.split(";")
    vars = vars.explode(column="HGVS_transcript")
    vars.loc[:, "HGVS_transcript"] = vars["HGVS_transcript"].str.replace("NM_000551.3:", "")
    vars = vars.replace('', np.nan)
    vars = vars.dropna()
    vars.to_csv("all_variants.csv")


    # ax = summary.plot(kind='bar')
    # plt.xticks(rotation=0)
    # for p in ax.patches:
    #     ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    # plt.savefig(f'summary.pdf')


