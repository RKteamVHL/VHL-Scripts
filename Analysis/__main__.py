import argparse

from .features.kimstudents_dataframe_preprocessing import *
from .features.kimstudents_dataframe_stats import run_stats
from .features.kimstudents_dataframe_views import *
from .fetching.KimStudents import KimStudents
from .validation.core import *
from .features.kimstudents_dataframe_summaries import *

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

def filter_phenotype_mutanttype(df):
    df[df[COMPUTED_COLUMNS["generalized_phenotype"]] == 0] = np.NaN
    df = df.dropna(subset=COMPUTED_COLUMNS["generalized_phenotype"], how='all')
    df = df.dropna(subset=COMPUTED_COLUMNS["generalized_mutant_type"], how='all')
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--cached', help="Load data from local cache", action="store_true")
    parser.add_argument('-figs', '--createfigs', help="Create all figures", action="store_true")
    parser.add_argument('-cl', '--cluster', help="Only redo clustering, rather than create all figures", action="store_true")
    parser.add_argument('-va', '--validation', help="Toggle true to run validation scripts", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
        datefmt="%H:%M:%S")


    # fetching all sheets and combining into a masterlist
    fetcher = KimStudents()

    if args.cached:
        fetcher.load_from_dsv((os.path.join(OUTPUT_DIR, "data*.csv")))
    else:
        fetcher.process()
        fetcher.save_raw_file(os.path.join(OUTPUT_DIR, "data.csv"))


    # generate all clean columns needed for further anaysis
    out_table = pd.DataFrame(fetcher.rows)
    out_table = out_table.pipe(kimstudents_preprocessing)
    # generating all pre-drop summaries
    raw_out_df = {
        "patient": out_table.pipe(groupby_patient),
        "kindred": out_table.pipe(groupby_kindred),
        "variant": out_table.pipe(groupby_variant)
    }
    create_predrop_summary_table(raw_out_df)
    create_predropsupplementary_table(out_table)
    #dropping all rows that dont have a single mutation type and phenotype
    out_table = out_table.dropna(subset=[*COMPUTED_COLUMNS["generalized_phenotype"],
                                         *COMPUTED_COLUMNS["generalized_mutant_type"]], how='all')

    create_filtered_table(out_table)
    create_postdropsupplementary_table(out_table)
    if args.validation:
        create_umd_validation_table(out_table)
        create_litvar_variant_table(out_table)

    filtered_out_df = {
        "patient": out_table.pipe(groupby_patient).pipe(filter_phenotype_mutanttype),
        "kindred": out_table.pipe(groupby_kindred).pipe(filter_phenotype_mutanttype),
        "variant": out_table.pipe(groupby_variant).pipe(filter_phenotype_mutanttype)
    }
    create_postdrop_summary_table(filtered_out_df)

    if args.createfigs:
        if not args.cluster:
            create_descriptive_figures(filtered_out_df)
        create_cluster_figures(filtered_out_df)


    create_type_summary_tables(filtered_out_df)


    plt.close('all')


    run_stats(os.path.join(STATS_DIR, "patient"))
    run_stats(os.path.join(STATS_DIR, "kindred"))

    create_refs_table(out_table)

    create_vars_table(out_table)



