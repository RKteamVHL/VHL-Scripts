import argparse

from .fetching.KimStudents import KimStudents

from .features.kimstudents_dataframe_preprocessing import *
from .features.kimstudents_dataframe_views import *
from .features.kimstudents_dataframe_clustering import *
from .features.kimstudents_dataframe_decisiontree import *

from .validation.core import get_umd_variants

OUTPUT_DIR = "output"

if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--cached', help="Load data from local cache", action="store_true")
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

    #dropping all rows that dont have a single mutation type or phenotype
    # see pmid 28388566: row is recorded as CHB, even though there is no mutation for this patient
    out_table = out_table.dropna(subset=COMPUTED_COLUMNS["generalized_phenotype"], how='all')
    out_table = out_table.dropna(subset=COMPUTED_COLUMNS["generalized_mutant_type"], how='all')
    summary = pd.DataFrame(columns=["Number of Observations"], index=["patient", "kindred", "variant"])

    if not os.path.isfile("umd.csv"):
        umd_variant_df = get_umd_variants()
        umd_variant_df.to_csv("umd.csv")

    umd_variant_df = pd.read_csv('umd.csv', dtype={'PMID':str})

    umd_variant_df['cdna_in_students'] = umd_variant_df['Mutation Event c.DNA.'].isin(out_table['Mutation Event c.DNA.'])
    umd_variant_df['pmid_in_students'] = umd_variant_df['PMID'].isin(out_table['PMID'])
    umd_variant_df.to_csv('umd_out.csv')


    # Patient resolution

    patient_df = out_table[out_table["Resolution"].str.casefold() == "patient"]
    # out_table = out_table.dropna(subset=COMPUTED_COLUMNS['generalized_mutant_type'], how='all')
    summary.loc['patient', "Number of Observations"] = len(patient_df)

    if not args.cluster:
        create_figures(patient_df, "by_patient")

    clustered = dataframe_snf(patient_df, "by_patient").fillna(0)
    decision = dataframe_decisiontree(patient_df, "by_patient")


    # Variant resolution
    variant_df_all = out_table[out_table["Resolution"].isin(["patient", "family", "tumour", "variant"])]
    variant_df = variant_df_all.set_index(["Mutation Event c.DNA."])

    variant_df_phens = variant_df[COMPUTED_COLUMNS["generalized_phenotype"]].groupby(variant_df.index).agg("sum")
    variant_df_phens[variant_df_phens >= 1] = 1
    variant_df_rest = variant_df.drop(columns=COMPUTED_COLUMNS["generalized_phenotype"]).groupby(variant_df.index).first()
    variant_df = variant_df_phens.join(variant_df_rest)

    summary.loc['variant', "Number of Observations"] = len(variant_df)
    if not args.cluster:
        create_figures(variant_df, "by_variant")

    clustered = dataframe_snf(variant_df, "by_variant").fillna(0)
    decision = dataframe_decisiontree(variant_df, "by_variant")

    # Kindred resolution
    kindred_df_all = out_table[out_table["Resolution"].isin(["patient", "family"])]
    kindred_df = kindred_df_all.set_index(["PMID", "Kindred Case", "Mutation Event c.DNA."])
    kindred_df_phens = kindred_df[COMPUTED_COLUMNS["generalized_phenotype"]].groupby(kindred_df.index).agg("sum")
    kindred_df_phens[kindred_df_phens >= 1] = 1
    kindred_df_rest = kindred_df.drop(columns=COMPUTED_COLUMNS["generalized_phenotype"]).groupby(kindred_df.index).first()
    kindred_df = kindred_df_phens.join(kindred_df_rest)
    #
    summary.loc['kindred', "Number of Observations"] = len(kindred_df)
    if not args.cluster:
        create_figures(kindred_df, "by_kindred")
    decision = dataframe_decisiontree(kindred_df, "by_kindred")
    clustered = dataframe_snf(kindred_df, "by_kindred").fillna(0)

    plt.close('all')
    ax = summary.plot(kind='bar')
    plt.xticks(rotation=0)
    for p in ax.patches:
        ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    plt.savefig(f'summary.pdf')
