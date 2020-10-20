import argparse

from .fetching.KimStudents import KimStudents

from .features.kimstudents_dataframe_preprocessing import *
from .features.kimstudents_dataframe_views import *
from .features.kimstudents_dataframe_clustering import *

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

    out_table = out_table[out_table["Resolution"].str.casefold() == "patient"]
    # out_table = out_table.dropna(subset=COMPUTED_COLUMNS['generalized_mutant_type'], how='all')
    summary = pd.DataFrame(columns=["Number of Observations"], index=["patient", "kindred"])
    summary.loc['patient', "Number of Observations"] = len(out_table)
    if not args.cluster:
        create_figures(out_table, "by_patient")

    # Clustering
    clustered = dataframe_snf(out_table, "by_patient").fillna(0)

    # grouping by kindred
    out_table = out_table.set_index(["PMID", "Kindred Case", "Mutation Event c.DNA."])
    out_table_phens = out_table[COMPUTED_COLUMNS["generalized_phenotype"]].groupby(out_table.index).agg("sum")
    out_table_phens[out_table_phens >= 1] = 1
    out_table_rest = out_table.drop(columns=COMPUTED_COLUMNS["generalized_phenotype"]).groupby(out_table.index).first()
    out_table = out_table_phens.join(out_table_rest)
    #
    summary.loc['kindred', "Number of Observations"] = len(out_table)
    if not args.cluster:
        create_figures(out_table, "by_kindred")
    clustered = dataframe_snf(out_table, "by_kindred").fillna(0)

    plt.close('all')
    ax = summary.plot(kind='bar')
    plt.xticks(rotation=0)
    for p in ax.patches:
        ax.annotate(str(p.get_height()), (p.get_x() * 1.005, p.get_height() * 1.005))
    plt.savefig(f'summary.pdf')
