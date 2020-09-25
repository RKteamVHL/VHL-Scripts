import argparse

from .features.KimStudentsFeatureTables import *
from .fetching.KimStudents import KimStudents

OUTPUT_DIR = "output"

if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--cached', help="Load data from the inputted glob files", default="")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
        datefmt="%H:%M:%S")

    tables = [
        SummaryTable(),
        EvaluatedAgeBinnedTable(),
        LastKnownAgeBinnedTable(),
        IsolatedPhenotypeTable(),
        EvaluatedAgeTable(),
        PhenotypeTable(),
        VariantTypeTable(),
        PhenotypeCoocurrenceTable(),
        LossOfFunctionTable(),
        MutationEventTable()
    ]

    # test all functions
    fetcher = KimStudents()

    if args.cached != "":
        fetcher.load_from_dsv(args.cached)
    else:
        fetcher.process()
        fetcher.save_raw_file(os.path.join(OUTPUT_DIR, "data.csv"))

    resolution = ResolutionFeature()
    for row in fetcher.rows:
        tables[0].add_row(row)
        resolution.add_row(row)

    # fh = logging.FileHandler('phenotypes.log')
    # # fh.setLevel(logging.DEBUG)
    # logging.getLogger("phenotypes").addHandler(fh)
    #
    # fh = logging.FileHandler('variant_type.log')
    # # fh.setLevel(logging.DEBUG)
    # logging.getLogger("variant_type").addHandler(fh)

    for row in resolution.rows["patient"]:
        for table in tables[1:]:
            table.add_row(row)

    for table in tables:
        table.make_dataframe(include_count=False)
        # table.normalize(norm_types=['zscore', 'cdf'])
        table.to_csv(os.path.join(OUTPUT_DIR, f'{table.name}.csv'))
        table.dataframe.corr(method="pearson").to_csv(os.path.join(OUTPUT_DIR, f'{table.name}_corr.csv'))
        table.cluster(out_dir=OUTPUT_DIR)

# this is for error checking generalized phenotypes / variant types
# with open('phenotypes.log', 'r') as file_in:
# 	unique = set()
# 	for line in file_in:
# 		unique.add(line)
# 	with open('phenotypes.log', 'w') as file_out:
# 		for item in unique:
# 			file_out.write(item)
#
# with open('variant_type.log', 'r') as file_in:
# 	unique = set()
# 	for line in file_in:
# 		unique.add(line)
# 	with open('variant_type.log', 'w') as file_out:
# 		for item in unique:
# 			file_out.write(item)
