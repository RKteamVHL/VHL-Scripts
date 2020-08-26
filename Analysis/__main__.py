from .fetching.KimStudents import KimStudents
from .features.KimStudentsFeatures import ResolutionFeature
from .features.KimStudentsFeatureTables import *
import logging
import argparse

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
		EvaluatedAgeTable(),
		LastKnownAgeTable(),
		IsolatedPhenotypeTable(),
		PhenotypeTable(),
		VariantTypeTable()
	]

	# test all functions
	fetcher = KimStudents()

	if args.cached != "":
		fetcher.load_from_dsv(args.cached)
	else:
		fetcher.process()
		fetcher.save_raw_file("data.csv")

	resolution = ResolutionFeature()
	for row in fetcher.rows:
		tables[0].add_row(row)
		resolution.add_row(row)

	for row in resolution.rows["patient"]:
		for table in tables[1:]:
			table.add_row(row)

	for table in tables:
		table.make_dataframe()
		table.normalize(norm_types=["zscore", "minmax", "cdf"])
		table.to_csv(f'{table.name}.csv')



