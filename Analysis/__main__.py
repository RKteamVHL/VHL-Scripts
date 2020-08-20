from .fetching.KimStudents import KimStudents
from .features.summary import Summary
from .features.age import Age
from .features.phenotype import Phenotype
import logging
if __name__ == '__main__':
	logging.basicConfig(
		level=logging.DEBUG,
		format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
		datefmt="%H:%M:%S")

	tables = [
		Summary(),
		Age(),
		Phenotype()
	]

	# test all functions
	fetcher = KimStudents()
	fetcher.process()
	for row in fetcher.rows:
		for table in tables:
			table.add_row(row)

	for table in tables:
		table.make_dataframe()
		table.normalize(norm_types=["zscore", "minmax"])
		table.to_csv(f'{table.name}.csv')



