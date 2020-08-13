from .fetching.KimStudents import KimStudents
from .features.summary import SummaryCalculator
import logging
if __name__ == '__main__':
	logging.basicConfig(
		level=logging.DEBUG,
		format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
		datefmt="%H:%M:%S")

	summary = SummaryCalculator()
	# test all functions
	fetcher = KimStudents()
	fetcher.process()
	for row in fetcher.rows:
		summary.add_row(row)
	summary.to_csv("summary.csv")

