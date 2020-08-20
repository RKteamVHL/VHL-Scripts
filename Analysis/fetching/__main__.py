from .KimStudents import KimStudents
import logging
if __name__ == '__main__':
	logging.basicConfig(
		level=logging.DEBUG,
		format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
		datefmt="%H:%M:%S")

	# test all functions
	fetcher = KimStudents()
	fetcher.process()
	fetcher.save_raw_file(fetcher.filename)
	fetcher.save_dsv_file("dsv_"+fetcher.filename)
