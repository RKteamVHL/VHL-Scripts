from .Fetcher import FETCHING_FACTORY
if __name__ == '__main__':

	# test all functions
	for key, value in FETCHING_FACTORY.items():
		fetcher = value()
		fetcher.process()
		# fetcher.save_raw_file(fetcher.filename)
		fetcher.save_dsv_file("dsv_"+fetcher.filename)
