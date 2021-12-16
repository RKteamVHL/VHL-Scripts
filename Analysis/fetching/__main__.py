import logging

from .KimStudents import KimStudents

#fetching can be run as a module to test the fetching classes/functionality on a higher level
if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
        datefmt="%H:%M:%S")

    # test all functions
    fetcher = KimStudents()
    fetcher.process()
    fetcher.save_raw_file(fetcher.filename)
    fetcher.save_dsv_file("dsv_" + fetcher.filename)
