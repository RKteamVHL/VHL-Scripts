import csv
import glob
import gzip
import io
import urllib.parse
import urllib.request
from itertools import chain

import certifi

from Analysis.kim_masterlist.constants import *


class Fetcher(object):
    """Base class for fetching data from external sources.

    Includes basic functionality for reading files from href links,
    extracting data (if href links to a zip file), and uploading the data to
    Subclasses of Fetcher should override class methods and properties depending
    on the specific database files.

    Attributes:
        name: name of the database. used to name the file.
        filename: name of the file (without url) being fetched.
        href: the url reference to the file to download.
        compressed: bytestream of the raw compressed data.
        decompressed: bytestream of data after decompression.
        data: the final stream of data that gets written to a file.
            Used to fix files, if needed
        request_data: dictionary of the data to be used in the body of the post request
    """

    def __init__(self):
        """Inits the Fetcher class and Attributes."""
        self.name = None
        self.filename = None
        self.href = None

        self.rows = None
        self.dsv_header = None

        # stores any data needed for post requests
        self.request_data = None

        # stores the compressed file byte data
        self.compressed = {}

        # stores the decompressed file byte data- this is the same as
        # the compressed data if it doesn't need extraction
        self.decompressed = {}

        # holds the final data, after fixes (eg, headers) have been
        # made to the file
        self.data = {}

        # flags if the file needs to be extracted
        self.needs_extraction = False

    def fetch(self):
        """Reads data from instance's href attribute."""

        # href is converted to list of 1 element if it's a string
        print(self.name)
        if isinstance(self.href, str):
            self.href = [self.href]

        if isinstance(self.href, list):
            for href in self.href:
                print(href)
                with urllib.request.urlopen(href, data=self.request_data, cafile=certifi.where()) as response:
                    compressed_bytes = io.BytesIO()
                    compressed_bytes.write(response.read())
                    compressed_bytes.seek(0)
                    self.compressed[href] = compressed_bytes

    def extract(self):
        """Extracts data in compressed attribute to decompressed using gzip.
        This method will do nothing to the data if the file doesn't need extraction
        """
        for href, compressed_bytes in self.compressed.items():
            decompressed_bytes = io.BytesIO()
            if self.needs_extraction:
                with gzip.open(compressed_bytes) as decompressed:
                    decompressed_bytes.write(decompressed.read())
            else:
                decompressed_bytes.write(compressed_bytes.read())

            decompressed_bytes.seek(0)
            self.decompressed[href] = decompressed_bytes

    def save_raw_file(self, filename):
        """Save the data file to the given filename

        Args:
            filename: filename (str) to save data
        """
        k = list(self.data.keys())
        for i in range(len(k)):
            with open(f'{filename[0:-4]}{i}{filename[-4:]}', 'w+', encoding='utf-8') as file:
                self.data[k[i]].seek(0)
                file.write(self.data[k[i]].read())
                self.data[k[i]].seek(0)

    def save_dsv_file(self, filename):
        """Saves the post-processed dictionary list
        """
        with open(filename, 'w', newline='', encoding='utf-8') as file:
            writer = csv.DictWriter(file, fieldnames=self.dsv_header, delimiter=ROW_DELIMITER)
            writer.writeheader()
            writer.writerows(self.rows)

    def to_dict_list(self):
        """Creates a dict list from output"""

        def empty_gen():
            yield from ()

        all_rows = empty_gen
        for href in self.data.keys():
            all_rows = chain(all_rows, filter(lambda row: not row.startswith('#'), self.data[href].read().splitlines()))

        self.rows = csv.DictReader(all_rows, delimiter=ROW_DELIMITER())
        self.dsv_header = list(self.rows.fieldnames)

    def fix_file(self):
        """Overriden by subclass if data needs to be fixed

        If no data is fetched into self.decompressed, this method
        does nothing
        """
        if self.decompressed:
            for href, decompressed_data in self.decompressed.items():
                self.data[href] = io.TextIOWrapper(decompressed_data, encoding='utf-8')

    def filter_rows(self):
        """Filters dictionary list
        Implemented by child classes, if needed
        """
        pass

    def process(self):
        """Performs all steps for file extraction + processing
        """
        self.fetch()
        self.extract()
        self.fix_file()
        self.to_dict_list()

    def load_from_dsv(self, file_glob):
        """Loads multiple csv files with identical headers and adds them to the fetcher row list

        """
        self.rows = []
        self.dsv_header = []
        for filename in glob.iglob(file_glob):
            with open(filename, 'r', encoding='utf-8') as file:
                self.rows.extend(list(csv.DictReader(file, delimiter=ROW_DELIMITER)))
            # self.dsv_header.extend(list(self.rows.fieldnames))
