from .Fetcher import Fetcher
from ..constants import *
from itertools import chain
import logging
import csv

STUDENTS_NAME = 'KimStudents'
STUDENTS_FILENAME = 'Final Masterlist of VHL Papers.csv'

# Garrett (2016): 1607707061
# Liam: 1618569887
# No PMID: 975723769
# Garrett (students): 64075395
# Kelly: 960822541
# Safa/Sam: 154810211
# Sean (VA): 410096135
# Andreea pre-2016: 1835654764
# Andreea ML: 859989440
STUDENTS_SHEETIDS = ['960822541', '1618569887', '154810211', '410096135', '975723769', '64075395', '1607707061', '859989440', '1835654764']
STUDENTS_HREF = [f'https://docs.google.com/spreadsheets/d/1evUglhZpGDUPm3GUXu0uIZSkZsHXnbbfQzl5_Nc8DAQ/export?format=csv&gid={name}' for name in STUDENTS_SHEETIDS]

STUDENTS_HEADER_NAMES = [
	"Checked on CIViC",
	"PMID",
	"cDNA_Position",
	"Multiple Mutants in Case",
	"Mutation Event c.DNA.",
	"Transcript Reference",
	"Predicted Consequence Protein Change",
	"variant_name",
	"Mutation Type",
	"Kindred Case",
	"Confirmed De Novo",
	"Phenotype (for reference)",
	"Phenotype",
	"Reference",
	"Age",
	"Sex",
	"Notes",
	"Evidence Statements",
	"Resolution",
	"associatedPhenotypes",
	"variantTypes",
	"Kindred Case (pedigree)"
]


class KimStudents(Fetcher):
	"""Fetcher for Raymond's custom student-procured database.
	"""

	def __init__(self):
		super().__init__()
		self.name = STUDENTS_NAME
		self.filename = STUDENTS_FILENAME
		self.href = STUDENTS_HREF
		self.needs_extraction = False
		self.logger = logging.getLogger(self.name)
		self.dsv_header = STUDENTS_HEADER_NAMES


	def to_dict_list(self):

		def empty_gen():
			yield from ()

		self.rows = empty_gen()
		for href in self.data.keys():
			new_rows = filter(lambda row: not row.startswith('#'), self.data[href].read().splitlines())

			self.rows = chain(self.rows, csv.DictReader(new_rows, delimiter=ROW_DELIMITER))

		self.filter_rows(lambda row:  row['PMID'].isdigit())


	def filter_rows(self, filt_fun):
		out_rows = []
		for row in self.rows:
			if filt_fun(row):
				out_rows.append(row)
		self.rows = out_rows
