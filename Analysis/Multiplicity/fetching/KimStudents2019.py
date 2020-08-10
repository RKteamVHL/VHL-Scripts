from .Fetcher import Fetcher
from .constants import *
from ... import variant_functions as vf
import logging
import csv
import re
import io

STUDENTS_NAME = 'KimStudents2019'
STUDENTS_FILENAME = 'CiVIC Extracted Data (Summer 2019-Present).txt'

# Civic extracted
STUDENTS_SHEETIDS = ['1826130937', '0', '763995397', '511845134', '301841530']
STUDENTS_HREF = [f'https://docs.google.com/spreadsheets/d/1YYXHLM9dOtGyC1l8edjQnfcAf0IM1EuTrz8ZKbQkQGU/export?format=tsv&gid={name}' for name in STUDENTS_SHEETIDS]
STUDENTS_HREF.append("https://docs.google.com/spreadsheets/d/1-4nc9VIYP2YQ2gNR7NuIYEnzcLjTpt5m3szz6HElL_Y/export?format=tsv&gid=0")


STUDENTS_HEADER = [
	"PMID",
	"cDNA_Position",
	"Multiple Mutants in Case",
	"Mutation Event c.DNA.",
	"Predicted Consequence Protein Change",
	"variant_name",
	"Mutation Type",
	"Kindred Case",
	"Familial/Non-familial",
	"Phenotype",
	"Reference",
	"Age",
	"Notes",
	"Evidence Statements"
]


class KimStudents2019(Fetcher):
	"""Fetcher for Raymond's custom student-procured database.
	"""

	def __init__(self):
		super().__init__()
		self.name = STUDENTS_NAME
		self.filename = STUDENTS_FILENAME
		self.href = STUDENTS_HREF
		self.needs_extraction = False
		self.logger = logging.getLogger(self.name)
		self.dsv_header = STUDENTS_HEADER

	def fix_file(self):
		# combining all sheets together; skip first header line in each worksheet
		for decompressed in self.decompressed:
			text_wrapper = io.TextIOWrapper(decompressed, encoding='utf-8')
			self.data.write(text_wrapper.read())
			self.data.write(LINE_DELIMITER)

		self.data.seek(0)

	def to_dict_list(self):

		all_rows = filter(lambda r: not r.startswith('#'), self.data.read().splitlines())

		self.rows = csv.DictReader(all_rows, fieldnames=self.dsv_header, delimiter=ROW_DELIMITER)

		self.filter_rows()
		for row in self.rows:

			# TODO: eventually this should be automatically verified
			row['cdnaChange'] = row['Mutation Event c.DNA.']

			hpo_list = re.split('[;,]', row['Phenotype'])
			so_list = re.split('[;,]', row['Mutation Type'])

			row['associatedPhenotypes'] = []
			for term in hpo_list:
				try:
					var_hpo = vf.get_valid_obo(term.strip())
					row['associatedPhenotypes'].append(var_hpo)
				except ValueError as e:
					self.logger.warning(repr(e))

			row['variantTypes'] = []
			for term in so_list:
				try:
					var_so = vf.get_valid_obo(term)
					row['variantTypes'].append(var_so)
				except ValueError as e:
					self.logger.warning(repr(e))

			row['PMID'] = [row['PMID']]

			row["proteinChange"] = [row["Predicted Consequence Protein Change"]]

	def filter_rows(self):
		out_rows = []
		for row in self.rows:
			out_row = {}
			for k, v in row.items():
				if k is not None:
					out_row[k] = v
			if row['PMID'].isdigit():
				out_rows.append(out_row)
		self.rows = out_rows
