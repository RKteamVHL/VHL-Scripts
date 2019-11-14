import urllib.request
import urllib.parse
#from . import variant_functions as vf
import certifi
import argparse
import gzip
import json
import csv
import os
import io

TAB_DELIMITER= '\t'
PIPE_DELIMITER= "|"


FILES_DIRECTORY = os.path.join('ApiScripts','files')

CLINVAR_NAME = 'ClinVar'
CLINVAR_FILENAME = "variant_summary.tsv"
CLINVAR_HREF = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_HEADER = "AlleleID{0}Type{0}Name{0}GeneID{0}GeneSymbol{0}HGNC_ID{0}ClinicalSignificance{0}ClinSigSimple{0}LastEvaluated{0}RS# (dbSNP){0}nsv/esv (dbVar){0}RCVaccession{0}PhenotypeIDS{0}PhenotypeList{0}Origin{0}OriginSimple{0}Assembly{0}ChromosomeAccession{0}Chromosome{0}Start{0}Stop{0}ReferenceAllele{0}AlternateAllele{0}Cytogenetic{0}ReviewStatus{0}NumberSubmitters{0}Guidelines{0}TestedInGTR{0}OtherIDs{0}SubmitterCategories{0}VariationID\n".format(
	TAB_DELIMITER)


STUDENTS_NAME = 'KimStudents2019'
STUDENTS_FILENAME = 'CiVIC Extracted Data (Summer 2019-Present).tsv'
STUDENTS_SHEETIDS = ['1826130937', '0', '763995397', '511845134', '301841530']
STUDENTS_HREF = [f'https://docs.google.com/spreadsheets/d/1YYXHLM9dOtGyC1l8edjQnfcAf0IM1EuTrz8ZKbQkQGU/export?format=tsv&gid={name}' for name in STUDENTS_SHEETIDS ]
STUDENTS_HEADER = "PMID{0}cDNA_Position{0}Multiple Mutants in Case{0}Mutation Event c.DNA.{0}Predicted Consequence Protein Change{0}variant_name{0}Mutation Type{0}Kindred Case{0}Familial/Non-familial{0}Phenotype{0}Reference{0}Age{0}Notes{0}Evidence Statements".format(
	TAB_DELIMITER)


GNOMAD_POST_DATA = {
	"query": '''{
		gene(gene_id: "ENSG00000134086", reference_genome: GRCh37) {
			variants(dataset: gnomad_r2_1) {
				consequence
				flags
				hgvs
				hgvsc
				hgvsp
				lof
				lof_filter
				lof_flags
				pos
				rsid
				variantId
				exome {
					ac
					ac_hemi
					ac_hom
					an
					af
					filters
					populations {
						id
						ac
						an
						ac_hemi
						ac_hom
					}
				}
				genome {
					ac
					ac_hemi
					ac_hom
					an
					af
					filters
					populations {
						id
						ac
						an
						ac_hemi
						ac_hom
					}
				}
			}
		}
	}'''
}

GNOMAD_NAME = 'gnomAD_v2_1'
GNOMAD_FILENAME = 'gnomAD_v2_1.tsv'
GNOMAD_HREF = 'http://gnomad.broadinstitute.org/api'

SO_NAME = 'SequenceOntology'
SO_FILENAME = 'so.obo'
SO_HREF = 'https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/so.obo'

HPO_NAME = 'HumanPhenotypeOntology'
SO_FILENAME = 'hp.obo'
SO_HREF = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'


class Fetcher(object):
	"""Base class for fetching data from external sources.

	Includes basic funtionality for reading files from href links,
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
		self.compressed = []

		# stores the decompressed file byte data- this is the same as
		# the compressed data if it doesn't need extraction
		self.decompressed = []


		# holds the final data, after fixes (eg, headers) have been
		# made to the file
		self.data = io.StringIO() 

		# flags if the file needs to be extracted
		self.needs_extraction = False

	def fetch(self):
		"""Reads data from instance's href attribute."""

		# href is converted to list of 1 element if it's a string
		if isinstance(self.href, str):
			self.href = [self.href]


		if isinstance(self.href, list):
			for href in self.href:
				with urllib.request.urlopen(href, data=self.request_data, cafile=certifi.where()) as response:
					compressed_bytes = io.BytesIO()
					compressed_bytes.write(response.read())
					compressed_bytes.seek(0)
					self.compressed.append(compressed_bytes)			

		print(self.href)


				
	def extract(self):
		"""Extracts data in compressed attribute to decompressed using gzip.
		This method will do nothing to the data if the file doesn't need extraction
		"""
		for compressed_bytes in self.compressed:
			decompressed_bytes = io.BytesIO()
			if self.needs_extraction:			
				with gzip.open(compressed_bytes) as decompressed:
					decompressed_bytes.write(decompressed.read())
			else:
				decompressed_bytes.write(compressed_bytes.read())

			decompressed_bytes.seek(0)
			self.decompressed.append(decompressed_bytes)


	def save_raw_file(self, filename):
		"""Save the data file to the given filename

		Args:
			filename: filename (str) to save data
		"""
		with open(filename, 'w+', encoding='utf-8') as file:
			file.write(self.data.read())
			self.data.seek(0)

	def save_dsv_file(self, filename):
		with open(filename, 'w', newline='', encoding='utf-8') as file:			
			writer = csv.DictWriter(file, fieldnames=self.dsv_header, delimiter=TAB_DELIMITER)
			writer.writeheader()
			writer.writerows(self.rows)

	def to_dict_list(self):
		"""Creates a dict list from output"""

		# find header list by reading the first line that
		# doesn't start with #
		header_list = filter(lambda row: not row.startswith('#'), self.data.read().splitlines())

		self.rows = csv.DictReader(header_list, delimiter=TAB_DELIMITER)

		self.dsv_header = list(self.rows.fieldnames)

		self.data.seek(0)

			

	def fix_file(self):
		"""Overriden by subclass if headers need to be fixed
		By default, only the first requested file is saved
		"""
		self.data = io.TextIOWrapper(self.decompressed[0], encoding='utf-8')
		pass



class ClinVar(Fetcher):
	"""Fetcher for ClinVar database.
	"""
	def __init__(self):
		super().__init__()
		self.name = CLINVAR_NAME
		self.filename = CLINVAR_FILENAME
		self.href = CLINVAR_HREF
		self.needs_extraction = True

	def to_dict_list(self):		
		# Clinvar's header line starts with a '#'

		self.dsv_header = list(self.data.readline().split(TAB_DELIMITER))	
		self.dsv_header[0] = self.dsv_header[0].replace("#", "")	

		self.rows = csv.DictReader(self.data, fieldnames=self.dsv_header, delimiter=TAB_DELIMITER)

		self.rows = filter(lambda row: row['GeneSymbol']=='VHL', self.rows)

		self.data.seek(0)

class KimStudents2019(Fetcher):
	"""Fetcher for Raymond's custom student-procured database.
	"""
	def __init__(self):
		super().__init__()
		self.name = STUDENTS_NAME
		self.filename = STUDENTS_FILENAME
		self.href = STUDENTS_HREF
		self.needs_extraction = False

	def fix_file(self):
		# combining all sheets together
		for decompressed in self.decompressed:
			text_wrapper = io.TextIOWrapper(decompressed, encoding='utf-8')
			self.data.write(text_wrapper.read())
			self.data.write('\n')
		
		self.data.seek(0)


class Gnomad(Fetcher):
	"""Fetcher for Gnomad database.

	"""
	def __init__(self):
		super().__init__()
		self.name = GNOMAD_NAME
		self.filename = GNOMAD_FILENAME
		self.href = GNOMAD_HREF
		self.needs_extraction = False
		self.request_data =  urllib.parse.urlencode(GNOMAD_POST_DATA).encode('ascii')

	def to_dict_list(self):
		self.rows = []		
		gnomad_dict = json.loads(self.data.read())
		variants = gnomad_dict['data']['gene']['variants']
		for variant in variants:
			new_row = variant

			#combining exome and genome data
			has_genome = new_row['genome'] is not None
			has_exome = new_row['exome'] is not None
			for key in ["ac", "ac_hemi", "ac_hom", "an", "af"]:
				new_row[key] = 0
				new_row[key] += new_row['genome'].get(key, 0) if has_genome else 0
				new_row[key] += new_row['exome'].get(key, 0) if has_exome else 0
			
			new_row['filters'] = []
			new_row['filters'].extend(new_row['genome']['filters'] if has_genome else [])
			new_row['filters'].extend(new_row['exome']['filters'] if has_exome else [])

			new_row['af'] = new_row['ac']/new_row['an']
			new_row.pop('exome', None)
			new_row.pop('genome', None)
			self.dsv_header = list(new_row.keys())

			# gnomad has a quality filter, eliminating variants that dont meet it
			if len(new_row['filters']) == 0:
				self.rows.append(new_row) 
		self.data.seek(0)
		


FETCHING_DICT = {
	'ClinVar': ClinVar,
	'KimStudents2019': KimStudents2019,
	'Gnomad': Gnomad
}

if __name__ == '__main__':

	for key, value in FETCHING_DICT.items():
		fetcher = value()
		fetcher.fetch()
		fetcher.extract()
		fetcher.fix_file()
		fetcher.to_dict_list()
		fetcher.save_raw_file(fetcher.filename)
		fetcher.save_dsv_file("dsv_"+fetcher.filename)
