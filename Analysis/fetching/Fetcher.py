import urllib.request
import urllib.parse
from civicpy import civic
from .. import variant_functions as vf
from .constants import *
import certifi
import argparse
import gzip
import json
import csv
import re
import os
import io

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

		print(self.name)


				
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
		"""Saves the post-processed dictionary list
		"""
		with open(filename, 'w', newline='', encoding='utf-8') as file:			
			writer = csv.DictWriter(file, fieldnames=self.dsv_header, delimiter=ROW_DELIMITER)
			writer.writeheader()
			writer.writerows(self.rows)

	def to_dict_list(self):
		"""Creates a dict list from output"""

		# find header list by reading the first line that
		# doesn't start with #
		header_list = filter(lambda row: not row.startswith('#'), self.data.read().splitlines())

		self.rows = csv.DictReader(header_list, delimiter=ROW_DELIMITER)

		self.dsv_header = list(self.rows.fieldnames)

			

	def fix_file(self):
		"""Overriden by subclass if data needs to be fixed
		
		By default, only the first requested file is saved.
		If no data is fetched into self.decompressed, this method
		does nothing
		"""
		if self.decompressed:
			self.data = io.TextIOWrapper(self.decompressed[0], encoding='utf-8')


	# TODO: currently in the subclasses, this method is being used for 
	# additional functionality that may not belong (namely, cdnaChange)
	# This should be changed in the future
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

		self.dsv_header = list(self.data.readline().split(ROW_DELIMITER))	
		self.dsv_header[0] = self.dsv_header[0].replace("#", "")
		self.dsv_header.append('cdnaChange')
		self.dsv_header.append('associatedPhenotypes')
		self.dsv_header.append('variantTypes')

		self.rows = csv.DictReader(self.data, fieldnames=self.dsv_header, delimiter=ROW_DELIMITER)

		self.filter_rows()

		for row in self.rows:
			row['cdnaChange'] = vf.get_valid_cdna(row['Name'], check_version=True)
			all_ids = re.split('[;,]', row['PhenotypeIDS'])
			hpo_re = re.compile('Human Phenotype Ontology:(?P<hpo>HP:[0-9]+)')

			hpo_full = filter(lambda aid: hpo_re.search(aid) is not None, all_ids)

			hpo_list = [hpo_re.search(hpoid).groupdict()['hpo'] for hpoid in hpo_full]

			row['associatedPhenotypes'] = [vf.get_valid_hpo(term) for term in hpo_list]
			row['variantTypes'] = [row['Type'].replace(' ', '_')]

	def filter_rows(self):
		# filter clinvar by VHL genesymbol
		self.rows = list(filter(lambda row: row['GeneSymbol']=='VHL', self.rows))

		# we might also want to filter by reviewer status, no conflicts 


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
			self.data.write(LINE_DELIMITER)
		
		self.data.seek(0)

	def to_dict_list(self):
		super().to_dict_list()
		self.dsv_header.append('cdnaChange')
		self.dsv_header.append('associatedPhenotypes')
		self.dsv_header.append('variantTypes')


		self.filter_rows()	
		for row in self.rows:
			row['cdnaChange'] = vf.get_valid_cdna(row['Mutation Event c.DNA.'])

			hpo_list  = re.split('[;,]', row['Phenotype'])
			row['associatedPhenotypes'] = [vf.get_valid_hpo(term) for term in hpo_list]

			row['variantTypes']  = re.split(';,', row['Mutation Type'])
			


	def filter_rows(self):
		#make sure row has an integer PMID
		self.rows = list(filter(lambda row: row['PMID'].isdigit(), self.rows))


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
		variants = gnomad_dict['data']['transcript']['variants']
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

			new_row['alleleCount'] = new_row.pop('ac', 0)
			new_row['alleleFrequency'] = new_row.pop('af', 0)
			new_row['cdnaChange'] = vf.get_valid_cdna(new_row['hgvsc'])
			# gnomad doesnt have phenotypes
			new_row['associatedPhenotypes'] = []
			new_row['variantTypes'] =[new_row['consequence']]

			self.dsv_header = list(new_row.keys())
			self.dsv_header.append('associatedPhenotypes')
			self.dsv_header.append('variantTypes')	

			self.rows.append(new_row)

		self.filter_rows()	



	def filter_rows(self):
		# gnomad has a quality filter, eliminating variants that dont meet it
		self.rows = list(filter(lambda row: len(row['filters'])==0, self.rows))		

			


class Civic(Fetcher):
	"""Fetcher for Civic database.
	"""
	def __init__(self):
		super().__init__()
		self.name = CIVIC_NAME
		self.filename = CIVIC_FILENAME



	def to_dict_list(self):
		self.rows = []
		for variant in civic.get_gene_by_id(CIVIC_VHL_GENE_ID).variants:
			new_row = {}

			# sets the next node id based on cdna change
			#this is used to automatically merge nodes across databases
			cdna_id = None

			# EXTRACTING CDS METHOD 1: using hgvs
			# TODO: cross-validate different hgvs' instead of overwrite 
			for hgvs in variant.hgvs_expressions:
				cdna_id = vf.get_valid_cdna(hgvs, check_version=True)

			
			# EXTRACTING CDS METHOD 2: using civic variant name
			cds_change = re.compile("\((c\..*)\)").search(variant.name)
			if cds_change is not None:
				cdna_id = cds_change.group(1)
				if cdna_id is None:
					cdna_id = cds_change.group(1)

			new_row['civicId'] = variant.id
			new_row['variantName'] = variant.name
			new_row['cdnaChange'] = cdna_id
			new_row['evidenceSubmitted'] = []
			new_row['evidenceAccepted'] = []
			new_row['phenotypesSubmitted'] = []
			new_row['phenotypesAccepted'] = []
			new_row['associatedPhenotypes'] = []
			new_row['variantTypes'] = []
			new_row['hgvsExpressions'] = []
			new_row['proteinChange'] = None


			#finding how many evidence items exist for the variant.
			#TODO: actually verify that evidence: supports, is germline, and is case study
			for evidence in variant.evidence_items:
				# get all phenotypes from evidences
				phenotypes = [phenotype.hpo_id for phenotype in evidence.phenotypes]

				# add evidence to appropriate list
				if evidence.status == "submitted":
					new_row['evidenceSubmitted'].append(evidence.id)
					new_row['phenotypesSubmitted'].extend(phenotypes)


				elif evidence.status == "accepted":
					new_row['evidenceAccepted'].append(evidence.id)	
					new_row['phenotypesAccepted'].extend(phenotypes)


			#finding the types of the variant and adding it to the node
			for variant_type in variant.variant_types:
				new_row['variantTypes'].append(variant_type.name)

			#finding the HGVS expressions
			for hgvs in variant.hgvs_expressions:
				new_row['hgvsExpressions'].append(hgvs)		

				# EXTRACTING AA CHANGE METHOD 1: via hgvs
				index = hgvs.find(vf.CURRENT_VHL_PROTEIN)
				if index >-1 :
					aa_change = hgvs[1+index+len(vf.CURRENT_VHL_PROTEIN):]
					new_row['proteinChange'] = aa_change	

			self.rows.append(new_row)
			self.dsv_header = list(new_row.keys())
		self.filter_rows()


	def filter_rows(self, accepted_only=False):
		# it's possible variants have no submitted or accepted evidence
		# in the case of a rejection; remove these
		self.rows = list(filter(lambda row: 
			len(row['evidenceAccepted'])>0 or len(row['evidenceSubmitted'])>0, 
			self.rows
		))

		# filter rows based on optional parameter
		if accepted_only:
			self.rows = filter(lambda row: len(row['evidenceAccepted'])>0, self.rows)

		for row in self.rows:
			valid_accepted = [vf.get_valid_hpo(term) for term in row['phenotypesAccepted']]
			row['associatedPhenotypes'].extend(valid_accepted)
			if not accepted_only:
				valid_submitted = [vf.get_valid_hpo(term) for term in row['phenotypesSubmitted']]
				row['associatedPhenotypes'].extend(valid_submitted)



	


FETCHING_FACTORY = {
	'ClinVar': ClinVar,
	'KimStudents2019': KimStudents2019,
	'Gnomad': Gnomad,
	'Civic': Civic
}

