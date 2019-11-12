import urllib.request
import certifi
import argparse
import gzip
import os
import io

TAB_DELIMITER= '\t'
PIPE_DELIMITER= "|"


FILES_DIRECTORY = os.path.join('ApiScripts','files')

STUDENTS_FILE = os.path.join(FILES_DIRECTORY, 'CiVIC Extracted Data (Summer 2019-Present).xlsx')

GNOMAD_FILE = os.path.join(FILES_DIRECTORY, 'gnomAD_v2.1.1_ENSG00000134086_2019_10_04_17_21_52.csv')



CLINVAR_NAME = 'ClinVar'
CLINVAR_FILENAME = "variant_summary.tsv"
CLINVAR_HREF = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_HEADER = "AlleleID{0}Type{0}Name{0}GeneID{0}GeneSymbol{0}HGNC_ID{0}ClinicalSignificance{0}ClinSigSimple{0}LastEvaluated{0}RS# (dbSNP){0}nsv/esv (dbVar){0}RCVaccession{0}PhenotypeIDS{0}PhenotypeList{0}Origin{0}OriginSimple{0}Assembly{0}ChromosomeAccession{0}Chromosome{0}Start{0}Stop{0}ReferenceAllele{0}AlternateAllele{0}Cytogenetic{0}ReviewStatus{0}NumberSubmitters{0}Guidelines{0}TestedInGTR{0}OtherIDs{0}SubmitterCategories{0}VariationID\n".format(
	TAB_DELIMITER)


STUDENTS_NAME = 'KimStudents2019'
STUDENTS_FILENAME = 'CiVIC Extracted Data (Summer 2019-Present).tsv'
STUDENTS_SHEETIDS = ['0', '1826130937', '763995397', '511845134', '301841530']
STUDENTS_HREF = [f'https://docs.google.com/spreadsheets/d/1YYXHLM9dOtGyC1l8edjQnfcAf0IM1EuTrz8ZKbQkQGU/export?format=tsv&gid={name}' for name in STUDENTS_SHEETIDS ]
STUDENTS_HEADER = "PMID{0}cDNA_Position{0}Multiple Mutants in Case{0}Mutation Event c.DNA.{0}Predicted Consequence Protein Change{0}variant_name{0}Mutation Type{0}Kindred Case{0}Familial/Non-familial{0}Phenotype{0}Reference{0}Age{0}Notes{0}Evidence Statements".format(
	TAB_DELIMITER)

GNOMAD_NAME = 'gnomAD_v2'
GNOMAD_FILENAME = 'gnomAD_v2.tsv'
GNOMAD_HREF = ''
GNOMAD_HEADERS = "Chromosome{0}Position{0}rsID{0}Reference{0}Alternate{0}Source{0}Filters - exomes{0}Filters - genomes{0}Consequence{0}Protein Consequence{0}Transcript Consequence{0}Annotation{0}Flags{0}Allele Count{0}Allele Number{0}Allele Frequency{0}Homozygote Count{0}Hemizygote Count{0}Allele Count African{0}Allele Number African{0}Homozygote Count African{0}Hemizygote Count African{0}Allele Count Latino{0}Allele Number Latino{0}Homozygote Count Latino{0}Hemizygote Count Latino{0}Allele Count Ashkenazi Jewish{0}Allele Number Ashkenazi Jewish{0}Homozygote Count Ashkenazi Jewish{0}Hemizygote Count Ashkenazi Jewish{0}Allele Count East Asian{0}Allele Number East Asian{0}Homozygote Count East Asian{0}Hemizygote Count East Asian{0}Allele Count European (Finnish){0}Allele Number European (Finnish){0}Homozygote Count European (Finnish){0}Hemizygote Count European (Finnish){0}Allele Count European (non-Finnish){0}Allele Number European (non-Finnish){0}Homozygote Count European (non-Finnish){0}Hemizygote Count European (non-Finnish){0}Allele Count Other{0}Allele Number Other{0}Homozygote Count Other{0}Hemizygote Count Other{0}Allele Count South Asian{0}Allele Number South Asian{0}Homozygote Count South Asian{0}Hemizygote Count South Asian".format(
	TAB_DELIMITER)

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
			Used to fix headers of tsv/csv files, if needed
	"""

	def __init__(self):
		"""Inits the Fetcher class and Attributes."""
		self.name = None
		self.filename = None
		self.href = None

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
		print(self.href)
		if isinstance(self.href, list):
			for href in self.href:
				with urllib.request.urlopen(href, cafile = certifi.where()) as response:
					compressed_bytes = io.BytesIO()
					compressed_bytes.write(response.read())
					compressed_bytes.seek(0)
					self.compressed.append(compressed_bytes)			

		elif isinstance(self.href, str):
			with urllib.request.urlopen(self.href, cafile = certifi.where()) as response:
				compressed_bytes = io.BytesIO()
				compressed_bytes.write(response.read())
				compressed_bytes.seek(0)
				self.compressed.append(compressed_bytes)
		

				
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


	def save_to_file(self, filename):
		"""Save the data file to the given filename

		Args:
			filename: filename (str) to save data
		"""
		with open(filename, 'w+', encoding='utf-8') as file:
			file.write(self.data.read())
			

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

	def fix_file(self):
		# clinvar's first character is a '#', and tossing it corrects its header 
		self.data = io.TextIOWrapper(self.decompressed[0], encoding='utf-8')
		self.data.read(1)

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


FETCHING_DICT = {
	'ClinVar': ClinVar,
	'KimStudents2019': KimStudents2019,
	#'Gnomad': Gnomad
}

if __name__ == '__main__':

	for key, value in FETCHING_DICT.items():
		fetcher = value()
		fetcher.fetch()
		fetcher.extract()
		fetcher.fix_file()
		fetcher.save_to_file(fetcher.filename)
