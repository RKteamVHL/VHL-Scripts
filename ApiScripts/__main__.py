from civicpy import civic
import argparse
import time
import json
import csv
import os 

import pprint

# 58 is gene VHL
GENE_ID =  58

if __name__ == '__main__':
	parser = argparse.ArgumentParser()


	# parser.add_argument('-exp', '--experiment_directory', help = '''Experiment Directory''')
	# parser.add_argument('-s', '--save', help = '''Save all databases localy.''', action="store_true")

	gene = civic.get_gene_by_id(GENE_ID)
	#making a dictionary of variant types
	variant_type_dict = {}

	for variant in gene.variants:
		for var_type in variant.types:
			if not var_type.name in variant_type_dict:
				variant_type_dict[var_type.name] = []
			variant_type_dict[var_type.name].append(variant.name)
			# using the actual civic ids instead of mutation type names
			# if not var_type.id in variant_type_dict:
			# 	variant_type_dict[var_type.id] = []
			# variant_type_dict[var_type.id].append(variant)

	pprint.pprint(variant_type_dict)

# QUESTIONS:
# - where do you find the keys each attribute type has?
# - is there a list of ids for mutation types (Civic Attribute)