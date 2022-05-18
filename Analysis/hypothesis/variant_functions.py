import os
import csv

CLINVAR_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files", "input", "clinvar_vhl.tsv")


def clinvarid_to_variant_dict():
    id_dict = {}
    with open(CLINVAR_FILE, newline='\n') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            id_dict[int(row["VariationID"])] = row["Name"]

    return id_dict