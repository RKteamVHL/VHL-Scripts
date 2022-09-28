import os
import csv
import urllib.request
import urllib.parse
import urllib.error

import json

from . import config
from .fetching.clinvar_variants import fetch_clinvar_vhl_variants

CAID_URL = "https://reg.genome.network/allele/"

CLINVAR_FILE = os.path.join(config.INPUT_DIR, "clinvar_vhl.tsv")


def clinvarid_to_variant_dict():
    id_dict = {}
    if not config.USE_CACHE:
        fetch_clinvar_vhl_variants(CLINVAR_FILE)
    with open(CLINVAR_FILE, newline='\n') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            id_dict[int(row["VariationID"])] = row["Name"]

    return id_dict


def get_variant_by_caid(caid: str):
    url = f"{CAID_URL}{caid}"
    request = urllib.request.Request(url, method="GET")
    response_dict = None
    try:
        with urllib.request.urlopen(request) as response:
            response_dict = json.loads(response.read().decode('utf-8'))
    except urllib.error.HTTPError:
        pass

    return response_dict


VHL_PHENOTYPES = {
        'asymptomatic': 'asymptomatic',
        'adrenalpheochromocytoma': 'neuroendocrineneoplasm',
        'cerebellarhemangioblastoma': 'hemangioblastoma',
        'clearcellrenalcellcarcinoma': 'renalcellcarcinoma',
        'hemangioblastoma': 'hemangioblastoma',
        'hemangioma': 'hemangioma',
        'neoplasmoftheliver': 'neoplasmoftheliver',
        'neuroendocrineneoplasm': 'neuroendocrineneoplasm',
        'pancreaticcysts': 'pancreaticcysts',
        'pancreaticendocrinetumor': 'pancreaticendocrinetumor',
        'pheochromocytoma': 'neuroendocrineneoplasm',
        'renalcellcarcinoma': 'renalcellcarcinoma',
        'renalcyst': 'renalcyst',
        'renalneoplasm': 'renalneoplasm',
        'retinalcapillaryhemangioma': 'hemangioblastoma',
        'spinalhemangioblastoma': 'hemangioblastoma'
}


