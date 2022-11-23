import urllib.request
import urllib
import json
import csv
import os
from .. import config
CAID_URL = "http://reg.genome.network/alleles?gene=VHL"
LIMIT = 5000
# CAID_URL = "https://reg.genome.network/allele/"

CAID_HEADERS = ['caid', 'communitystandardtitle']

CAID_FILE = os.path.join(config.DIRS['lib'], "caid_vhl.tsv")


def caid_to_variant():
    id_dict = {}
    if not config.USE_CACHE:
        fetch_caid_vhl_alleles(CAID_FILE)
    with open(CAID_FILE, newline='\n') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            id_dict[row['caid']] = row["communitystandardtitle"]

    while True:
        yield id_dict


def fetch_caid_vhl_alleles(fileout):
    total_read = 0
    num_read = LIMIT
    with open(fileout, 'w', newline='') as fout:

        caid_writer = csv.DictWriter(fout, fieldnames=CAID_HEADERS, delimiter='\t')
        caid_writer.writeheader()

        while num_read == LIMIT:
            num_read = 0
            h_request = urllib.request.Request(f'{CAID_URL}&limit={LIMIT}&skip={total_read}', method="GET")
            with urllib.request.urlopen(h_request) as response:
                response_dict = json.loads(response.read().decode('utf-8'))

                for allele_dict in response_dict:
                    if all([x in allele_dict for  x in ['@id', 'communityStandardTitle']]):
                        out_dict = {
                            'caid': allele_dict["@id"].split('/')[-1],
                            'communitystandardtitle': allele_dict['communityStandardTitle'][0]
                        }

                        caid_writer.writerow(out_dict)
                    num_read += 1

                total_read += num_read


caid_to_variant_generator = caid_to_variant()