import os
import io
import gzip
import urllib.request
import urllib.parse
import csv

csv.field_size_limit(1310720)

CLINVAR_SUMMARY_URI = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz'


def fetch_clinvar_vhl_variants(fileout):

    h_request = urllib.request.Request(CLINVAR_SUMMARY_URI, method="GET")

    with urllib.request.urlopen(h_request) as response:
        with gzip.open(response, mode='rt') as clinvar_tsv:
            reader = csv.DictReader(clinvar_tsv, delimiter='\t', quoting=csv.QUOTE_NONE)
            with open(fileout, 'w', newline='') as fout:
                writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter='\t')
                writer.writeheader()
                for row in reader:

                    if 'VHL' in row['GeneSymbol']:
                        writer.writerow(row)

