import os

BASE_DIR = os.path.abspath(os.path.dirname(__file__))

INPUT_DIR = os.path.join(BASE_DIR, 'files', 'input')
OUTPUT_DIR = os.path.join(BASE_DIR, 'files', 'output')
SUMMARY_DIR = os.path.join(OUTPUT_DIR, "summary")
VALIDATION_DIR = os.path.join(OUTPUT_DIR, "validation")

# make directories if they don't exist
for dirname in [INPUT_DIR, OUTPUT_DIR, SUMMARY_DIR, VALIDATION_DIR]:
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

for a_type in ['patient', 'kindred', 'variant']:
    directory = os.path.join(OUTPUT_DIR, a_type)
    if not os.path.isdir(directory):
        os.makedirs(directory)
