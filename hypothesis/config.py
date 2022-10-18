import os

# constants for directories and files

INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files", "input")
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files", "output")
LIB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'files', 'lib')
ANNOTATION_OUTPUT = "hypothesis_annotations.json"

PROBLEM_ANNOTATIONS = "problem_annotations.csv"
ANNOTATION_SUMMARY = "annotation_summary.csv"

# constants relating to the hypothesis 'VHL Annotations' group
GROUP_ID = "dKymJJpZ" # VHL annotations group
GROUP_EPOCH = "2019-08-27T00:00:00" # timestamp before all annotations

# constants relating to VHL annotation processing
NULL_TERMS = ['n/a', 'N/A', 'NA']
NON_STANDARD_REFS = ['NM_000551.2', 'AF010238.1']
STANDARD_REFS = ['NM_000551.3', 'NM_000551.4']

# 'global' variables that change on runtime
USE_CACHE = False


# pre-startup code
if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

if not os.path.isdir(INPUT_DIR):
    os.makedirs(INPUT_DIR)

if not os.path.isdir(LIB_DIR):
    os.makedirs(LIB_DIR)
