import os

# constants
INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files", "input")

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files", "output")
ANNOTATION_OUTPUT = "hypothesis_annotations.json"

PROBLEM_ANNOTATIONS = "problem_annotations.csv"
ANNOTATION_SUMMARY = "annotation_summary.csv"

GROUP_ID = "dKymJJpZ" #VHL annotations group

GROUP_EPOCH = "2019-08-27T00:00:00" # timestamp before all annotations

# 'global' variables that change on runtime
USE_CACHE = False

if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

if not os.path.isdir(INPUT_DIR):
    os.makedirs(INPUT_DIR)
