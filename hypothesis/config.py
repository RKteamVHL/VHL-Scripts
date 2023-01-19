import os

## constants for directories and files
# base directories
DIRS = {
    'input': os.path.join(os.path.dirname(os.path.abspath(__file__)), "files", "input"),
    'output': os.path.join(os.path.dirname(os.path.abspath(__file__)), "files", "output"),
    'lib': os.path.join(os.path.dirname(os.path.abspath(__file__)), "files", "lib"),
}
# derrived dirs
DIRS['summary'] = os.path.join(DIRS['output'], 'summary')


# annotation file names
ANNOTATION_OUTPUT = "hypothesis_annotations.json"
PROBLEM_ANNOTATIONS = "problem_annotations.csv"
ANNOTATION_SUMMARY = "annotation_summary.csv"
RAW_ANNOTATION_DF = "all_annotations.csv"
CASE_ANNOTATION_DF = "case_annotations.csv"

# constants relating to the hypothesis 'VHL Annotations' group
GROUP_ID = "dKymJJpZ" # VHL annotations group
GROUP_EPOCH = "2019-08-27T00:00:00" # timestamp before all annotations

# constants relating to VHL annotation processing
NULL_TERMS = ['n/a', 'N/A', 'NA']
NON_STANDARD_REFS = ['NM_000551.2', 'AF010238.1']
STANDARD_REFS = ['NM_000551.3', 'NM_000551.4']

# constant for determining what gets returned by HPO/SO
OBO_RETURN_TYPE = 'name_spaceless'

# 'global' variables that change on runtime
USE_CACHE = False
CASEFOLD_TAG_NAMES = False



# pre-startup code
for directory_type, directory in DIRS.items():
    if not os.path.isdir(directory):
        os.makedirs(directory)
