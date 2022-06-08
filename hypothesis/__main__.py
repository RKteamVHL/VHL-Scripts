import logging
import os
import json
from .fetching.hypothesis_api import get_annotations_by_group
from .features.statistics import get_all_statistics


OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "files", "output")
ANNOTATION_OUTPUT = "hypothesis_annotations.json"

if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

PROBLEM_ANNOTATIONS = "problem_annotations.csv"
ANNOTATION_SUMMARY = "annotation_summary.csv"

GROUP_ID = "dKymJJpZ" #VHL annotations group

GROUP_EPOCH = "2019-08-27T00:00:00" # timestamp before all annotations

if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
        datefmt="%H:%M:%S")


    annotations = get_annotations_by_group(GROUP_ID, GROUP_EPOCH)
    with open(os.path.join(OUTPUT_DIR, ANNOTATION_OUTPUT), "w") as file:
        json.dump([a.as_dict() for a in annotations], file, indent=4)

    output_df, output_stats = get_all_statistics(annotations)
    output_df = output_df.reindex(sorted(output_df.columns), axis=1)
    output_df.to_csv(os.path.join(OUTPUT_DIR, "all_annotations.csv"))

    for stat in output_stats:
        stat.to_csv(os.path.join(OUTPUT_DIR, f"{stat.name}.csv"))

