import logging
import os
import json
from .hypothesis_api import get_annotations_by_group
from .Annotation import get_invalid_dataframe


OUTPUT_DIR = "output"
OUTPUT_NAME = "hypothesis_annotations.json"

PROBLEM_ANNOTATIONS = "problem_annotations.csv"

GROUP_EPOCH = "2019-08-27T00:00:00"

if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
        datefmt="%H:%M:%S")

    annotations = get_annotations_by_group("dKymJJpZ", GROUP_EPOCH)
    with open(os.path.join(OUTPUT_DIR, OUTPUT_NAME), "w") as file:
        json.dump([a.as_dict() for a in annotations], file, indent=4)

    get_invalid_dataframe(annotations).to_csv(os.path.join(OUTPUT_DIR, PROBLEM_ANNOTATIONS))
