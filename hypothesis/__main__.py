import logging
import argparse
import os
import json
from .fetching.hypothesis_api import get_annotations_by_group, get_annotations_from_json
from .features.summary import get_all_summaries
from . import config
from .annotations.Annotation import AugmentedAnnotation

if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s",
        datefmt="%H:%M:%S")

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--cached', help="Load data from local cache", action="store_true")

    args = parser.parse_args()

    config.USE_CACHE = args.cached

    if config.USE_CACHE:
        # load the stored annotation json file
        annotations = get_annotations_from_json(os.path.join(config.DIRS['output'], config.ANNOTATION_OUTPUT))

    else:
        # dump the raw annotations to a json file
        annotations = get_annotations_by_group(config.GROUP_ID, config.GROUP_EPOCH)
        with open(os.path.join(config.DIRS['output'], config.ANNOTATION_OUTPUT), "w") as file:
            json.dump([a.as_dict() for a in annotations], file, indent=4)

    # merging article information into all other annotations
    AugmentedAnnotation.merge_across_document_title_and_source(annotations)

    # converting to csv and computing summary
    get_all_summaries(annotations)



