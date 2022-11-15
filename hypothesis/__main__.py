import logging
import argparse
import os
import json
from .fetching.hypothesis_api import get_annotations_by_group, get_annotations_from_json
from .features.summary import get_all_statistics
from . import config
from .annotations.Annotation import AugmentedAnnotation
from .features.filter import keep_only_cases

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
        annotations = get_annotations_from_json(os.path.join(config.OUTPUT_DIR, config.ANNOTATION_OUTPUT))

    else:
        # dump the raw annotations to a json file
        annotations = get_annotations_by_group(config.GROUP_ID, config.GROUP_EPOCH)
        with open(os.path.join(config.OUTPUT_DIR, config.ANNOTATION_OUTPUT), "w") as file:
            json.dump([a.as_dict() for a in annotations], file, indent=4)

    # merging article information into all other annotations
    AugmentedAnnotation.merge_across_source(annotations)

    # converting to csv and computing summary
    output_df, output_stats = get_all_statistics(annotations)
    output_df = output_df.reindex(sorted(output_df.columns), axis=1)
    output_df.to_csv(os.path.join(config.OUTPUT_DIR, "all_annotations.csv"))

    for stat in output_stats:
        stat.to_csv(os.path.join(config.OUTPUT_DIR, f"{stat.name}.csv"))

    # keeping only CASE annotations
    cases_df = keep_only_cases(output_df)
    cases

