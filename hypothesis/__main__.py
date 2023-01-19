import logging
import argparse
import os
import json
from .fetching.hypothesis_api import get_annotations_by_group, get_annotations_from_json
from .features.summary import get_all_summaries
from .features.preprocess import preprocess, unique_cases
from .figures.views import create_descriptive_figures
from . import config
from .annotations.Annotation import AugmentedAnnotation, AnnotationHeader

if __name__ == '__main__':
    logging.basicConfig(
        level=logging.DEBUG,
        format="[%(asctime)s] %(levelname)s [%(filename)s.%(funcName)s:%(lineno)d] %(message)s",
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

    # creating a df from the annotations
    annotation_df = AugmentedAnnotation.df_from_annotations(annotations)

    # preprocessing the df
    annotation_df = preprocess(annotation_df)

    # getting summary dfs
    summary_dfs = get_all_summaries(annotation_df)

    # outputting the annotation df and summary stats
    annotation_df.to_csv(os.path.join(config.DIRS['output'], config.RAW_ANNOTATION_DF))
    for stat in summary_dfs:
        stat.to_csv(os.path.join(config.DIRS['summary'], f"{stat.name}.csv"))

    # TODO: if this is done before the preprocessing some time will be saved
    case_df = unique_cases(annotation_df)
    case_df.to_csv(os.path.join(config.DIRS['output'], config.CASE_ANNOTATION_DF))

    create_descriptive_figures(config.DIRS['output'], {'case': case_df})

    # do clustering

    # do plotting

    # do stats




