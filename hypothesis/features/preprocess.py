from typing import List
from ..annotations.Annotation import BODY_TAGS_NAME, TEXT_TAGS_NAME, AugmentedAnnotation, AnnotationType
from ..fetching.caid_variants import get_variant_by_caid
from ..fetching.clinvar_variants import clinvarid_to_variant_dict
from ..variant_functions import DISEASE_ENTITY_TO_HPO
import pandas as pd

import re

from ..annotations.Annotation import AnnotationType


def unique_cases(df):
	out_df = df[df['type'] == AnnotationType.CASE]
	return out_df


def preprocess_variant_annotations(annotations: List[AugmentedAnnotation]):
    annotation_df = AugmentedAnnotation.df_from_annotations(annotations)
    variant_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    output_df_list.append(df)

    return raw_df, output_df_list