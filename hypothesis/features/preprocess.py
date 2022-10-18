from typing import List
from ..annotations.Annotation import BODY_TAGS_NAME, TEXT_TAGS_NAME, AugmentedAnnotation, AnnotationType
from ..variant_functions import get_variant_by_caid, clinvarid_to_variant_dict, VHL_PHENOTYPES
import pandas as pd
import numpy as np
import re

def _mutant_string_to_list(x, generalize=True):
    """
    Converts a string of mutation type sequence ontology (SO) terms separated by a comma or semicolon into a list
    @param x: string of mutation types
    @param generalize: if True, the function takes the listed mutation type and converts it to a context-relevant
        mutation based on the SO tree. The conversion is usually an abstraction to a higher node on the SO tree
    @return:
    """
    so_list = re.split('[;,]', x)
    output_categories = []
    type_ = "mutant_type"
    if generalize:
        type_ = "generalized_" + type_
    for term in so_list:
        term = term.casefold().strip()
        try:
            if term not in NULL_TERMS:
                var_obo = None
                if generalize:
                    var_obo = vf.generalized_so_terms(term)
                else:
                    var_obo = vf.get_valid_obo(term)

                output_categories.append(f"{type_}.{var_obo}")
        except ValueError as e:
            logging.getLogger("mutant_type").warning(repr(e))
    return output_categories

def add_phenotype_columns(df):
    """
    Adds non-generalized phenotype data to a copy of the inputted dataframe and returns it
    @param df:
    @return:
    """
    hpo_series = df['Phenotype'].apply(_phenotype_string_to_list, generalize=False)

    pheno_counts = hpo_series.apply(collections.Counter)
    pheno_featurized = pd.DataFrame.from_records(pheno_counts)
    COMPUTED_COLUMNS["phenotype"].extend(pheno_featurized.columns.to_list())
    df = df.join(pheno_featurized)
    return df

def add_mutant_type_columns(df):
    """
    Adds generalized mutation-type data to a copy of the inputted dataframe and returns it
    @param df:
    @return df:
    """
    muttype_series = df[f'{BODY_TAGS_NAME}.MutationType']

    mutant_type_counts = variant_series.apply(collections.Counter)
    mutants_featurized = pd.DataFrame.from_records(mutant_type_counts)
    COMPUTED_COLUMNS["generalized_mutant_type"].extend(mutants_featurized.columns.to_list())
    df = df.join(mutants_featurized)
    return df


def preprocess_variant_annotations(annotations: List[AugmentedAnnotation]):
    annotation_df = AugmentedAnnotation.df_from_annotations(annotations)
    variant_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    output_df_list.append(df)

    return raw_df, output_df_list