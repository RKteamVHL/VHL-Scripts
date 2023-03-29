from typing import List
from ..annotations.Annotation import AnnotationHeader
from ..fetching.caid_variants import caid_to_variant_generator
from ..fetching.clinvar_variants import clinvarid_to_variant_generator

from .. import config
from .. import variant_functions
import pandas as pd
import numpy as np
import collections

import re

from ..annotations.Annotation import AnnotationType


# TODO: this file will become very large, and should probably separated into a python module at some point
# TODO: describe some of the pandas magic


def unique_cases(df: pd.DataFrame):
    """
    Get all unique cases from input annotation df
    @param df:
    @return:
    """
    out_df = df[df[str(AnnotationHeader.TYPE)] == "CASE"]
    # TODO: create index based on case, pmid, and kindred. grouping by 'first' may not be the right solution
    out_df = out_df.set_index([str(AnnotationHeader.TYPE), str(AnnotationHeader.PMID_CLEAN),
                               str(AnnotationHeader.KINDRED_ID_CLEAN)])
    out_df = out_df.groupby(out_df.index).first()
    return out_df


def fix_na(df: pd.DataFrame):
    """Given a df, replace all NULL terms with np.nan
    """

    def _fix_na(_df):
        return_df = _df.replace(config.NULL_TERMS, [np.nan] * len(config.NULL_TERMS))
        return return_df

    out_df = df.pipe(_fix_na)
    return out_df


def fix_obos(series: pd.Series):
    """	Given a series of strings with hpo or so terms, get their valid terms if they exist
    """
    out_series = series.apply(variant_functions.get_valid_obo, return_as=config.OBO_RETURN_TYPE)
    return out_series


def clean_pmid(df: pd.DataFrame):
    pmid_colname = str(AnnotationHeader.PMID)
    # this df magic explodes the lists into separate rows, splits them based on spaces, then tries to convert it to int
    pmid_series = pd.to_numeric(df[pmid_colname].explode().str.split(' ', n=1, expand=True)[0], errors='coerce')
    pmid_series = pmid_series.groupby(pmid_series.index).first()

    df[str(AnnotationHeader.PMID_CLEAN)] = pmid_series
    return df

def _caid_to_variant(caid):
    """Given a caid, return its full variant name
    """
    variant_dict = next(caid_to_variant_generator)
    variant_name = None
    if caid in variant_dict:
        # caid is guaranteed to be in the dict, but this checks for None
        variant_name = variant_dict[caid]
    return variant_name


def _fix_caid(caid_list):
    """Given a list in the unclean format: ['CAID, yes'], get a valid caid for the CAID in the unclean format """
    if isinstance(caid_list, list):
        caid = caid_list[0].split(', ')[0]
        variant_dict = next(caid_to_variant_generator)
        caid_clean = None
        if caid in variant_dict:
            caid_clean = caid
        return caid_clean


def caid_to_variant(df: pd.DataFrame):
    """
    Fix the caid column in the given df, and append two columns: the cleaned caid and the caid variant
    @param df:
    @return:
    """

    _caid = str(AnnotationHeader.CA_ID)
    _caid_clean = str(AnnotationHeader.CA_ID_CLEAN)
    _caid_variant = str(AnnotationHeader.CA_ID_VARIANT)

    caid_id_series = df[_caid].apply(_fix_caid)
    caid_id_series.name = _caid_clean

    caid_variant_series = caid_id_series.apply(_caid_to_variant)
    caid_variant_series.name = _caid_variant

    df = pd.concat([df, caid_id_series, caid_variant_series], axis=1)

    return df


def _fix_clinvar(clinvar_list):
    """Given a list of potential clinvar IDs, validate and return it
    """
    clinvar_id = None
    if isinstance(clinvar_list, list):
        for item in clinvar_list:
            if str.isdigit(item):
                clinvar_id = int(item)
    return clinvar_id


def _clinvar_to_variant(clinvar_id):
    """Given a clinvar id, return its full variant name
    """
    variant_dict = next(clinvarid_to_variant_generator)
    variant_name = None
    if clinvar_id in variant_dict:
        variant_name = variant_dict[clinvar_id]
    return variant_name


def clinvar_to_variant(df: pd.DataFrame):
    """
    Fix the clinvar column in the given df, and append two columns: the cleaned clinvar and the clinvar variant
    @param df:
    @return:
    """

    _clinvar = str(AnnotationHeader.CLINVAR_ID)
    _clinvar_clean = str(AnnotationHeader.CLINVAR_ID_CLEAN)
    _clinvar_variant = str(AnnotationHeader.CLINVAR_ID_VARIANT)

    clinvar_id_series = df[_clinvar].apply(_fix_clinvar)
    clinvar_id_series.name = _clinvar_clean

    clinvar_variant_series = clinvar_id_series.apply(_clinvar_to_variant)
    clinvar_variant_series.name = _clinvar_variant

    df = pd.concat([df, clinvar_id_series, clinvar_variant_series], axis=1)

    return df


def _separate_phenotype_age_of_presentation(aop_series: pd.Series):
    """Given a series where each cell is a list of dicts with keys as phenotypes and values as age of presentations,
    create a df where columns are phenotypes and cells are the first age of presentation for each phenotype
    """
    # separate the age dicts onto separate lines, convert to float
    aop_df = aop_series.explode().apply(pd.Series).drop(columns=0).astype(float)
    # rename the columns according to valid obos
    aop_df = aop_df.rename(columns=aop_df.columns.to_series().pipe(fix_obos))
    aop_df = aop_df.drop(columns=[None])
    # we have to merge across both columns and index, since columns can also be duplicated
    # TODO: if we need more than just the FIRST age of presentation, this will be the part to change:
    aop_df = aop_df.groupby(aop_df.index).min()
    aop_df = aop_df.groupby(aop_df.columns, axis=1).min()
    return aop_df


def clean_phenotype_age_of_presentation(df: pd.DataFrame):
    """Adds the first age of presentation columns for each phenotype
    @param df:
    @return:
    """
    aop_col = str(AnnotationHeader.AGE_OF_PRESENTATION)
    aop_df = df[aop_col].pipe(_separate_phenotype_age_of_presentation)
    col_map = {x: f'{str(AnnotationHeader.AGE_OF_PRESENTATION_CLEAN)}.{x}' for x in aop_df.columns}
    aop_df = aop_df.rename(columns=col_map)
    df = pd.concat([df, aop_df], axis=1)
    return df


def _separate_disease_entity(disease_series: pd.Series):
    """Given a series of HPO terms where each cell is a list, create a df where columns are phenotypes and each cell
    is a true or false for if the phenotype is present"""
    exploded_series = disease_series.explode().pipe(fix_obos)
    disease_df = pd.get_dummies(exploded_series.apply(pd.Series).stack(dropna=False))

    disease_df = disease_df.groupby(level=0).sum().astype(int)
    return disease_df


def clean_disease_entity(df: pd.DataFrame):
    """	Adds cleaned disease entity phenotypes to the df
    @param df:
    @return:
    """
    disease_col = str(AnnotationHeader.DISEASE_ENTITY)

    disease_df = df[disease_col].pipe(_separate_disease_entity)
    col_map = {x: f'{str(AnnotationHeader.DISEASE_ENTITY_CLEAN)}.{x}' for x in disease_df.columns}
    disease_df = disease_df.rename(columns=col_map)

    df = pd.concat([df, disease_df], axis=1)

    return df


def _separate_mutation_type(muttype_series: pd.Series):
    exploded = muttype_series.explode().str.split(';', n=1, expand=True)[0].pipe(fix_obos)
    muttype_df = pd.get_dummies(exploded.apply(pd.Series).stack(dropna=False))

    muttype_df = muttype_df.groupby(level=0).sum().astype(int)
    return muttype_df


def clean_mutation_type(df: pd.DataFrame):
    mutation_type_col = str(AnnotationHeader.MUTATION_TYPE)
    mutation_type_df = df[mutation_type_col].pipe(_separate_mutation_type)

    col_map = {x: f'{str(AnnotationHeader.MUTATION_TYPE_CLEAN)}.{x}' for x in mutation_type_df.columns}
    mutation_type_df = mutation_type_df.rename(columns=col_map)

    df = pd.concat([df, mutation_type_df], axis=1)
    return df


def clean_kindred_id(df: pd.DataFrame):
    kindred_colname = str(AnnotationHeader.KINDRED_ID)
    kindred_df = df[kindred_colname].explode()
    kindred_df = kindred_df.groupby(kindred_df.index).first()

    df[str(AnnotationHeader.KINDRED_ID_CLEAN)] = kindred_df
    return df


def _clean_vhl_columns(df: pd.DataFrame, column):
    # find all columns that start with clean_disease_col
    colname_series = df.columns
    disease_colnames = colname_series[colname_series.str.startswith(column)]
    clean_columns = df[disease_colnames]

    clean_vhl_columns = clean_columns.rename(lambda x: variant_functions.generalized_vhl_phenotype(x.split('.')[-1],
                                                                                                   return_as=config.OBO_RETURN_TYPE),
                                             axis='columns')
    return clean_vhl_columns


def add_generalized_disease_entity(df: pd.DataFrame):
    """
    Adds VHL-generalized disease entity phenotypes
    This has to be run after clean_disease_entity()
    @param df:
    @return:
    """
    clean_vhl_columns = _clean_vhl_columns(df, str(AnnotationHeader.DISEASE_ENTITY_CLEAN))

    clean_grouped_columns = clean_vhl_columns.groupby(by=clean_vhl_columns.columns, axis='columns').sum()
    clean_grouped_columns[clean_grouped_columns >= 1] = 1
    clean_grouped_columns = clean_grouped_columns.rename(
        lambda x: f'{AnnotationHeader.GENERALIZED_VHL_DISEASE_ENTITY}.{x}', axis=1)

    df[clean_grouped_columns.columns] = clean_grouped_columns

    return df


def add_generalized_phenotype_age_of_presentation(df: pd.DataFrame):
    """
    Adds VHL-generalized phenotypes to age of presentation
    This has to be run after clean_phenotype_age_of_presentation()
    @param df:
    @return:
    """
    clean_vhl_columns = _clean_vhl_columns(df, str(AnnotationHeader.AGE_OF_PRESENTATION_CLEAN))
    # using min here implies first age of presentation
    clean_grouped_columns = clean_vhl_columns.groupby(by=clean_vhl_columns.columns, axis='columns').min()
    clean_grouped_columns = clean_grouped_columns.rename(
        lambda x: f'{AnnotationHeader.GENERALIZED_VHL_AGE_OF_PRESENTATION}.{x}', axis=1)

    df[clean_grouped_columns.columns] = clean_grouped_columns

    return df


def clean_protein_position(df: pd.DataFrame):
    """
    Adds a cleaned protein position column to the dataframe
    @param df:
    @return:
    """
    protein_position_col = str(AnnotationHeader.PROTEIN_POSITION)
    pp_numeric = pd.to_numeric(df[protein_position_col].explode(), errors='coerce')
    cleaned_pp = pp_numeric.groupby(by=pp_numeric.index).first()
    df[str(AnnotationHeader.PROTEIN_POSITION_CLEAN)] = cleaned_pp

    return df


# TODO: bring up cases with multiple cdna / cdna with underlines
def clean_cdna_position(df: pd.DataFrame):
    """
    Adds a cleaned cdna column to the dataframe
    @param df:
    @return:
    """
    cdna_col = str(AnnotationHeader.CDNA_POSITION)
    cdna_numeric = pd.to_numeric(df[cdna_col].explode(), errors='coerce')
    cleaned_cdna = cdna_numeric.groupby(by=cdna_numeric.index).first()
    df[str(AnnotationHeader.CDNA_POSITION_CLEAN)] = cleaned_cdna

    return df

def add_region_columns(df: pd.DataFrame):
    """
    Adds a column corresponding to each function region, with values as 1 if the region was affected
    This must be run after clean_protein_position()
    @param df:
    @return:
    """
    protein_position = df[str(AnnotationHeader.PROTEIN_POSITION_CLEAN)]

    for region, codon in variant_functions.VHL_FUNCTIONAL_REGIONS.items():
        rows_in_region = protein_position.isin(codon)
        df.loc[:, f'{str(AnnotationHeader.FUNCTIONAL_REGION)}.{region}'] = 0
        df.loc[rows_in_region, f'{str(AnnotationHeader.FUNCTIONAL_REGION)}.{region}'] = 1

    return df


def _clean_so_columns(df: pd.DataFrame, column):
    # find all columns that start with column
    colname_series = df.columns
    disease_colnames = colname_series[colname_series.str.startswith(column)]
    clean_columns = df[disease_colnames]

    clean_so_columns = clean_columns.rename(lambda x: variant_functions.generalized_so_terms(x.split('.')[-1],
                                            return_as=config.OBO_RETURN_TYPE),
                                            axis='columns')
    return clean_so_columns


def add_generalized_mutation_type(df: pd.DataFrame):
    """
    Adds generalized sequence ontology terms
    This has to be run after clean_mutation_type()
    @param df:
    @return:
    """
    clean_so_columns = _clean_so_columns(df, str(AnnotationHeader.MUTATION_TYPE_CLEAN))

    clean_grouped_columns = clean_so_columns.groupby(by=clean_so_columns.columns, axis='columns').sum()
    clean_grouped_columns[clean_grouped_columns >= 1] = 1
    clean_grouped_columns = clean_grouped_columns.rename(
        lambda x: f'{AnnotationHeader.GENERALIZED_MUTATION_TYPE}.{x}', axis=1)

    df[clean_grouped_columns.columns] = clean_grouped_columns

    return df


def clean_civic_name_column(df: pd.DataFrame):
    """
    Adds a cleaned civic name column to the dataframe
    @param df:
    @return:
    """
    civic_name_col = str(AnnotationHeader.CIVIC_NAME)
    civic_name_string = df[civic_name_col].explode()
    cleaned_civic_name = civic_name_string.groupby(by=civic_name_string.index).first()
    df[str(AnnotationHeader.CIVIC_NAME_CLEAN)] = cleaned_civic_name

    return df


def add_exon_deletion_columns(df: pd.DataFrame):
    """
    Adds a column corresponding to each exon deletion type, with values as 1 if the variant is of this type.
    This has to be run after clean_civic_name_column()
    @param df:
    @return:
    """
    civic_name = df[str(AnnotationHeader.CIVIC_NAME_CLEAN)]
    
    for exon_deletion_type in variant_functions.EXON_DELETION_TERMS:
        rows_with_deletion = civic_name.eq(exon_deletion_type)
        df.loc[:, f'{str(AnnotationHeader.EXON_DELETION_TERM)}.{exon_deletion_type}'] = 0
        df.loc[rows_with_deletion, f'{str(AnnotationHeader.EXON_DELETION_TERM)}.{exon_deletion_type}'] = 1
        
    return df


def add_grouped_mutation_type(df: pd.DataFrame):
    """
    Adds truncating/non-truncating columns
    This has to be run after add_generalized_mutation_type()
    @param df:
    @return:
    """

    generalized_so_columns = _clean_so_columns(df, str(AnnotationHeader.GENERALIZED_MUTATION_TYPE))

    for mutgroup, so_type in variant_functions.SO_TERM_TYPES.items():
        rows_in_mutgroup = generalized_so_columns.iloc[:, generalized_so_columns.columns.isin(so_type)].any(axis=1)
        df.loc[:, f'{str(AnnotationHeader.GROUPED_MUTATION_TYPE)}.{mutgroup}'] = 0
        df.loc[rows_in_mutgroup, f'{str(AnnotationHeader.GROUPED_MUTATION_TYPE)}.{mutgroup}'] = 1
    return df


def preprocess(df):
    out_df = (
        df.copy()
        .pipe(fix_na)
        .pipe(clinvar_to_variant)
        .pipe(caid_to_variant)
        .pipe(clean_pmid)
        .pipe(clean_kindred_id)
        .pipe(clean_phenotype_age_of_presentation)
        .pipe(clean_disease_entity)
        .pipe(clean_mutation_type)
        .pipe(clean_cdna_position)
        .pipe(add_generalized_disease_entity)
        .pipe(add_generalized_phenotype_age_of_presentation)
        .pipe(clean_protein_position)
        .pipe(add_region_columns)
        .pipe(add_generalized_mutation_type)
        .pipe(clean_civic_name_column)
        .pipe(add_exon_deletion_columns)
        .pipe(add_grouped_mutation_type)

    )

    # columns in alphabetical order
    out_df = out_df.reindex(sorted(out_df.columns), axis=1)

    return out_df
