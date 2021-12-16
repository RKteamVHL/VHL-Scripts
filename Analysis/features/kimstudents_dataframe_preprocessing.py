import collections
import logging
import re

import numpy as np
import pandas as pd

from .. import variant_functions as vf
from ..constants import *

# this file contains all code relevant for cleaning up the raw KimStudents masterlist dataframe. generally, the strategy
# is to keep all of the input columns and append additional, analysis-specific columns
# The majority of the functions in this file are written to be used by the pandas.DataFrame.pipe function, and the rest
# are helper-functions used by the piping functions


# patterns for extracted the year/month age numbers from the Evaluated and Last Known Age columns
EVALUATED_AGE_REGEX = re.compile("E((?P<Y>[0-9]+)Y)?((?P<M>[0-9]+)M)?")
LASTKNOWN_AGE_REGEX = re.compile("lk((?P<Y>[0-9]+)Y)?((?P<M>[0-9]+)M)?")

# CDNA_REGEX = re.compile("c\.([0-9]+)[ATCG]")

# this dictionary keeps lists of all additional computed columns grouped together roughly by the type of
# analysis the columns will be used for
COMPUTED_COLUMNS = {
    "phenotype": [],
    "denovo": [],
    "generalized_phenotype": [],
    "mutant_type": [],
    "generalized_mutant_type": [],
    "cdna": [],
    "codon": [],
    "age": [],
    "sex": [],
    "resolution": [],
    "domain": [],
    "region": [],
    "cdna_change": [],
    "aa_change": [],
    "grouped_mutation_type": []
}


def _phenotype_string_to_list(x, generalize=True):
    """
    Converts a string of phenotype hpo terms separated by a comma or semicolon into a list
    @param x: string of phenotypes
    @param generalize: if True, the function takes the listed phenotype and converts it to a context-relevant phenotype
        based on the HPO tree. The conversion is usually an abstraction to a higher node on the HPO tree
    @return:
    """
    hpo_list = re.split('[;,]', x)
    output_categories = []
    type_ = "phenotype"
    if generalize:
        type_ = "generalized_" + type_

    for term in hpo_list:
        term = term.casefold().strip()
        try:

            if term not in NULL_TERMS:

                var_hpo = None
                if generalize:
                    var_hpo = vf.generalized_vhl_phenotype(term)
                else:
                    var_hpo = vf.get_valid_obo(term)
                output_categories.append(f"{type_}.{var_hpo}")
            elif term.casefold() is 'none':
                output_categories.append(f"{type_}.none")

        except ValueError as e:
            logging.getLogger(type_).warning(repr(e))
    return list(set(output_categories))


def add_generalized_phenotype_columns(df):
    """
    Adds generalized phenotype data to a copy of the inputted dataframe and returns it
    @param df:
    @return:
    """
    hpo_series = df['Phenotype'].apply(_phenotype_string_to_list, generalize=True)

    pheno_counts = hpo_series.apply(collections.Counter)
    pheno_featurized = pd.DataFrame.from_records(pheno_counts)
    COMPUTED_COLUMNS["generalized_phenotype"].extend(pheno_featurized.columns.to_list())
    df = df.join(pheno_featurized)
    return df


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


def add_generalized_mutant_type_columns(df):
    """
    Adds generalized mutation-type data to a copy of the inputted dataframe and returns it
    @param df:
    @return:
    """
    variant_series = df['Mutation Type'].apply(_mutant_string_to_list, generalize=True)

    mutant_type_counts = variant_series.apply(collections.Counter)
    mutants_featurized = pd.DataFrame.from_records(mutant_type_counts)
    COMPUTED_COLUMNS["generalized_mutant_type"].extend(mutants_featurized.columns.to_list())
    df = df.join(mutants_featurized)
    return df


def _age_to_records(age_str):
    """
    Converts an Age string into an Age dict, calculating total age of the patient at the same time
    @param age_str:
    @return:
    """
    record = {}
    if isinstance(age_str, str):
        e_match = EVALUATED_AGE_REGEX.search(age_str)
        lk_match = LASTKNOWN_AGE_REGEX.search(age_str)

        if lk_match is None:
            lk_match = e_match

        if e_match is not None:
            var = e_match.groupdict()
            months = 0
            if var['Y'] is not None:
                months += int(var['Y']) * 12
            if var['M'] is not None:
                months += int(var['M'])

            if not (var['Y'] is None and var['M'] is None):
                years = months / 12
                record["evaluated_age"] = years

        if lk_match is not None:
            var = lk_match.groupdict()
            months = 0
            if var['Y'] is not None:
                months += int(var['Y']) * 12
            if var['M'] is not None:
                months += int(var['M'])

            if not (var['Y'] is None and var['M'] is None):
                years = months / 12
                record["last_known_age"] = years

    return record


def add_age_columns(df):
    """
    Adds age data to a copy of the inputted dataframe and returns it
    @param df:
    @return:
    """

    age_series = df['Age'].apply(_age_to_records)

    age_featurized = pd.DataFrame(age_series.to_list())
    COMPUTED_COLUMNS["age"].extend(age_featurized.columns.to_list())
    df = df.join(age_featurized)
    return df


def _start_cdna_change(x):
    cdna_list = re.split('[;,]', x)
    output = {}
    for cdna in cdna_list:

        match = vf.get_cdna_start(cdna)

        if match is not None:
            output["cdna_start"] = match

    return output

def add_cdna_start_columns(df):

    cdna_series = df['Mutation Event c.DNA.'].apply(_start_cdna_change)

    cdna_featurized = pd.DataFrame(cdna_series.to_list())
    COMPUTED_COLUMNS["cdna"].extend(cdna_featurized.columns.to_list())
    df = df.join(cdna_featurized)
    return df

def add_codon_columns(df):

    codon_series = np.ceil(df['cdna_start']/3)
    df["codon_start"] = codon_series[(codon_series >= 0) & (codon_series <= 213)]
    COMPUTED_COLUMNS["codon"].append("codon_start")
    return df


def _sex_to_list(x):
    output_categories = []
    if x:
        output_categories.append(f"sex.{x.casefold()}")

    return output_categories

def add_sex_columns(df):
    sex_series = df['Sex'].apply(_sex_to_list)

    sexcounts = sex_series.apply(collections.Counter)
    sex_featurized = pd.DataFrame.from_records(sexcounts)
    COMPUTED_COLUMNS["sex"].extend(sex_featurized.columns.to_list())
    df = df.join(sex_featurized)
    return df

def add_resolution_columns(df):
    series = df['Resolution'].apply(lambda x: [f"resolution.{x.casefold()}"])

    counts = series.apply(collections.Counter)
    featurized = pd.DataFrame.from_records(counts)
    COMPUTED_COLUMNS["resolution"].extend(featurized.columns.to_list())
    df = df.join(featurized)
    return df

def _functional_regions_from_cdna(x):
    mutlist = re.split('[;]', x)

    domain_val = {}
    for term in mutlist:
        term = term.strip()
        if term not in NULL_TERMS:
            for domain in vf.affected_domains(term):
                col_name = f"region.{domain}"
                domain_val[col_name] = domain_val.get(col_name, 0) + 1

    return domain_val


def add_region_columns(df):
    series = df['Mutation Event c.DNA.'].apply(_functional_regions_from_cdna)

    featurized = pd.DataFrame(series.to_list())
    COMPUTED_COLUMNS["region"].extend(featurized.columns.to_list())
    COMPUTED_COLUMNS["domain"].extend([f'region.{dom}' for dom in vf.VHL_DOMAIN_NAMES])
    df = df.join(featurized)
    return df


def add_grouped_mutation_type_columns(df):
    for grouptype, muttype_list in vf.SO_TERM_TYPES.items():
        groupcol = f'grouped_mutation_type.{grouptype}'
        df[groupcol] = 0
        for mtype in muttype_list:
            df[groupcol] = df[groupcol] + df[f'generalized_mutant_type.{mtype}'].fillna(0)

    colnames = [f'grouped_mutation_type.{sotype}' for sotype in vf.SO_TERM_TYPES.keys()]
    COMPUTED_COLUMNS["grouped_mutation_type"].extend(colnames)

    return df




def add_aa_change_columns(df):
    series = df['Predicted Consequence Protein Change'].apply(vf.get_aa_from_predicted_consequence)
    series = series.apply(lambda x: f"aa_change.{x}")
    featurized = pd.get_dummies(series).drop(columns=['aa_change.None'])

    COMPUTED_COLUMNS["aa_change"].extend(featurized.columns.to_list())
    df = df.join(featurized)
    return df

def add_denovo_column(df):
    series = df["Confirmed De Novo"].str.strip()
    series = series.str.casefold()
    series = series.apply(lambda x: f"denovo.{x}")

    featurized = pd.get_dummies(series)
    COMPUTED_COLUMNS["denovo"].extend(featurized.columns.to_list())

    df = df.join(featurized)
    return df


def kimstudents_preprocessing(df):
    df = (df
          .pipe(add_generalized_phenotype_columns)
          # .pipe(add_phenotype_columns)
          .pipe(add_age_columns)
          .pipe(add_generalized_mutant_type_columns)
          .pipe(add_cdna_start_columns)
          .pipe(add_codon_columns)
          .pipe(add_sex_columns)
          .pipe(add_denovo_column)
          .pipe(add_resolution_columns)
          .pipe(add_region_columns)
          .pipe(add_grouped_mutation_type_columns)
          .pipe(add_aa_change_columns)
          )
    df["Resolution"] = df["Resolution"].str.casefold()
    return df



