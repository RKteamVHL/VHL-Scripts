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
	out_df = df[df['type'] == AnnotationType.CASE]
	# TODO: create index based on case, pmid, and kindred
	return out_df


def fix_na(df: pd.DataFrame):
	def _fix_na(_df):
		return_df = _df.replace(config.NULL_TERMS, [np.nan] * len(config.NULL_TERMS))
		return return_df

	out_df = df.pipe(_fix_na)
	return out_df


def fix_obos(series: pd.Series):
	out_series = series.apply(variant_functions.get_valid_obo, return_as='name_spaceless')
	pass
	return out_series


def _caid_to_variant(caid):
	variant_dict = next(caid_to_variant_generator)
	variant_name = None
	if caid in variant_dict:
		# caid is guaranteed to be in the dict, but this checks for None
		variant_name = variant_dict[caid]
	return variant_name


def _fix_caid(caid_list):
	if isinstance(caid_list, list):
		caid = caid_list[0].split(', ')[0]
		variant_dict = next(caid_to_variant_generator)
		caid_clean = None
		if caid in variant_dict:
			caid_clean = caid
		return caid_clean


def caid_to_variant(df: pd.DataFrame, replace=False):
	"""
	Fix the caid column in the given df, and append two columns: the cleaned caid and the caid variant
	@param df:
	@param replace: if true, replaces the original caid column
	@return:
	"""

	out_df = df.copy()

	_caid = str(AnnotationHeader.CA_ID)
	_caid_clean = str(AnnotationHeader.CA_ID_CLEAN)
	_caid_variant = str(AnnotationHeader.CA_ID_VARIANT)

	caid_id_series = out_df[_caid].apply(_fix_caid)
	caid_id_series.name = _caid_clean

	caid_variant_series = caid_id_series.apply(_caid_to_variant)
	caid_variant_series.name = _caid_variant

	out_df = pd.concat([out_df, caid_id_series, caid_variant_series], axis=1)
	if replace:
		out_df = out_df.drop(columns=_caid)
	return out_df


def _fix_clinvar(clinvar_list):
	clinvar_id = None
	if isinstance(clinvar_list, list):
		for item in clinvar_list:
			if str.isdigit(item):
				clinvar_id = int(item)
	return clinvar_id


def _clinvar_to_variant(clinvar_id):
	variant_dict = next(clinvarid_to_variant_generator)
	variant_name = None
	if clinvar_id in variant_dict:
		variant_name = variant_dict[clinvar_id]
	return variant_name


def clinvar_to_variant(df: pd.DataFrame, replace=False):
	"""
	Fix the clinvar column in the given df, and append two columns: the cleaned clinvar and the clinvar variant
	@param df:
	@param replace: if true, replaces the original clinvar column
	@return:
	"""

	_clinvar = str(AnnotationHeader.CLINVAR_ID)
	_clinvar_clean = str(AnnotationHeader.CLINVAR_ID_CLEAN)
	_clinvar_variant = str(AnnotationHeader.CLINVAR_ID_VARIANT)

	out_df = df.copy()
	clinvar_id_series = out_df[_clinvar].apply(_fix_clinvar)
	clinvar_id_series.name = _clinvar_clean

	clinvar_variant_series = clinvar_id_series.apply(_clinvar_to_variant)
	clinvar_variant_series.name = _clinvar_variant

	out_df = pd.concat([out_df, clinvar_id_series, clinvar_variant_series], axis=1)
	if replace:
		out_df = out_df.drop(columns=_clinvar)

	return out_df


def _separate_phenotype_age_of_presentation(aop_series: pd.Series):


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


def add_phenotype_age_of_presentation(df: pd.DataFrame, replace=False):
	"""
	Adds the first age of presentation columns for each phenotype
	@param df:
	@param replace: if true, replaces the original ageofpresentation column
	@return:
	"""
	out_df = df.copy()
	aop_col = str(AnnotationHeader.AGE_OF_PRESENTATION)
	aop_df = out_df[aop_col].pipe(_separate_phenotype_age_of_presentation)
	col_map = {x: f'{str(AnnotationHeader.AGE_OF_PRESENTATION_CLEAN)}.{x}' for x in aop_df.columns}
	aop_df = aop_df.rename(columns=col_map)
	out_df = pd.concat([out_df, aop_df], axis=1)

	return out_df


def _separate_disease_entity(disease_series: pd.Series):
	exploded_series = disease_series.explode().pipe(fix_obos)
	disease_df = pd.get_dummies(exploded_series.apply(pd.Series).stack(dropna=False))

	disease_df = disease_df.groupby(level=0).sum().astype(int)
	return disease_df

def add_disease_entity(df: pd.DataFrame, replace=False):
	"""
	Adds disease entity phenotypes
	@param df:
	@param replace: if true, replaces the original disease column
	@return:
	"""
	out_df = df.copy()
	disease_col = str(AnnotationHeader.DISEASE_ENTITY)

	disease_df = out_df[disease_col].pipe(_separate_disease_entity)
	col_map = {x: f'{str(AnnotationHeader.DISEASE_ENTITY_CLEAN)}.{x}' for x in disease_df.columns}
	disease_df = disease_df.rename(columns=col_map)

	out_df = pd.concat([out_df, disease_df], axis=1)

	return out_df


def preprocess(df):
	out_df = df\
		.pipe(fix_na)\
		.pipe(clinvar_to_variant)\
		.pipe(caid_to_variant)\
		.pipe(add_phenotype_age_of_presentation)\
		.pipe(add_disease_entity)

	# columns in alphabetical order
	out_df = out_df.reindex(sorted(out_df.columns), axis=1)

	return out_df
