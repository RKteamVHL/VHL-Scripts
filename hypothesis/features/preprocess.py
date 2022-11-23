from typing import List
from ..annotations.Annotation import AnnotationHeader
from ..fetching.caid_variants import caid_to_variant_generator
from ..fetching.clinvar_variants import clinvarid_to_variant_generator

from .. import config
import pandas as pd
import numpy as np

import re

from ..annotations.Annotation import AnnotationType


# TODO: this file will become very large, and should probably separated into a python module at some point

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

	_caid = str(AnnotationHeader.CAID)
	_caid_clean = str(AnnotationHeader.CAID_CLEAN)
	_caid_variant = str(AnnotationHeader.CAID_VARIANT)

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
	Fix the clinvar column in the given df
	@param df:
	@param replace: if true, replaces the original clinvar column
	@return:
	"""

	_clinvar = str(AnnotationHeader.CLINVAR)
	_clinvar_clean = str(AnnotationHeader.CLINVAR_CLEAN)
	_clinvar_variant = str(AnnotationHeader.CLINVAR_VARIANT)

	out_df = df.copy()
	clinvar_id_series = out_df[_clinvar].apply(_fix_clinvar)
	clinvar_id_series.name = _clinvar_clean

	clinvar_variant_series = clinvar_id_series.apply(_clinvar_to_variant)
	clinvar_variant_series.name = _clinvar_variant

	out_df = pd.concat([out_df, clinvar_id_series, clinvar_variant_series], axis=1)
	if replace:
		out_df = out_df.drop(columns=_clinvar)

	return out_df


def fix_singular_lists(df: pd.DataFrame, cols: List[str] = None):
	"""
	Fix cells that have lists in the input annotation df
	This function operates on columns that only (should) have a single value in each column. If there are
	@param df: Dataframe to fix
	@param cols: Columns to operate on. If not specified, all columns will be selected
	@return:

	"""

	if not cols:
		cols = list(df.columns)

	for col in cols:
		s = df[col]
		s_unwrapped = s.apply(lambda x: x[0])


def fix_df_lists(df: pd.DataFrame, cols: List[str] = None):
	"""
	Fix cells that have lists in the input annotation df

	This function replaces lists with dummy encoded columns based on what's in the list; a separate column will be
	created for each string in the list. List of dicts will create a separate column for each unique key in the column.
	Columns that have mostly boolean or numeric values will be unchanged
	@param df: Dataframe to fix
	@param cols: Columns to operate on. If not specified, all columns will be fixed
	@return:

	"""

	if not cols:
		cols = list(df.columns)

	for col in cols:
		series = df[col]


def preprocess(df):
	out_df = df\
		.pipe(fix_na)\
		.pipe(clinvar_to_variant)\
		.pipe(caid_to_variant)

	# columns in alphabetical order
	out_df = out_df.reindex(sorted(out_df.columns), axis=1)

	return out_df
