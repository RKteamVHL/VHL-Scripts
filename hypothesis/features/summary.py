from typing import List
from ..annotations.Annotation import AugmentedAnnotation, AnnotationType, AnnotationHeader
from ..fetching.clinvar_variants import clinvarid_to_variant_dict
from ..fetching.caid_variants import get_variant_by_caid
from ..variant_functions import DISEASE_ENTITY_TO_HPO
from .. import config
import pandas as pd
import numpy as np
import os

def _fix_na(df):
    return_df = df.replace(config.NULL_TERMS, [np.nan]*len(config.NULL_TERMS))
    return return_df


def _fix_clinvar(clinvar_list):
    clinvar_id = None
    if isinstance(clinvar_list, list):
        for item in clinvar_list:
            if str.isdigit(item):
                clinvar_id = int(item)
    return clinvar_id


def _caid_to_variant(caid):
    variant_name = None
    variant = get_variant_by_caid(caid)
    if variant is not None:
        variant_name = variant['communityStandardTitle'][0]
    return variant_name

# aliases to the header types in annotation for brevity
_clinvar = str(AnnotationHeader.CLINVAR)
_caid = str(AnnotationHeader.CAID)
_civic = str(AnnotationHeader.CIVICNAME)
_type = str(AnnotationHeader.TYPE)
_variant = str(AnnotationHeader.VARIANT)
_pmid = str(AnnotationHeader.PMID)
_assay = str(AnnotationHeader.EXPERIMENTALASSAY)
_unreg = str(AnnotationHeader.UNREGISTEREDVARIANT)
_refseq = str(AnnotationHeader.ARTICLEREFERENCESEQUENCE)
_ped = str(AnnotationHeader.FAMILYPEDIGREE)
_pub = str(AnnotationHeader.PREVIOUSLYPUBLISHED)
_mut = str(AnnotationHeader.MUTATIONTYPE)
_aa = str(AnnotationHeader.AMINOACIDCHANGE)
_pp = str(AnnotationHeader.PROTEINPOSITION)
_pheno = str(AnnotationHeader.DISEASEENTITY)

# paper stats
def get_unique_clinvar_variants(annotation_df: pd.DataFrame):

    valid_df = annotation_df[annotation_df[_type].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    clinvar_id_series = valid_df[_clinvar].apply(_fix_clinvar)

    variants_df = pd.concat([valid_df[[_type, _caid, _civic]], clinvar_id_series], axis=1)


    variants_df.loc[:, _caid] = variants_df[_caid].map(lambda x: x[0] if isinstance(x, list) else np.nan)
    variants_df.loc[:, _civic] = variants_df[_civic].map(lambda x: x[0] if isinstance(x, list) else np.nan)
    variants_df = variants_df.dropna(subset=[_clinvar, _caid, _civic], how='all')

    clinvar_mapped_series = variants_df[_clinvar].map(clinvarid_to_variant_dict())
    clinvar_mapped_series.name = _variant

    variants_df = pd.concat([variants_df, clinvar_mapped_series], axis=1)

    counts_df = variants_df.groupby([_clinvar, _variant], as_index=False).count()
    return variants_df


def get_unique_caid_variants(annotation_df: pd.DataFrame):

    valid_df = annotation_df[annotation_df[_type].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    clinvar_id_series = valid_df[_clinvar].apply(_fix_clinvar)
    caid_series = valid_df[_caid].apply(_caid_to_variant)

    variants_df = pd.concat([valid_df[_type], clinvar_id_series, caid_series], axis=1)
    variants_df = variants_df.dropna(subset=[_caid])

    counts_df = variants_df
    return counts_df


def get_unique_variants(annotation_df: pd.DataFrame):

    valid_df = annotation_df[annotation_df[_type].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]
    variants_df = valid_df.dropna(subset=[_variant])
    counts_df = variants_df.groupby(_variant, as_index=False).count()
    counts_df = counts_df[[_variant, _type]]
    return counts_df


def get_papers(annotation_df: pd.DataFrame):
    information_df = annotation_df[annotation_df[_type] == AnnotationType.INFORMATION.name]
    information_df = information_df.dropna(subset=[_pmid])
    information_df.loc[:, _pmid] = information_df[_pmid].map(lambda x: x[0])
    information_df = information_df[[_pmid, _type]]
    counts_df = information_df.groupby(_pmid, as_index=False).count()
    return counts_df

# #Questions:
# # Right now, this is finding all unique variants then dropping the ones that dont have a clinvar id- is this right?
# def get_clinvar_variants(annotation_df: pd.DataFrame):
#     variant_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]
#     clinvar_df = variant_df.dropna(subset=["ClinVarID"])
#     counts_df = clinvar_df.groupby('Variant', as_index=False).count()
#     counts_df = counts_df[["Variant", "type"]]
#     return counts_df


# Questions:
# Do we count just the evidence / assay tags, or count all instances of ExperimentalAssay?
def get_experimental_assays(annotation_df: pd.DataFrame):
    assay_df = annotation_df[annotation_df[_type].isin([AnnotationType.ASSAY.name])]
    assay_df = assay_df.dropna(subset=[_assay])
    return assay_df


def get_unregistered_variants(annotation_df: pd.DataFrame):
    valid_df = annotation_df[annotation_df[_type].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    valid_df = valid_df.dropna(subset=[_unreg])
    valid_df.loc[:, _unreg] = valid_df[_unreg].map(lambda x: x[0])
    valid_df.loc[:, _variant] = valid_df[_variant].map(lambda x: x[0])


    counts_df = valid_df[[_variant, _unreg, _type]]
    return counts_df

#Questions:
# right now this is only using the refseq column, which includes refseqs that have been assumed to be standard, but
# not outrighted stated
def get_nonstandard_refseq(annotation_df: pd.DataFrame):

    valid_df = annotation_df[annotation_df[_type].isin([AnnotationType.METHODOLOGY.name])]
    refseq_df = valid_df.dropna(subset=[_refseq])
    refseq_df.loc[:, _refseq] = refseq_df[_refseq].map(lambda x: x[0])

    refseq_df.loc[:, "StandardRef"] = False
    refseq_df.loc[:, "NonStandardRef"] = False

    refseq_df.loc[:, "StandardRef"] = refseq_df[_refseq].str.contains('|'.join(config.STANDARD_REFS))
    refseq_df.loc[:, "NonStandardRef"] = refseq_df[_refseq].str.contains('|'.join(config.NON_STANDARD_REFS))

    nonstandard_df = refseq_df
    return nonstandard_df


def get_family_pedigree_variants(annotation_df: pd.DataFrame):

    valid_df = annotation_df[annotation_df[_type].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    valid_df = valid_df.dropna(subset=[_ped])
    valid_df.loc[:, _ped] = valid_df[_ped].map(lambda x: x[0])
    valid_df.loc[:, _variant] = valid_df[_variant].map(lambda x: x[0])

    counts_df = valid_df[[_variant, _ped, _type]]
    return counts_df

def get_previously_published_variants(annotation_df: pd.DataFrame):
    valid_df = annotation_df[annotation_df[_type].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    valid_df = valid_df.dropna(subset=[_pub])
    valid_df.loc[:, _pub] = valid_df[_pub].map(lambda x: x[0])
    valid_df.loc[:, _variant] = valid_df[_variant].map(lambda x: x[0])
    valid_df = valid_df.pipe(_fix_na)
    valid_df = valid_df.dropna(subset=[_pub])


    counts_df = valid_df[[_variant, _pub, _type]]
    return counts_df

# harder metrics


def get_problem_variants(annotation_df: pd.DataFrame):
    pass

# maybe use PMIDs?
def get_previously_published(annotation_df: pd.DataFrame):
    pass

# what do we have tag-wise that can be used for determining tumor vs case report?
def get_case_tumors(annotation_df: pd.DataFrame):
    pass

# this can be done using AgeOfPresentation
def get_penetrance(annotation_df: pd.DataFrame):
    pass

def get_missense_variants(annotation_df: pd.DataFrame):
    x_labels = ['from_aa', 'to_aa', 'pos']
    y_labels = list(set(DISEASE_ENTITY_TO_HPO.values()))

    valid_df = annotation_df[annotation_df[_type].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    # clean up the features
    feature_df = valid_df[[_mut, _aa, _pp]].dropna().applymap(lambda x: x[0]).apply(_fix_na)
    feature_df = feature_df.dropna(how='any')
    feature_df = feature_df[feature_df[_mut].str.contains('missense_variant')]

    feature_df[['from_aa', 'to_aa']] = feature_df[_aa].str.split('to', expand=True)
    feature_df[feature_df['to_aa'] == '*'] = 'TER'
    feature_df = feature_df.rename(columns={_pp: 'pos'})

    # clean up the labels
    pheno_df = valid_df.loc[feature_df.index, [_pheno]]

    def _fix_pheno(cell):
        to_return = []
        if isinstance(cell, list):
            for ele in cell:
                if isinstance(ele, str):
                    if ele.casefold() in DISEASE_ENTITY_TO_HPO:
                        to_return.append(DISEASE_ENTITY_TO_HPO[ele.casefold()])
        else:
            to_return = []

        if not to_return:
            to_return.append('asymptomatic')

        return to_return

    pheno_s = pheno_df[_pheno].map(_fix_pheno)
    pheno_s = pheno_s.rename('phenotype')
    # s = pheno_s.explode()

    x = feature_df[x_labels]
    #y = (pd.crosstab(s.index, s)).clip(0, 1)
    y = pheno_s

    return x.join(y)


# these are statistics related to error-checking
def get_invalid_dataframe(annotations: List[AugmentedAnnotation]):
    dict_list = []
    for annotation in annotations:
        if annotation.type == AnnotationType.INVALID.name:
            dict_list.append({
                "link": annotation.links["html"],
                "type": annotation.type,
                "author": annotation.user,
            })
    out_df = pd.DataFrame.from_records(dict_list)
    return out_df


def get_annotation_summary(annotations: List[AugmentedAnnotation]):
    record_list = []
    for annotation in annotations:
        record = {"type": annotation.type}
        record.update({k: v[0] for k, v in annotation.tag_dictionary.items()})
        record_list.append(record)

    df = pd.DataFrame.from_records(record_list)
    df["count"] = 1
    out_df = df.groupby("type").count()
    return out_df["count"]


def get_all_summaries(annotations: List[AugmentedAnnotation]):
    output_df_list = []
    statistic_fns = [
        get_missense_variants,
        get_unique_clinvar_variants,
        get_previously_published_variants,
        get_papers,
        get_experimental_assays,
        get_unregistered_variants,
        get_nonstandard_refseq,
        get_family_pedigree_variants
    ]
    raw_df = AugmentedAnnotation.df_from_annotations(annotations)
    for fn in statistic_fns:
        df = raw_df.pipe(fn)
        df = df.dropna(axis="columns", how="all")

        df.name = fn.__name__.replace("get_", "")

        output_df_list.append(df)

    # columns in alphabetical order
    raw_df = raw_df.reindex(sorted(raw_df.columns), axis=1)
    raw_df.to_csv(os.path.join(config.DIRS['output'], config.RAW_ANNOTATION_DF))

    for stat in output_df_list:
        stat.to_csv(os.path.join(config.DIRS['summary'], f"{stat.name}.csv"))

    return raw_df, output_df_list
