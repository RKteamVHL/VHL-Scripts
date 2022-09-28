from typing import List
from ..fetching.Annotation import BODY_TAGS_NAME, TEXT_TAGS_NAME, AugmentedAnnotation, AnnotationType
from ..variant_functions import get_variant_by_caid, clinvarid_to_variant_dict, VHL_PHENOTYPES
import pandas as pd
import numpy as np

NULL_TERMS = ['n/a', 'N/A', 'NA']
NON_STANDARD_REFS = ['NM_000551.2', 'AF010238.1']
STANDARD_REFS = ['NM_000551.3', 'NM_000551.4']




def _fix_na(df):
    return_df = df.replace(NULL_TERMS, [np.nan]*len(NULL_TERMS))
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


# paper stats
def get_unique_clinvar_variants(annotation_df: pd.DataFrame):
    clinvar_col = f"{BODY_TAGS_NAME}.ClinVarID"
    caid_col = f'{BODY_TAGS_NAME}.CAID'
    civic_col = f'{BODY_TAGS_NAME}.CivicName'


    valid_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    clinvar_id_series = valid_df[clinvar_col].apply(_fix_clinvar)

    variants_df = pd.concat([valid_df[["type", caid_col, civic_col]], clinvar_id_series], axis=1)


    variants_df.loc[:, caid_col] = variants_df[caid_col].map(lambda x: x[0] if isinstance(x, list) else np.nan)
    variants_df.loc[:, civic_col] = variants_df[civic_col].map(lambda x: x[0] if isinstance(x, list) else np.nan)
    variants_df = variants_df.dropna(subset=[clinvar_col, caid_col, civic_col], how='all')

    clinvar_mapped_series = variants_df[clinvar_col].map(clinvarid_to_variant_dict())
    clinvar_mapped_series.name = "Variant"

    variants_df = pd.concat([variants_df, clinvar_mapped_series], axis=1)


    counts_df = variants_df.groupby([clinvar_col, "Variant"], as_index=False).count()
    return variants_df


def get_unique_caid_variants(annotation_df: pd.DataFrame):
    caid_col = f'{BODY_TAGS_NAME}.CAID'
    clinvar_col = f"{BODY_TAGS_NAME}.ClinVarID"

    valid_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    clinvar_id_series = valid_df[clinvar_col].apply(_fix_clinvar)
    caid_series = valid_df[caid_col].apply(_caid_to_variant)

    variants_df = pd.concat([valid_df["type"], clinvar_id_series, caid_series], axis=1)
    variants_df = variants_df.dropna(subset=[caid_col])


    counts_df = variants_df
    return counts_df

def get_unique_variants(annotation_df: pd.DataFrame):
    valid_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]
    variants_df = valid_df.dropna(subset=["Variant"])
    counts_df = variants_df.groupby('Variant', as_index=False).count()
    counts_df = counts_df[["Variant", "type"]]
    return counts_df


def get_papers(annotation_df: pd.DataFrame):
    pmid_col = f'{TEXT_TAGS_NAME}.PMID'
    type_col = "type"

    information_df = annotation_df[annotation_df["type"] == AnnotationType.INFORMATION.name]
    information_df = information_df.dropna(subset=[pmid_col])
    information_df.loc[:, pmid_col] = information_df[pmid_col].map(lambda x: x[0])
    information_df = information_df[[pmid_col, type_col]]
    counts_df = information_df.groupby(pmid_col, as_index=False).count()
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
    assay_col = f'{BODY_TAGS_NAME}.ExperimentalAssay'
    assay_df = annotation_df[annotation_df["type"].isin([AnnotationType.ASSAY.name])]
    assay_df = assay_df.dropna(subset=[assay_col])
    return assay_df

def get_unregistered_variants(annotation_df: pd.DataFrame):
    unreg_col = f'{BODY_TAGS_NAME}.UnregisteredVariant'
    var_col = f'{TEXT_TAGS_NAME}.Variant'
    valid_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    valid_df = valid_df.dropna(subset=[unreg_col])
    valid_df.loc[:, unreg_col] = valid_df[unreg_col].map(lambda x: x[0])
    valid_df.loc[:, var_col] = valid_df[var_col].map(lambda x: x[0])


    counts_df = valid_df[[var_col, unreg_col, "type"]]
    return counts_df

#Questions:
# right now this is only using the refseq column, which includes refseqs that have been assumed to be standard, but
# not outrighted stated
def get_nonstandard_refseq(annotation_df: pd.DataFrame):
    refseq_col = f'{TEXT_TAGS_NAME}.ArticleReferenceSequence'

    valid_df = annotation_df[annotation_df["type"].isin([AnnotationType.METHODOLOGY.name])]
    refseq_df = valid_df.dropna(subset=[refseq_col])
    refseq_df.loc[:, refseq_col] = refseq_df[refseq_col].map(lambda x: x[0])

    refseq_df.loc[:, "StandardRef"] = False
    refseq_df.loc[:, "NonStandardRef"] = False

    refseq_df.loc[:, "StandardRef"] = refseq_df[refseq_col].str.contains('|'.join(STANDARD_REFS))
    refseq_df.loc[:, "NonStandardRef"] = refseq_df[refseq_col].str.contains('|'.join(NON_STANDARD_REFS))

    nonstandard_df = refseq_df
    return nonstandard_df


def get_family_pedigree_variants(annotation_df: pd.DataFrame):
    ped_col = f'{BODY_TAGS_NAME}.FamilyPedigree'
    var_col = f'{TEXT_TAGS_NAME}.Variant'
    valid_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    valid_df = valid_df.dropna(subset=[ped_col])
    valid_df.loc[:, ped_col] = valid_df[ped_col].map(lambda x: x[0])
    valid_df.loc[:, var_col] = valid_df[var_col].map(lambda x: x[0])

    counts_df = valid_df[[var_col, ped_col, "type"]]
    return counts_df

def get_previously_published_variants(annotation_df: pd.DataFrame):
    pub_col = f'{TEXT_TAGS_NAME}.PreviouslyPublished'
    var_col = f'{TEXT_TAGS_NAME}.Variant'
    valid_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    valid_df = valid_df.dropna(subset=[pub_col])
    valid_df.loc[:, pub_col] = valid_df[pub_col].map(lambda x: x[0])
    valid_df.loc[:, var_col] = valid_df[var_col].map(lambda x: x[0])
    valid_df = valid_df.pipe(_fix_na)
    valid_df = valid_df.dropna(subset=[pub_col])


    counts_df = valid_df[[var_col, pub_col, "type"]]
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
    y_labels = list(set(VHL_PHENOTYPES.values()))

    muttype_col = f'{BODY_TAGS_NAME}.MutationType'

    # features
    aa_change_col = f'{BODY_TAGS_NAME}.AminoAcidChange'
    codon_col = f'{BODY_TAGS_NAME}.ProteinPosition'

    # labels
    pheno_col = f'{BODY_TAGS_NAME}.DiseaseEntity'

    valid_df = annotation_df[annotation_df["type"].isin([AnnotationType.COHORT.name, AnnotationType.CASE.name])]

    # clean up the features
    feature_df = valid_df[[muttype_col, aa_change_col, codon_col]].dropna().applymap(lambda x: x[0]).apply(_fix_na)
    feature_df = feature_df.dropna(how='any')
    feature_df = feature_df[feature_df[muttype_col].str.contains('missense_variant')]

    feature_df[['from_aa', 'to_aa']] = feature_df[aa_change_col].str.split('to', expand=True)
    feature_df[feature_df['to_aa'] == '*'] = 'TER'
    feature_df = feature_df.rename(columns={codon_col: 'pos'})

    # clean up the labels
    pheno_df = valid_df.loc[feature_df.index, [pheno_col]]

    def _fix_pheno(cell):
        to_return = []
        if isinstance(cell, list):
            for ele in cell:
                if isinstance(ele, str):
                    if ele.casefold() in VHL_PHENOTYPES:
                        to_return.append(VHL_PHENOTYPES[ele.casefold()])
        else:
            to_return = []

        if not to_return:
            to_return.append('asymptomatic')

        return to_return

    pheno_s = pheno_df[pheno_col].map(_fix_pheno)
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


def get_all_statistics(annotations: List[AugmentedAnnotation]):
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

        #df = df.reindex(sorted(df.columns), axis=1)

        df.name = fn.__name__.replace("get_", "")

        output_df_list.append(df)

    return raw_df, output_df_list
