# the 'annotations' import (used for type-hinting in static methods)
# should not be confused with our custom Annotation class
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from collections import deque
from typing import Dict, List
from .. import config
import numpy as np
import pandas as pd
import re
import copy
import enum

NULL_TERMS = ["unknown", "none", "", "N/A"]

# BODY, TEXT, and COMPUTED need to be separated so that we can discern where the tags came from. It might make sense to
# separate these even more in the future
BODY_TAGS_NAME = "BODY"
TEXT_TAGS_NAME = "TEXT"
COMPUTED_TAGS_NAME = "COMPUTED"

# TODO: discuss separating annotation headers and types into separate file
# TODO: consider listing expected cell types in this enum


class AnnotationHeader(enum.Enum):
    """
    The body and text annotation headers are taken directly from the 2022 Hypothes.is VHL Annotation Protocol, and are in
    the same order as Table 1 in that document
    These were made an enum (as compared to a dict or list) for re-usability and type-hinting in IDEs.
    If a computed column is needed / being created, it should be put here, and ONLY this enum should ever be used
    while referencing columns dataframes later on
    """

    ## properties of the hypothesis annotations, not specific to the VHL annotations
    TYPE = ('type',)
    URI = ("uri",)
    SOURCE = ("source",)
    TEXT = ("text",)

    ## article information
    PMID = (TEXT_TAGS_NAME, "PMID")
    GENE = (TEXT_TAGS_NAME, "Gene")
    STANDARDIZED_REFERENCE_SEQUENCE = (TEXT_TAGS_NAME, "StandardizedReferenceSequence")

    ## methodology
    ARTICLE_REFERENCE_SEQUENCE = (TEXT_TAGS_NAME, "ArticleReferenceSequence")
    GENOTYPING_METHOD = (TEXT_TAGS_NAME, "GenotypingMethod")
    SAMPLING_METHOD = (TEXT_TAGS_NAME, "SamplingMethod")

    ## case-individual
    PATIENT_ID = (TEXT_TAGS_NAME, "PatientID")
    DISEASE_ASSERTION = (TEXT_TAGS_NAME, "DiseaseAssertion")
    FAMILY_INFO = (TEXT_TAGS_NAME, "FamilyInfo")
    CASE_PRESENTING_HPOS = (TEXT_TAGS_NAME, "CasePresentingHPOs")
    CASE_HPO_FREE_TEXT = (TEXT_TAGS_NAME, "CaseHPOFreeText")
    CASE_NOT_HPOS = (TEXT_TAGS_NAME, "CaseNotHPOs")
    CASE_NOT_HPO_FREE_TEXT = (TEXT_TAGS_NAME, "CaseNotHPOFreeText")
    CASE_PREVIOUS_TESTING = (TEXT_TAGS_NAME, "CasePreviousTesting")
    PREVIOUSLY_PUBLISHED = (TEXT_TAGS_NAME, "PreviouslyPublished")
    SUPPLEMENTAL_DATA = (TEXT_TAGS_NAME, "SupplementalData")
    VARIANT = (TEXT_TAGS_NAME, "Variant")
    LEGACY_VARIANT = (TEXT_TAGS_NAME, "LegacyVariant")
    CASE_PROBLEM_VARIANT_FREE_TEXT = (TEXT_TAGS_NAME, "CaseProblemVariantFreeText")
    CLINVAR_ID = (BODY_TAGS_NAME, "ClinVarID")  # NOTE: there is also a text tag for ClinVarID
    CA_ID = (BODY_TAGS_NAME, "CAID")    # NOTE: there is also a text tag for CAID
    GNOMAD = (TEXT_TAGS_NAME, "gnomAD")
    VARIANT_EVIDENCE = (TEXT_TAGS_NAME, "VariantEvidence")
    MUTATION_TYPE = (BODY_TAGS_NAME, "MutationType")    # NOTE: there is also a text tag for mutation type
    CIVIC_NAME = (BODY_TAGS_NAME, "CivicName")  # NOTE: there is also a text tag for civic name
    MULTIPLE_GENE_VARIANTS = (BODY_TAGS_NAME, "MultipleGeneVariants")   # NOTE: there is also a text tag for civic name

    ## group report- some of the group annotation columns are also in case reports
    GROUP_ID = (TEXT_TAGS_NAME, "GroupID")
    KINDRED_ID = (TEXT_TAGS_NAME, "KindredID")
    GROUP_SIZE = (TEXT_TAGS_NAME, "GroupSize")
    GROUP_NUMBER_MALES = (TEXT_TAGS_NAME, "GroupNumberMales")
    GROUP_NUMBER_FEMALES = (TEXT_TAGS_NAME, "GroupNumberFemales")
    GROUP_NUMBER_ETHNICITY = (TEXT_TAGS_NAME, "GroupNumberEthnicity")
    GROUP_AGE_RANGE = (TEXT_TAGS_NAME, "GroupAgeRange")
    GROUP_INFO = (TEXT_TAGS_NAME, "GroupInfo")
    GROUP_PRESENTING_HPOS = (TEXT_TAGS_NAME, "GroupPresentingHPOs")
    GROUP_HPO_FREE_TEXT = (TEXT_TAGS_NAME, "GroupHPOFreeText")
    GROUP_NOT_HPOS = (TEXT_TAGS_NAME, "GroupNotHPOs")
    GROUP_NOT_HPO_FREE_TEXT = (TEXT_TAGS_NAME, "GroupNotHPOFreeText")
    GROUP_PREVIOUS_TESTING = (TEXT_TAGS_NAME, "GroupPreviousTesting")

    ## case-report annotation tags
    # inheritance pattern
    INHERITANCE_PATTERN = (BODY_TAGS_NAME, "InheritancePattern")
    # disease entity
    DISEASE_ENTITY = (BODY_TAGS_NAME, "DiseaseEntity")
    # age of presentation
    AGE_OF_PRESENTATION = (BODY_TAGS_NAME, "AgeOfPresentation")
    # mutation
    MUTATION = (BODY_TAGS_NAME, "Mutation")
    # variant standardization
    REF_SEQ = (BODY_TAGS_NAME, "RefSeq")
    ASSUMED_REF_SEQ = (BODY_TAGS_NAME, "AssumedRefSeq")
    PROBLEM_VARIANT = (BODY_TAGS_NAME, "ProblemVariant")
    # variant
    CDNA_POSITION = (BODY_TAGS_NAME, "cDNAposition")
    PROTEIN_POSITION = (BODY_TAGS_NAME, "ProteinPosition")
    AMINO_ACID_CHANGE = (BODY_TAGS_NAME, "AminoAcidChange")
    UNREGISTERED_VARIANT = (BODY_TAGS_NAME, "UnregisteredVariant")
    MULTIPLE_VHL_VARIANTS = (BODY_TAGS_NAME, "MultipleVHLVariants")
    # mutation type, civic name, already covered above
    # previous testing
    PREVIOUS_TESTING = (BODY_TAGS_NAME, "PreviousTesting")
    # family information
    FAMILY_PEDIGREE = (BODY_TAGS_NAME, "FamilyPedigree")
    FAMILIAL = (BODY_TAGS_NAME, "Familial")
    NON_FAMILIAL = (BODY_TAGS_NAME, "NonFamilial")
    NO_FAMILY_INFO = (BODY_TAGS_NAME, "NoFamilyInfo")
    DE_NOVO_COMFIRMED = (BODY_TAGS_NAME, "deNovoConfirmed")
    COMPOUND_HETEROZYGOUS = (BODY_TAGS_NAME, "CompoundHeterozygous")
    FAMILY_COHORT = (BODY_TAGS_NAME, "FamilyCohort")
    # supplemental data covered above

    ## experimental assay
    EXPERIMENTAL_ASSAY = (BODY_TAGS_NAME, "ExperimentalAssay")
    # clinvar and caid are already covered above
    # TODO: CAiD in the document is inconsistently capitalized

    ## Evidence statement
    EVIDENCE_STATEMENT = (BODY_TAGS_NAME, "EvidenceStatement")
    # clinvarid, caid, civicname, and problem variant are all listed above

    ## computed columns
    # any columns, new or derrived from the above, should be put here
    PMID_CLEAN = (COMPUTED_TAGS_NAME, "PMID")
    CLINVAR_ID_CLEAN = (COMPUTED_TAGS_NAME, "ClinVarID")
    CA_ID_CLEAN = (COMPUTED_TAGS_NAME, "CAid")
    CLINVAR_ID_VARIANT = (COMPUTED_TAGS_NAME, "ClinVarIDVariant")
    CA_ID_VARIANT = (COMPUTED_TAGS_NAME, "CAidVariant")

    PROTEIN_POSITION_CLEAN = (COMPUTED_TAGS_NAME, "ProteinPosition")
    CDNA_POSITION_CLEAN = (COMPUTED_TAGS_NAME, "cDNAposition")

    AGE_OF_PRESENTATION_CLEAN = (COMPUTED_TAGS_NAME, "AgeOfPresentation")
    DISEASE_ENTITY_CLEAN = (COMPUTED_TAGS_NAME, "DiseaseEntity")
    KINDRED_ID_CLEAN = (COMPUTED_TAGS_NAME, "KindredID")
    MUTATION_TYPE_CLEAN = (COMPUTED_TAGS_NAME, "MutationType")

    GENERALIZED_VHL_AGE_OF_PRESENTATION = (COMPUTED_TAGS_NAME, "GeneralizedVHLAgeOfPresentation")
    GENERALIZED_VHL_DISEASE_ENTITY = (COMPUTED_TAGS_NAME, "GeneralizedVHLDiseaseEntity")
    GENERALIZED_MUTATION_TYPE = (COMPUTED_TAGS_NAME, "GeneralizedMutationType")
    GROUPED_MUTATION_TYPE = (COMPUTED_TAGS_NAME, "GroupedMutationType")

    FUNCTIONAL_REGION = (COMPUTED_TAGS_NAME, "FunctionalRegion")

    def __str__(self):
        """
        Function called when converting this enum to a sting (i.e., when called with str()).
        This functionality HAS to be included in a method rather than a constant,
        because the need for capitalization is determined at runtime (via config.CASEFOLD_TAG_NAMES)
        """
        # if the enum tuple only has one string i.e.,  type, uri, source, text
        if len(self.value) == 1:
            return self.value[0].casefold() if config.CASEFOLD_TAG_NAMES else self.value[0]
        else:
            a_list = [x.casefold() if config.CASEFOLD_TAG_NAMES else x for x in self.value[1:]]
            return '.'.join([self.value[0], *a_list])


class AnnotationType(enum.IntEnum):
    INVALID = -1
    REPLY = enum.auto()
    EVIDENCE = enum.auto()
    INFORMATION = enum.auto()
    METHODOLOGY = enum.auto()
    CASE = enum.auto()
    COHORT = enum.auto()
    ASSAY = enum.auto()


BODY_TAG_REGEX = re.compile("^(?P<name>\w+): *(?P<body>.*)$", flags=re.MULTILINE)

# for tag lists, only mandatory tags are used to determine the annotation type

EVIDENCE_TAGS = [
    AnnotationHeader.EVIDENCE_STATEMENT,
]

INFORMATION_TAGS = [
    AnnotationHeader.PMID,
    AnnotationHeader.GENE,
    AnnotationHeader.STANDARDIZED_REFERENCE_SEQUENCE
]

METHODOLOGY_TAGS = [
    AnnotationHeader.ARTICLE_REFERENCE_SEQUENCE,
    AnnotationHeader.GENOTYPING_METHOD,
    AnnotationHeader.SAMPLING_METHOD
]

CASE_TAGS = [
    AnnotationHeader.CASE_PRESENTING_HPOS
]

COHORT_TAGS = [
    AnnotationHeader.GROUP_PRESENTING_HPOS
]

ASSAY_TAGS = [
    AnnotationHeader.EXPERIMENTAL_ASSAY
]


def _fix_df_nan(df: pd.DataFrame):
    out_df = df
    for term in NULL_TERMS:
        out_df = out_df.replace(term, np.NaN)
    return out_df

@dataclass
class HypothesisAnnotation:
    """The base class for a hypothesis annotation

    This class only holds the data that a hypothesis annotation holds, with no additional functionality
    """
    id: str
    created: str
    updated: str
    user: str
    uri: str
    text: str
    tags: List[str]
    group: str
    permissions: Dict[str, List[str]]
    target: List[Dict[str, List]]
    document: Dict[str, List[str]]
    links: Dict[str, str]
    flagged: bool
    hidden: bool
    user_info: Dict[str, str]
    references: List[str] = field(default_factory=list)

    @staticmethod
    def from_dict(d):
        new_annotation = HypothesisAnnotation(**copy.deepcopy(d))
        return new_annotation


@dataclass
class AugmentedAnnotation(HypothesisAnnotation):
    """Augmented Hypothesis Annotation class

    This class contains additional properties and methods for parsing and storing tags embedded into the annotation
    body, along with converting all tags to a dictionary (rather than a string)
    """
    text_tags: Dict[str, str] = field(default_factory=dict)
    body_tags: Dict[str, str] = field(default_factory=dict)
    type: str = AnnotationType.INVALID.name

    def _get_tags_from_text(self, casefold_tag_names=True):
        tag_match = re.finditer(BODY_TAG_REGEX, self.text)
        self.text_tags = {}
        if tag_match:
            # TODO: some assertion/error checking on tag_dict
            for match in tag_match:
                tag_dict = match.groupdict()
                tag_name = tag_dict["name"]
                tag_value = tag_dict["body"].strip()

                # if we want the case folded
                if casefold_tag_names:
                    tag_name = tag_name.casefold()
                # if the tag name is already in the tag_dict, fetch its corresponding list; else, return an empty list
                new_list = self.text_tags.get(tag_name, [])

                # if the value is already in the tag_dict, ignore it (removes tags that get duplicated in body
                # and tag list)
                if tag_value not in new_list:
                    new_list.append(tag_value)
                self.text_tags[tag_name] = new_list
        else:  # error or log something here
            pass

    def _get_tags_from_tags_list(self, casefold_tag_names=True):
        self.body_tags = {}
        for tag in self.tags:
            tag_split = tag.split(":")
            tag_name = ""
            tag_value = ""

            # double tags, in the format TagName1:TagName2:TagValue i.e., AgeOfPresentation
            if len(tag_split) == 3:
                tag_name = tag_split[0].strip()
                tag_name2 = tag_split[1].strip()

                # if we want the case folded
                if casefold_tag_names:
                    tag_name2 = tag_name2.casefold()
                tag_value = {tag_name2: tag_split[2].strip()}

            # normal tags in the format TagName:TagValue
            elif len(tag_split) == 2:
                tag_name = tag_split[0].strip()
                tag_value = tag_split[1].strip()

            # flag tags, in the format TagName
            elif len(tag_split) == 1:
                tag_name = tag_split[0].strip()
                tag_value = True

            # if we want the case folded
            if casefold_tag_names:
                tag_name = tag_name.casefold()

            # if the tag name is already in the tag_dict, fetch its corresponding list; else, return an empty list
            new_list = self.body_tags.get(tag_name, [])

            # if the value is already in the tag_dict, ignore it (removes tags that get duplicated in body
            # and tag list)
            if tag_value not in new_list:
                new_list.append(tag_value)
            self.body_tags[tag_name] = new_list

    def _assign_type(self):
        if len(self.references) > 0:
            self.type = AnnotationType.REPLY.name
        # TODO: refactor code duplication
        # checking for evidence statement annotation
        elif any([str(key).split('.')[-1] in self.body_tags for key in EVIDENCE_TAGS]):
            self.type = AnnotationType.EVIDENCE.name

        # checking for article info annotation
        elif any([str(key).split('.')[-1] in self.text_tags for key in INFORMATION_TAGS]):
            self.type = AnnotationType.INFORMATION.name

        # checking for methodology annotation
        elif any([str(key).split('.')[-1] in self.text_tags for key in METHODOLOGY_TAGS]):
            self.type = AnnotationType.METHODOLOGY.name

        # checking for case annotation
        elif any([str(key).split('.')[-1] in self.text_tags for key in CASE_TAGS]):
            self.type = AnnotationType.CASE.name

        # checking for COHORT annotation
        elif any([str(key).split('.')[-1] in self.text_tags for key in COHORT_TAGS]):
            self.type = AnnotationType.COHORT.name

        # checking for experiment assay annotation
        elif any([str(key).split('.')[-1] in self.body_tags for key in ASSAY_TAGS]):
            self.type = AnnotationType.ASSAY.name

    def as_dict(self):
        return asdict(self)

    @staticmethod
    def from_dict(d, casefold_tag_names=True):
        new_annotation = AugmentedAnnotation(**copy.deepcopy(d))
        new_annotation._get_tags_from_text(casefold_tag_names)
        new_annotation._get_tags_from_tags_list(casefold_tag_names)
        new_annotation._assign_type()
        return new_annotation

    @staticmethod
    def df_from_annotations(annotations: List[AugmentedAnnotation]):
        """
        this static function returns a dataframe with a couple of caveats:
        1.    if the tag came from the body, it will be prepended with the BODY_TAGS_NAME constant,
              otherwise if it comes from the text, it will be prepended with TEXT_TAGS_NAME
        2.    individual cells will contain lists, which won't get formatted properly if exported to a csv file
        3.    some data (i.e., user info) is excluded from the output csv file
        """
        record_list = []
        column_set = set()
        for annotation in annotations:
            record = {
                str(AnnotationHeader.TYPE): annotation.type,
                str(AnnotationHeader.URI): annotation.links['html'],
                str(AnnotationHeader.SOURCE): annotation.target[0]['source'],
                str(AnnotationHeader.TEXT): annotation.text
            }
            record.update({f'{TEXT_TAGS_NAME}.{k}': v for k, v in annotation.text_tags.items()})
            record.update({f'{BODY_TAGS_NAME}.{k}': v for k, v in annotation.body_tags.items()})
            # record.update({f'{TEXT_TAGS_NAME}.{k}': v[0] if len(v) == 1 else v for k, v in annotation.text_tags.items()})
            # record.update({f'{BODY_TAGS_NAME}.{k}': v[0] if len(v) == 1 else v for k, v in annotation.body_tags.items()})
            record_list.append(record)
            column_set.update(record.keys())

        df = pd.DataFrame.from_records(record_list)
        df = df.pipe(_fix_df_nan)
        return df

    # NOTE: there seems to be an issue with how annotations were attached to each source; there can be multiple
    # sources for a single PMID, where only a single INFORMATION annotation exists. This means that for some PMIDs,
    # there will be some CASE or COHORT annotations without associated INFORMATION annotations (and therefore PMIDs).
    # The best solution so far is to use both document titles and source to merge across annotations
    @staticmethod
    def merge_across_document_title_and_source(annotations: List[AugmentedAnnotation]):
        '''
        Take a list of all annotations and copy tags across annotations related through their document title OR source
        @param annotations:
        @return: List[AugmentedAnnotation]
        '''

        source_title_dict = {}
        for a in annotations:
            # source always exists, no need to check
            a_source = a.target[0]['source']
            a_source_deque = source_title_dict.get(a_source, deque())
            if a.type == AnnotationType.INFORMATION.name:
                a_source_deque.appendleft(a)
            else:
                a_source_deque.append(a)
            source_title_dict[a_source] = a_source_deque

            # title might not exist, needs check
            if a.document:
                a_title = a.document['title'][0]
                a_title_deque = source_title_dict.get(a_title, deque())
                if a.type == AnnotationType.INFORMATION.name:
                    a_title_deque.appendleft(a)
                else:
                    a_title_deque.append(a)
                source_title_dict[a_title] = a_title_deque

        # for each unique title, find the information annotation and copy its pmid to other annotations
        for source, a_deque in source_title_dict.items():
            a_info = a_deque.popleft()
            # all INFORMATION annotations will have been appended to the left if there was one, but some sources
            # might not have INFORMATION so a check is still needed
            if a_info.type == AnnotationType.INFORMATION.name:
                while len(a_deque) > 0:
                    annotation = a_deque.popleft()
                    annotation.body_tags.update(a_info.body_tags)
                    annotation.text_tags.update(a_info.text_tags)
