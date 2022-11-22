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

# BODY, TEXT, and COMPUTED need to be separated so that we can discern where the tags came from
BODY_TAGS_NAME = "BODY"
TEXT_TAGS_NAME = "TEXT"
COMPUTED_TAGS_NAME = "COMPUTED"


# 'header' is used here to indicate that these are the string headers for the csv tabular format of the annotations
# these were made an enum (as compared to a dict or list) for re-usability and type-hinting in IDEs
class AnnotationHeader(enum.Enum):
    CLINVAR = (BODY_TAGS_NAME, "ClinVarID")
    CAID = (BODY_TAGS_NAME, "CAID")
    CIVICNAME = (BODY_TAGS_NAME, "CivicName")
    EXPERIMENTALASSAY = (BODY_TAGS_NAME, "ExperimentalAssay")
    UNREGISTEREDVARIANT = (BODY_TAGS_NAME, "UnregisteredVariant")
    FAMILYPEDIGREE = (BODY_TAGS_NAME, "FamilyPedigree")
    MUTATIONTYPE = (BODY_TAGS_NAME, "MutationType")
    AMINOACIDCHANGE = (BODY_TAGS_NAME, "AminoAcidChange")
    PROTEINPOSITION = (BODY_TAGS_NAME, "ProteinPosition")
    DISEASEENTITY = (BODY_TAGS_NAME, "DiseaseEntity")
    VARIANT = (TEXT_TAGS_NAME, "Variant")
    ARTICLEREFERENCESEQUENCE = (TEXT_TAGS_NAME, "ArticleReferenceSequence")
    PREVIOUSLYPUBLISHED = (TEXT_TAGS_NAME, "PreviouslyPublished")
    PMID = (TEXT_TAGS_NAME, "PMID")

    TYPE = ('type',)

    def __str__(self):
        if len(self.value) == 1:
            return self.value[0]
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
    "EvidenceStatement",
]

INFORMATION_TAGS = [
    "PMID",
    "Gene",
    "StandardizedReferenceSequence"
]

METHODOLOGY_TAGS = [
    "ArticleReferenceSequence",
    "GenotypingMethod",
    "SamplingMethod"
]

CASE_TAGS = [
    "CasePresentingHPOs"
]

COHORT_TAGS = [
    "GroupPresentingHPOs"
]

ASSAY_TAGS = [
    "ExperimentalAssay"
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
        # checking for evidence statement annotation
        elif any([key in self.body_tags for key in EVIDENCE_TAGS]):
            self.type = AnnotationType.EVIDENCE.name

        # checking for article info annotation
        elif any([key in self.text_tags for key in INFORMATION_TAGS]):
            self.type = AnnotationType.INFORMATION.name

        # checking for methodology annotation
        elif any([key in self.text_tags for key in METHODOLOGY_TAGS]):
            self.type = AnnotationType.METHODOLOGY.name

        # checking for case annotation
        elif any([key in self.text_tags for key in CASE_TAGS]):
            self.type = AnnotationType.CASE.name

        # checking for COHORT annotation
        elif any([key in self.text_tags for key in COHORT_TAGS]):
            self.type = AnnotationType.COHORT.name

        # checking for experiment assay annotation
        elif any([key in self.body_tags for key in ASSAY_TAGS]):
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

    # this function returns a dataframe with a couple of caveats:
    # 1.    if the tag came from the body, it will be prepended with the BODY_TAGS_NAME constant, otherwise if it comes from
    #       the text, it will be prepended with TEXT_TAGS_NAME
    # 2.    individual cells will contain lists, which won't get formatted properly if exported to a csv file
    # 3.    some data (i.e., user info) is excluded from the output csv file
    @staticmethod
    def df_from_annotations(annotations: List[AugmentedAnnotation]):
        record_list = []
        column_set = set()
        for annotation in annotations:
            record = {
                "type": annotation.type,
                'uri': annotation.links['html'],
                "source": annotation.target[0]['source'],
                "text": annotation.text
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
