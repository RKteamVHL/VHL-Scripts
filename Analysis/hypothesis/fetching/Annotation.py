from dataclasses import dataclass, field, asdict
from typing import Dict, List
import numpy as np
import pandas as pd
import re
import copy
import enum

NULL_TERMS = ["unknown", "none", "", "N/A"]

BODY_TAGS_NAME = "BODY"
TEXT_TAGS_NAME = "TEXT"


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

    def get_tags_from_text(self):
        tag_match = re.finditer(BODY_TAG_REGEX, self.text)
        if tag_match:
            # TODO: some assertion/error checking on tag_dict
            for match in tag_match:
                tag_dict = match.groupdict()
                tag_name = tag_dict["name"]
                tag_value = tag_dict["body"].strip()

                # if the tag name is already in the tag_dict, fetch its corresponding list; else, return an empty list
                new_list = self.text_tags.get(tag_name, [])

                # if the value is already in the tag_dict, ignore it (removes tags that get duplicated in body
                # and tag list)
                if tag_value not in new_list:
                    new_list.append(tag_value)
                self.text_tags[tag_name] = new_list
        else:  # error or log something here
            pass

    def get_tags_from_tags_list(self):
        for tag in self.tags:
            tag_split = tag.split(":")
            tag_name = ""
            tag_value = ""

            # double tags, in the format TagName1:TagName2:TagValue i.e., AgeOfPresentation
            if len(tag_split) == 3:
                tag_name = tag_split[0].strip()
                tag_name2 = tag_split[1].strip()
                tag_value = {tag_name2: tag_split[2].strip()}

            # normal tags in the format TagName:TagValue
            elif len(tag_split) == 2:
                tag_name = tag_split[0].strip()
                tag_value = tag_split[1].strip()

            # flag tags, in the format TagName
            elif len(tag_split) == 1:
                tag_name = tag_split[0].strip()
                tag_value = True

            # if the tag name is already in the tag_dict, fetch its corresponding list; else, return an empty list
            new_list = self.body_tags.get(tag_name, [])

            # if the value is already in the tag_dict, ignore it (removes tags that get duplicated in body
            # and tag list)
            if tag_value not in new_list:
                new_list.append(tag_value)
            self.body_tags[tag_name] = new_list

    def assign_type(self):
        # checking for evidence statement annotation
        if any([key in self.body_tags for key in EVIDENCE_TAGS]):
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

        elif len(self.references) > 0:
            self.type = AnnotationType.REPLY.name

    def as_dict(self):
        return asdict(self)

    @staticmethod
    def from_dict(d):
        new_annotation = AugmentedAnnotation(**copy.deepcopy(d))
        new_annotation.get_tags_from_text()
        new_annotation.get_tags_from_tags_list()
        new_annotation.assign_type()
        return new_annotation


def _fix_df_nan(df: pd.DataFrame):
    out_df = df
    for term in NULL_TERMS:
        out_df = out_df.replace(term, np.NaN)
    return out_df

# this function returns a dataframe with a couple of caveats:
# 1.    only the tag type and properties in the annotation tag_dictionary are included
# 2.    if the tag is one which can have multiple values simultaneously ie. DiseaseEntity or AgeOfPresentation,
#       then only the first of the values is returned
def get_df_from_annotations(annotations: List[AugmentedAnnotation]):
    record_list = []
    column_set = set()
    for annotation in annotations:
        record = {"type": annotation.type}
        record.update({f'{TEXT_TAGS_NAME}.{k}': v for k, v in annotation.text_tags.items()})
        record.update({f'{BODY_TAGS_NAME}.{k}': v for k, v in annotation.body_tags.items()})
        record_list.append(record)
        column_set.update(record.keys())

    df = pd.DataFrame.from_records(record_list)
    df = df.pipe(_fix_df_nan)
    return df
