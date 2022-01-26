from dataclasses import dataclass, field, asdict
from typing import Dict, List
import pandas as pd
import re
import copy
import enum


class AnnotationType(enum.IntEnum):
    INVALID = -1
    REPLY = enum.auto()
    EVIDENCE = enum.auto()
    INFORMATION = enum.auto()
    METHODOLGY = enum.auto()
    CASE = enum.auto()
    COHORT = enum.auto()


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
    # "Case",
    # "CasePresentingHPOs",
    # "CaseHPOFreeText",
    # "CaseNotHPOs",
    # "CaseNotHPOFreeText",
    # "CasePreviousTesting"
]

COHORT_TAGS = [
    "GroupPresentingHPOs"
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
    body_tags: List[str] = field(default_factory=list)
    tag_dictionary: Dict[str, str] = field(default_factory=dict)
    annotation_type: str = AnnotationType.INVALID.name

    def get_tags_from_body(self):
        tag_match = re.finditer(BODY_TAG_REGEX, self.text)
        if tag_match:
            # TODO: some assertion/error checking on tag_dict
            for match in tag_match:
                tag_dict = match.groupdict()
                self.tag_dictionary[tag_dict["name"]] = tag_dict["body"].strip()
        else:  # error or log something here
            pass

    def get_tags_from_tags_list(self):
        for tag in self.tags:
            tag_split = tag.split(":")
            tag_name = ""
            tag_value = ""

            # normal tags in the format TagName:TagValue
            if len(tag_split) > 1:
                tag_name = ":".join(tag_split[0:-1]).strip()
                tag_value = tag_split[-1].strip()

            # flag tags, in the format TagName
            elif len(tag_split) == 1:
                tag_name = tag_split[0].strip()

            self.tag_dictionary[tag_name] = tag_value

    def assign_type(self):
        # checking for evidence statement annotation
        if any([key in self.tag_dictionary for key in EVIDENCE_TAGS]):
            self.annotation_type = AnnotationType.EVIDENCE.name

        # checking for article info annotation
        elif any([key in self.tag_dictionary for key in INFORMATION_TAGS]):
            self.annotation_type = AnnotationType.INFORMATION.name

        # checking for methodology annotation
        elif any([key in self.tag_dictionary for key in METHODOLOGY_TAGS]):
            self.annotation_type = AnnotationType.METHODOLGY.name

        # checking for case annotation
        elif any([key in self.tag_dictionary for key in CASE_TAGS]):
            self.annotation_type = AnnotationType.CASE.name

        # checking for COHORT annotation
        elif any([key in self.tag_dictionary for key in COHORT_TAGS]):
            self.annotation_type = AnnotationType.COHORT.name

        elif len(self.references) > 0:
            self.annotation_type = AnnotationType.REPLY.name

    def as_dict(self):
        return asdict(self)

    def as_dataframe(self):
        """Returns the annotation as a pandas Dataframe

        This strips the annotation of all it's metadata
        """
        pass


    @staticmethod
    def from_dict(d):
        new_annotation = AugmentedAnnotation(**copy.deepcopy(d))
        new_annotation.get_tags_from_body()
        new_annotation.get_tags_from_tags_list()
        new_annotation.assign_type()
        return new_annotation


"""Notes:
GroupID/KindredID and GroupID/ KindredID -> both appear in body tags
INVALID annotations -> 
    - evidencestatement instead of EvidenceStatement for some annotations
    - most are caused by GroupID/KidredID inconsistencies
    - some are replies to other annotations
    - "Evidence": "", in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5503545/#annotations:piBWZEfzEeyuGbdGPpU0qw
    -  GroupID: 2, where the spaces mess up the regex. https://hyp.is/P5YEylOrEeyZUH-REJFQqQ/onlinelibrary.wiley.com/doi/abs/10.1002/humu.20385
    -  GroupID: 20, same as above. https://hyp.is/Qsb8-FftEeyKeXeNmFIZ-A/onlinelibrary.wiley.com/doi/abs/10.1002/humu.20385
    -  GroupID: 36, same as above. https://hyp.is/4YVg3l1ZEey-rK9ETACNew/onlinelibrary.wiley.com/doi/abs/10.1002/humu.20385
    - the following seems to be EvidenceStatement, but is missing the tag. https://hyp.is/2hjYGmgAEeyxyyeAzATdxg/iovs.arvojournals.org/article.aspx?articleid=2162471
    - https://hyp.is/fvvwcmgSEey8rcsikkfNuA/jmg.bmj.com/content/39/7/e38
    - https://hyp.is/HophImj6Eey5L0efg0JjHQ/jmg.bmj.com/content/39/7/e38
    - "" https://hyp.is/qd0YnHrnEey7X5esN467Aw/onlinelibrary.wiley.com/doi/abs/10.1002/humu.21219
    
ClinVarID:Yes vs ClinVarID:XXX -> having both makes the second get overwritten

Last known age vs age of presentation -> do we care about last known?

The following is an example where age of presentation is 'mid 40s', but no exact date is given. Does this pop up a lot, and if so where is the age range usually put?

"""