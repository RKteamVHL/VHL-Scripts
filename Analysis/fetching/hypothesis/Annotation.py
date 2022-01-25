from dataclasses import dataclass, field, asdict
from typing import Dict, List
import re
import copy
import enum


class AnnotationType(enum.IntEnum):
    INVALID = -1
    EVIDENCE = enum.auto()
    INFORMATION = enum.auto()
    METHODOLGY = enum.auto()
    CASE = enum.auto()


BODY_TAG_REGEX = re.compile("^(?P<name>\w+):(?P<body>.*)$", flags=re.MULTILINE)

# for tag lists, only mandatory tags are used to determine the annotation type

EVIDENCE_TAGS = [
    "EvidenceStatement",
]

INFORMATION_TAGS = [
    "PMID",
    "Gene"
]

METHODOLOGY_TAGS = [
    "ArticleReferenceSequence",
    "GenotypingMethod",
    "SamplingMethod"
]

CASE_TAGS = [
    "PatientID",
    "KindredID",
    "Case",
    "DiseaseAssertion",
    "FamilyInfo",
    "CasePresentingHPOs",
    "CaseHPOFreeText",
    "CaseNotHPOs",
    "CaseNotHPOFreeText",
    "CasePreviousTesting"
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
            self.tag_dictionary[tag_split[0]] = ",".join(tag_split[1:]).strip() if len(tag_split) > 1 else ""
            pass

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

    def as_dict(self):
        return asdict(self)

    @staticmethod
    def from_dict(d):
        new_annotation = AugmentedAnnotation(**copy.deepcopy(d))
        new_annotation.get_tags_from_body()
        new_annotation.get_tags_from_tags_list()
        new_annotation.assign_type()
        return new_annotation
