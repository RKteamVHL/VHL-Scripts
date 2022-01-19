from dataclasses import dataclass
from typing import Set, Type, Dict, List
import re
import copy

BODY_TAG_REGEX = re.compile("^(?P<name>\w+):(?P<body>.*)$")


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
    tags: Set[str]
    group: str
    permissions: Dict[str, List[str]]
    target: List[Dict[str, List]]
    document: Dict[str, List[str]]
    links: Dict[str, str]
    flagged: bool
    hidden: bool
    references: List[str]
    user_info: Dict[str, str]

    @staticmethod
    def from_dict(d):
        if "references" not in d:
            d["references"] = []
        new_annotation = HypothesisAnnotation(**copy.deepcopy(d))

        return new_annotation

@dataclass
class AugmentedAnnotation(HypothesisAnnotation):
    """Augmented Hypothesis Annotation class

    This class contains additional properties and methods for parsing and storing tags embedded into the annotation
    body, along with converting all tags to a dictionary (rather than a string)
    """
    body_tags: Set[str]
    tag_dictionary: Dict[str, str]

    def create_tag_dict(self):
        self.assign_tags_to_dict()
        self.get_body_tags()

    def get_body_tags(self):
        tag_search = re.search(BODY_TAG_REGEX, self.text)
        if tag_search:
            tag_dict = tag_search.groupdict()
            #TODO: some assertion/error checking on tag_dict
            self.tag_dictionary[tag_dict["name"]] = tag_dict["body"]
        else:  # error or log something here
            pass

    def assign_tags_to_dict(self):
        for tag in self.tags:
            tag_split = tag.split(":")
            self.tag_dictionary[tag_split[0]] = tag_split[1:] if len(tag_split) > 1 else ""








