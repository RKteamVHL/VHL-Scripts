import os
import urllib.request
import urllib.parse

import json

from .Annotation import AugmentedAnnotation

INPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../files", "input")
TOKEN_FILE = os.path.join(INPUT_DIR, "secret_token.txt")
SECRET_TOKEN = ""

if not os.path.isdir(INPUT_DIR):
    os.makedirs(INPUT_DIR)

with open(TOKEN_FILE, "r", encoding="utf-8") as file:
    SECRET_TOKEN = file.readline().strip()


BASE_HEADER = {
    "Host": "hypothes.is",
    "Accept": "application/json",
    "Authorization": f"Bearer {SECRET_TOKEN}"
}

ANNOTATION_LIMIT = 200


def get_annotations_by_group(group_id, search_after):
    annotations = []
    reached_end = False
    next_checkpoint = search_after
    total_annotations = 0
    while not reached_end:
        h_data = {
            "limit": ANNOTATION_LIMIT,
            "sort": "created",
            "order": "asc",
            "search_after": next_checkpoint,
            "group": group_id
        }
        h_header = dict(BASE_HEADER)
        h_url = f"https://api.hypothes.is/api/search?{urllib.parse.urlencode(h_data)}"
        h_request = urllib.request.Request(h_url, headers=h_header, method="GET")
        with urllib.request.urlopen(h_request) as response:
            response_dict = json.loads(response.read().decode('utf-8'))

            for row in response_dict['rows']:
                new_annotation = AugmentedAnnotation.from_dict(row)
                annotations.append(new_annotation)

            total_annotations = response_dict['total']
            if len(annotations) == total_annotations:
                reached_end = True
            next_checkpoint = response_dict['rows'][-1]['created']

    return annotations
