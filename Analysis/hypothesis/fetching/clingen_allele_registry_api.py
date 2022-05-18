import urllib.request
import urllib.parse
import urllib.error

import json

BASE_URL = "https://reg.genome.network/allele/"


def get_variant_by_caid(caid: str):
    url = f"{BASE_URL}{caid}"
    request = urllib.request.Request(url, method="GET")
    response_dict = None
    try:
        with urllib.request.urlopen(request) as response:
            response_dict = json.loads(response.read().decode('utf-8'))
    except urllib.error.HTTPError:
        pass

    return response_dict
