import urllib.request
import urllib.parse
import json

SERVICE_ROOT = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/'

EMAIL = 'sean.delong@ryerson.ca'
TOOLNAME = 'VHL-Scripts'


def validate_pmid(pmids):
	"""
	Looks up input PMID to make sure it exists on PMC
	pmids (str) or List[str]: single pmid or list of pmids
	"""
	id_str = ''
	if isinstance(pmids, list):
		id_str = ','.join(pmids)
	elif isinstance(pmids, str):
		id_str = pmids

	data = {
		'email': EMAIL,
		'tool': TOOLNAME,
		'ids': id_str,
		'format': 'json'
	}
	request = urllib.request.Request(url=f'{SERVICE_ROOT}?{urllib.parse.urlencode(data, safe=",")}', method='GET')

	with urllib.request.urlopen(request) as response:
		response_dict = json.loads(response.read().decode('utf-8'))
		pass
	#TODO: finish this utility function if it ever becomes important

#NOTE: this function only works for articles on PMC- other articles that have PMIDs or DOIs but
#are not on PMC come back as 'invalid article id'
if __name__ == '__main__':
	validate_pmid('28697140')
	validate_pmid(['28697140', '28169069'])
