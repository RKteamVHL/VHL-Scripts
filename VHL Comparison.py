#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Install a pip package in the current Jupyter kernel
import sys
get_ipython().system('{sys.executable} -m pip install RISparser')

import os
import json
import functools
from pprint import pprint
from RISparser import readris

#the ris file given to Marta for the new imports in May
ris_file = 'mendeley_vhl.ris'

#the json file here contains ALL entries that were entered on May 20th
json_file = 'mendeley_vhl.json'


out_file = "mendeley_unimported.json"


# In[2]:


ris_without_doi = []
#ALL RIS entries, with the doi serial being the unique key
ris_entries = {}
ris_list = None
with open(ris_file, 'r') as bibliography_file:
    ris_list = list(readris(bibliography_file))


    for entry in ris_list:
        if "doi" not in entry:
            ris_without_doi.append(entry)
        else:
            doi_split = entry["doi"].split("/")
            serial= "/".join([doi_split[-2], doi_split[-1]]) 

            if serial not in ris_entries:
                ris_entries[serial] = []
            ris_entries[serial].append(entry)
            
print("Out of {0} entries in the ris file, there are {1} without dois".format(len(ris_list), len(ris_without_doi)))


# In[3]:


json_without_doi = []
#ALL json entries, with the doi serial being the unique key
json_entries = {}
json_list = None
with open(json_file, 'r') as bibliography_file:
    json_list = list(json.load(bibliography_file))
    

    
    for entry in json_list:
        if "doi" not in entry["identifiers"]:
            json_without_doi.append(entry)
        else:
            doi_split = entry["identifiers"]["doi"].split("/")
            serial= "/".join([doi_split[-2], doi_split[-1]]) 
            
            if serial not in json_entries:
                json_entries[serial] = []
            json_entries[serial].append(entry)
            
print("Out of {0} entries in the json file, there are {1} without dois".format(len(json_list), len(json_without_doi)))


# In[4]:


#seeing how many unique dois there are in each set of entries
ris_total = 0
for value in ris_entries.values():
    ris_total += len(value)
print("Out of the {0} ris entries with dois, there are {1} unique dois".format(ris_total, len(ris_entries.keys())))

json_total = 0
for value in json_entries.values():
    json_total += len(value)
print("Out of the {0} json entries with dois, there are {1} unique dois".format(json_total, len(json_entries.keys())))


# In[5]:


#theres a discrepancy here with how many unique dois there are...
#here are all of the non-unique dois


# In[6]:


#for the mendeley ris file
ris_repeat_dois = set()
for doi in ris_entries:
    if(len(ris_entries[doi])>1):
        ris_repeat_dois.add(doi)
        
print("There are {0} dois with repeats in the ris file:\n {1}".format(len(ris_repeat_dois), ris_repeat_dois))


# In[7]:


#for the mendeley json file
json_repeat_dois = set()
for doi in json_entries:
    if(len(json_entries[doi])>1):
        json_repeat_dois.add(doi)
      
print("There are {0} dois with repeats in the json file:\n {1}".format(len(json_repeat_dois), json_repeat_dois))


# In[8]:


#and the non-unique dois in both files:
common_repeat_dois = ris_repeat_dois.intersection(json_repeat_dois)
print("There are {0} repeated dois common in both files:\n {1}".format(len(common_repeat_dois), common_repeat_dois))


# In[9]:


#it turns out all of the dois in the json file are in the RIS file


# In[10]:


#now, finding which entries did not get added to the json file on May 20th
json_set = set(json_entries.keys())
ris_set = set(ris_entries.keys())

#finding which serials are not inside the json (May 20th) set
json_extras = json_set.difference(ris_set)
ris_extras = ris_set.difference(json_set)
# print(json_extras) #all of the json entries were in the ris file
print("There are {0} dois in the ris file that don't appear in the json file: {1}\n".format(len(ris_extras),  ris_extras))


# In[11]:


out_set = {key: ris_entries[key] for key in ris_extras}


# In[12]:


with open(out_file, 'w') as bibliography_file:
    json.dump(out_set, bibliography_file)


# In[13]:


### manual inspection results ###

## The following are some examples of work duplication
# 10.18632/oncotarget.23470 (eg. Reviewed - Andrea / Irrelevant, Reviewed - Clarissa / Irrelevant)
# 10.1297/cpe.27.87
# 10.1021/acs.jmedchem.7b00635
# 10.1038/s41598-018-21524-5
# 10.1016/j.semnephrol.2018.01.006
# 10.1016/j.clgc.2018.01.013
# 10.1172/jci.insight.92193
# 10.1021/acs.jmedchem.6b01781
# 10.1177/1078155218790342
# 10.3892/ol.2018.9328
# 10.2217/pgs-2017-0160 (Andrea Irrelevant, Kelly No access )
# 10.1038/s41598-017-11035-0 (Yasser, Clariss Irrelevant)
# 10.1158/1535-7163.MCT-17-1076 (Pre 2016 Payal, Clarissa Irrelevant)
# 10.1016/j.bbadis.2018.07.016 (Payal Irrelevant, Veronica Unsure)
# 10.1136/jmedgenet-2018-105567 (Veronica - Not Relevant)

## The following are some examples where there are duplicate entries
# 10.1021/acs.jmedchem.7b00635
# 10.1080/07347332.2018.1450320
# 10.1038/s41598-018-21524-5
# 10.1016/j.semnephrol.2018.01.006
# 10.2967/jnumed.118.216796
# 10.1038/s41598-018-27220-8
# 10.1016/j.humpath.2018.07.033
# 10.2217/pgs-2017-0160 (Andrea Irrelevant, Kelly No access )
# 10.1194/jlr.M083311
# 10.1158/1535-7163.MCT-17-1076 (Pre 2016 Payal, Clarissa Irrelevant)

## The following, from the ris file extras, do not appear in the private vhl repo
# 10.1186/s12935-018-0530-2 

