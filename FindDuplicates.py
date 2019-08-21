#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Install a pip package in the current Jupyter kernel
import sys

import os
import json
import functools

from pprint import pprint
from RISparser import readris

ris_files = os.listdir("files")


ris_without_doi = []
#ALL RIS entries, with the doi serial being the unique key
ris_entries = {}
ris_list = []

for ris_file in ris_files:
    with open(os.path.join("files", ris_file), 'r', encoding="utf8") as bibliography_file:
        new_list = list(readris(bibliography_file))
        ris_list.extend(new_list) 


        for entry in new_list:
            if "doi" not in entry:
                ris_without_doi.append(entry)
            else:
                doi_split = entry["doi"].split("/")
                serial= "/".join([doi_split[-2], doi_split[-1]]) 

                if serial not in ris_entries:
                    ris_entries[serial] = []
                ris_entries[serial].append(entry)
            
print("Out of {0} entries in the ris file, there are {1} without dois".format(len(ris_list), len(ris_without_doi)))




#seeing how many unique dois there are in each set of entries
ris_total = 0
for value in ris_entries.values():
    ris_total += len(value)
print("Out of the {0} ris entries with dois, there are {1} unique dois".format(ris_total, len(ris_entries.keys())))




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



