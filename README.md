# VHL-Scripts
## Requirements
Python version 3.7 or later is required for running all scripts. The libraries used by the scripts include:
```requirements.txt
numpy~=1.21.2           # numpy, pandas, and scipy are used for all data manipulation, including  
pandas~=1.3.3           # acquistion, cleaning/filtering, and metric calculations
scipy~=1.7.1 
obonet~=0.3.0           # obonet is used for Human phenotype and sequennce ontologies
biopython~=1.79         # biopython is used for pVHL protein changes
scikit-learn~=0.24.2    # scikit-learn is used with snfpy and networkx for the spectral clustering of patient, kindred, and variant graphs
networkx~=2.6.3
snfpy~=0.2.2
matplotlib~=3.4.3       # matplotlib and cycler are used in the creation of all figures
cycler~=0.10.0
seaborn~=0.11.2         # seaborn is used for the heatmap visualizations
certifi~=2021.5.30
requests~=2.26.0        # requests is used for the fetching of data from remote online sources
lxml~=4.6.3             # lxml is used in the html-scrubbing portion of cross-validation with the UMD database
```
## Installation
### Virtual Environment
To run the scripts, first all requirements must be installed. We recommend first setting up a Python virtual environment.
[The documentation for venv](https://docs.python.org/3/library/venv.html) fully details how to set one up, but for Windows 10
the commands might look like:
```commandline
C:\\...\\VHL-Scripts> py -m venv env
C:\\...\\VHL-Scripts> .\test_env\Scripts\activate
```

These should be entered from a cmd or powershell window inside VHL-Scripts directory. (pressing shift+rightclick in the
root directory, then clicking "Open PowerShell window here" while in the root of the repository opens a powershell 
window already scoped to the right location)

After activating the virtual environment, the name of the environment will prefix the command line:
```commandline
(test_env) C:\\...\\VHL-Scripts>
```
### Libraries
With the virtual environment activated, all depenencies can be installed by running:
```commandline
(test_env) C:\\...\\VHL-Scripts>pip install -r requirements.txt
```
If there are no errors, the environment should now be set up to run all scripts properly

##Hypothes.is
The VHL Hypothesis Annotation group on hypothes.is contains 5000+ expert-curated annotations for ~370 papers on germline
VHL variants. The hypothesis python package contains script related to downloading, cleaning, and summarizing statistics
 of these annotations. Further information on running these scripts can be found in 
[the hypothesis package](hypothesis/README.md). Moving forward, the VHL Hypothesis Annotation group will maintain our 
most up-to-date data.

##Kim Student Masterlist
The Kim Student VHL variant Masterlist contains expert-curated variant information extracted from 427 papers. 
A comprehensive description of the collection, screening, analysis, and discussion of the data can be found our paper:
[Large scale genotype- and phenotype-driven machine learning in Von Hippel-Lindau disease](https://doi.org/10.1002/humu.24392).
More information on running the scripts and understanding their outputs can be found in [the kim masterlist package](kim_masterlist/README.md).