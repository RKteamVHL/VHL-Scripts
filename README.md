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

## Running all Scripts
The core functions are implemented in the Analysis module, which can be run with the following command:
```commandline
(test_env) C:\\...\\VHL-Scripts>py -m Analysis --createfigs
```
Running the command without the "createfigs" argument will run all of the data analysis and statistical tests without creating
all of the figures, which can take some time.

Additionally, running the scripts with the following command:
```commandline
(test_env) C:\\...\\VHL-Scripts>py -m Analysis --createfigs --validation
```
Will also run the scripts that validate our data against the [UMD VHL](http://www.umd.be/VHL/) database. 
## Figures
All relevant output figures are saved in the "statistics" folder that gets created after running the scripts. The folders
are separated by patient, kindred, and variant. Inside these separate folders are the "data" and "figures" subfolders; the "figures"
folder contains the matplotlib-generated .pdf figures, and the "data" folder contains the raw data used to generate the figures.

Of interest are:

##### Fig 1. Age-related penetrance for patients that present with a single phenotype
```commandline
File: statistics\patient\figures\penetrance.pdf
Code: Analysis\features\kimstudents_dataframe_views.py (line 382)
```
##### Fig. 2. Phenotype co-occurrence ratios for (A) patient-, (B) family- and (C) variant-based data
```commandline
File: statistics\(patient, family, or variant)\figures\phenotype_correlation_ratio.pdf
Code: Analysis\features\kimstudents_dataframe_views.py (line 438)
```
##### Fig 3. Frequency of missense variants along the VHL gene for (A) patient-, (B) family- and (C) variant-based data. 
```commandline
Files: statistics\(patient, family, or variant)\figures\codon_histogram.pdf
Code: Analysis\features\kimstudents_dataframe_views.py (line 282)
```
##### Fig 4. Distribution of truncating and non-truncating variant by phenotype for (A) patient-, (B) family- and (C) variant-based data
```commandline
Files: statistics\(patient, family, or variant)\figures\grouped_mutant_type_counts.pdf
Code: Analysis\features\kimstudents_dataframe_views.py (line 240)
```
##### Fig 5. Frequency of coding variants in protein and functional domains for (A) patient-, (B) family- and (C) variant-based data
```commandline
Files: statistics\(patient, family, or variant)\figures\regions.pdf
Code: Analysis\features\kimstudents_dataframe_views.py (line 96)
```
##### Fig 6. Cluster phenotype, variant type, variant domain and codon distribution for two patient clusters.
```commandline
Files:  statistics\patient\cluster\cluster_labels_best\clustered_generalized_phenotype_counts.pdf
        statistics\patient\cluster\cluster_labels_best\clustered_domain_counts.pdf
        statistics\patient\cluster\cluster_labels_best\clustered_grouped_mutation_type_counts.pdf
        statistics\patient\cluster\cluster_labels_best\clustered_codon_start.pdf

Code: Analysis\features\kimstudents_dataframe_views.py (line 632)
``` 
##### Fig 7. Cluster phenotype (A), variant type (B), variant domain (C) and codon statistics (D) for each of the 4 patient clusters
```commandline
Files:  statistics\patient\cluster\cluster_labels_second\clustered_generalized_phenotype_counts.pdf
        statistics\patient\cluster\cluster_labels_second\clustered_domain_counts.pdf
        statistics\patient\cluster\cluster_labels_second\clustered_grouped_mutation_type_counts.pdf
        statistics\patient\cluster\cluster_labels_second\clustered_codon_start.pdf

Code: Analysis\features\kimstudents_dataframe_views.py (line 632)
``` 
## Files
The scripts output the fetched individual input files to the /output directory, which get combined into a single masterlist.
The first step of filtering is to drop all entries that do not have either a phenotype or a mutation- the output of this step is
saved in:
```
statistics/summary/postdropsupplementary_1.csv
``` 
The next step is to create additional annotated columns with clean values to 
be used in further analysis. The resulting file is:
```
statistics/summary/filtered_out.csv
``` 
#### Table 1
The summary statistics for the total number of patients, kindreds, and variants are saved to:
```
statistics/summary/postdropsummary.csv
``` 
The breakdown of these summaries by phenotype and mutation type are saved in :
```
statistics/summary/summary_by_type.csv
```
Together, these two files are used to create Table 1.