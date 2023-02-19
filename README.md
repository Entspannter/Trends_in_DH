# Analyzing Trends in the Use of Digital Health Technologies in Clinical Research and Beyond – Evidence from ClinicalTrials.gov and Other Sources
Files and code for the first two parts of the master thesis "Analyzing Trends in the Use of Digital Health Technologies in Clinical Research and Beyond – Evidence from ClinicalTrials.gov and Other Sources" submitted by Lars Masanneck at Hasso-Plattner Institute in February 2023. 


Parts of the third part of the thesis have been published. Read the [manuscript](https://www.nature.com/articles/s41746-023-00767-1#citeas), which was published in Nature npj Digital Medicine on February 10th.
Code and files belonging to this part are available in the respective [repository](https://github.com/Entspannter/DHTs-in-neurology-trials).


## Code structure
All code in is the folder 'Basic_API_scripts'. 
An essential part of this thesis has been to analyze trends related to digital health. For this, ClinicalTrials.gov, PubMed and Twitter have been queried. The overarching class of the written scripts is the Multiquery class in the multi_query.py. It can be used for gathering any keyword or list of keywords on ClinicalTrials and Pubmed. As Twitter has revoked research access to its API, querying it has gotten much harder. In case the user still has access, querying single queries (1024 characters) remains possible. The structure for different analyses throughout this work can be seen in the accompanying Jupyter notebooks (.ipynb files). In order to recreate any of these, please fill in your API tokens for Pubmed and, if available, Twitter in the api_core_functions.py file (see comments in there).

Complete queries using the Multiquery package are stored in the 'Graph' folder. The 'CTGOV_DTX' folder contains all graphs of the second part of this work, with the code for these graphs available in the Jupyter notebook 'jup_analyze_CTGOV_DTX.ipynb'.


## Data
As not all datasets could be pushed to this repository, some .csv or .xlsx files might be missing in the 'Data' folder. As all data has been pulled from publicly available sources, one should be able to recreate the data using the 'Methods' section of the thesis. In case someone is interested in the exact datasets, feel free to reach out to the author on GitHub.

## Why are there packages and graphs mentioned in the code but not included in the thesis?
While most files were thorougly cleaned up before publishing, I intentionally left some additional analyses (e.g., timeseries analyses) in the Jupyter notebooks as show different strategies and foci that were not part of the final thesis. These chunks of code also explain many of the packages included in the requirements that might be surprising when just reading the thesis.



