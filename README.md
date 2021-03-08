# cfdrug
Code for manuscript in CPT: Pharmacometrics & Systems Pharmacology, "Integrative chemogenomic analysis identifies small molecules that partially rescue ΔF508-CFTR for cystic fibrosis"

The code covers the following key aspects of the manuscript:
- compound scoring, ranking, and final selection for testing
- in silico validation of scoring strategy using known CF corrector compounds
- post hoc transcriptional analysis to infer potential mechanisms 

To get started, download data from the following location: https://drive.google.com/file/d/1QqSKSmbhC1j291pd-OENByAOrYccxd5a/view?usp=sharing

Then update the config file to point to the corresponding code, input and output directories.

Then (after installing necessary R packages) you should be able to run the scripts. The main entrypoint is scripts/run_all_analysis.R.
