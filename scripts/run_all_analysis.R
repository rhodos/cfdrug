setwd('/Users/rhodos/Desktop/Projects/cf/cfdrug/cfdrug')

library(magrittr)
library(dplyr)

source('src/Utils.R')
source('src/DrugRank.R')
source('src/CMAP.R')

# compound signature ranking
#source('scripts/run_all_drug_scores.R')

# final compound selection
source('scripts/get_top_drugs.R')

# in silico validation of drug ranking
source('scripts/in_silico_validation_of_drug_ranking.R')
source('scripts/in_silico_validation_dotplot.R')

# post hoc mechanistic analysis
source('scripts/post_hoc_mechanistic_analysis.R')