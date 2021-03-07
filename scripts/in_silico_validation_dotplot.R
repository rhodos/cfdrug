library(magrittr)
library(ggbio)
library(tidyr)
library(forcats)

set.seed(123)

# This script performs "in silico evaluation" of the connectivity scores using 15 corrector compounds that are represented in these databases. This version of the script looks at the individual pathway and signatures scores as well as the overall score, to try to make a case for data integration

out = list()
for(database in c('cmap','lincs')){
  for(metric in c('ks','xsum')){
nPerm=1000000

##### Load scores data frame #####
scores = LoadScores()
if(database == 'cmap'){
  scores = list(ks=scores$cmap_ks, xsum=scores$cmap_xsum)
  cfDrugs = GetKnownCorrectors()
  compoundName = 'compound'
}else if(database == 'lincs'){
  scores = list(ks=scores$lincs_ks, xsum=scores$lincs_xsum)
  compoundName = 'pert_iname'
  scores$ks = CollapseDuplicateRows(scores$ks, compoundName)
  scores$xsum = CollapseDuplicateRows(scores$xsum, compoundName)
  cfDrugs = toupper(GetKnownCorrectors())
}
nTotal = nrow(scores$ks)
S = scores[[metric]]

##### Extract relevant score columns #####
PW = S[,c('normScore.folding.pathway',
            'normScore.nitro.pathway',
            'normScore.degradation.pathway',
            'normScore.er2golgi.pathway',
            'normScore.sorting.pathway',
            'normScore.lysosome.pathway',
            'normScore.membrane.pathway',
            'normScore.clathrin.pathway',
            'normScore.traffic.pathway')]
names(PW) = sapply(names(PW), function(x) gsub(x, pattern='normScore[.]|[.]pathway', replacement='') %>% toupper)

# Extract signature scores
dx_raw = -apply(S[, c('normScore.50.pighuman.sig', 'normScore.100.pighuman.sig','normScore.150.pighuman.sig')], 1, mean)
SG = cbind(dx_raw, S[,c('normScore.pighuman.sig','normScore.cold.sig','normScore.rnai.sig')])
names(SG) = c('DISEASE_RAW', 'DISEASE_ABS', 'COLD_SHOCK','RNAi')

allScores = cbind(
  CellConsensus=S$multi_cell_score,
  CellSpecific=S$top_cell_score,
  SG, abs(PW))
S = cbind(compound=S[[compoundName]], allScores)

ES = sapply(names(allScores), function(score){
  print(score)
  ComputeCFEnrichmentScore(S[,c('compound', score)], scoreName=score, cfDrugs=cfDrugs, compoundName='compound')$out$ES
}) %>% setNames(names(allScores))

### Compute empirical p-values of cell specific enrichment scores:
nTotal = nrow(S)
nIntersect = length(intersect(cfDrugs, S[,'compound']))
x = sapply(1:nPerm, function(i) GSEA.EnrichmentScore2(1:nTotal, sample(nTotal,nIntersect,replace=FALSE)))
p = sapply(ES, function(es) length(which(x >= es)) / nPerm)

### Compute empirical p-value of cell consensus score
### (computed separately since there is a different number of drugs from which to choose)
multiCellDrugs = S[!is.na(S$CellConsensus),'compound']
nDrugsMultiCell = length(multiCellDrugs)
nCF = length(intersect(multiCellDrugs, cfDrugs))
x = sapply(1:nPerm, function(i) GSEA.EnrichmentScore2(1:nDrugsMultiCell, sample(nDrugsMultiCell, nCF, replace=FALSE)))
p[1] = length(which(x >= ES[[1]])) / nPerm

out[[metric]][[database]] = cbind(ES,p)
}}

# Make a dotplot
M = melt(out) %>% spread(X2, value) %>% setNames(c('score_column', 'database','metric','ES','pval'))
M$score_column %<>% factor(levels=c('CellConsensus','CellSpecific','DISEASE_RAW','DISEASE_ABS',
                                 'COLD_SHOCK', 'RNAi', 'CLATHRIN','DEGRADATION','ER2GOLGI',
                                 'FOLDING','LYSOSOME','MEMBRANE','NITRO','SORTING','TRAFFIC'))
M$score_column %<>% fct_recode('CELL_CONSENSUS'='CellConsensus', 'CELL_SPECIFIC'='CellSpecific')
M$db_metric = factor(paste(M$database %>% toupper, M$metric %>% toupper, sep=', '),
                     levels=rev(c('LINCS, XSUM', 'LINCS, KS', 'CMAP, XSUM', 'CMAP, KS')))
load('~/Desktop/Box/CFDR/cfdr/code/plot/2018_04_25/ES_all_drug_rankings.RData')

ggplot(M, aes(x=score_column, y=db_metric)) +
  geom_point(aes(colour=ES, size=-log10(pval))) +
  scale_color_gradient2(midpoint=0, low="darkgreen", mid="yellow", high="red") +
  scale_size(range = c(1, 10)) +
  xlab('') + ylab('') +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(PlotDir('ES_all_drug_rankings.pdf'), width=6.5, height=4)
#save(M, file=PlotDir('ES_all_drug_rankings.RData'))
