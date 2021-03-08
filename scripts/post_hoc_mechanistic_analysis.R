library(tidyr)
library(magrittr)
library(forcats)
library(data.table)

load(DataDir('lincs_signatures_bing_genes_tested_compounds.RData'))
load(DataDir('lincs_info_bing_genes_tested_compounds.RData'))

##### Load drug metadata
load(DataDir('tested_compound_data_and_metadata.RData'))
map = drugData$metadata

##### Identify 8 CFBE hits from this data
hit_names = c("15-DELTA PROSTAGLANDIN J2",
              "AZACITIDINE",
              "BRD-K94991378",
              "CD1530",
              "LDN-193189",
              "MD-II-008-P",
              "STROPHANTHIDIN",
              "WITHAFERIN-A")
m_hits = map[match(hit_names, map$name),]
hitSelect = sigInfo$pert_id %in% m_hits$pert_id

##### Identify non-hits (strict interval for 1 and 3uM, and allow decrease at 10uM)
lowDose = as.matrix(drugData$outcomes[,c('1uM_mean.c18','3uM_mean.c18','1uM_mean.dmso','3uM_mean.dmso')])
highDose = as.matrix(drugData$outcomes[,c('10uM_mean.c18','10uM_mean.dmso')])

idx1 = which(apply(lowDose, 1, function(x) InInterval(x, interval=c(0.85,1.15))))
idx2 = which(apply(highDose, 1, function(x) InInterval(x, interval=c(0.2, 1.15))))
idxNohit = intersect(idx1, idx2)
nohit_names = drugData$name[idxNohit]
m_nohits = map[match(nohit_names, map$name),]
nohitSelect = sigInfo$pert_id %in% m_nohits$pert_id

idxSelect = c(which(hitSelect), which(nohitSelect))
info = sigInfo[idxSelect,]
sigs = sigBing[idxSelect,]
info$hit = c(rep('hit', length(which(hitSelect))), rep('non-hit', length(which(nohitSelect))))

# Average per drug/cell combination
nGenes = ncol(sigs)
sigs$pert_id = info$pert_id
sigs$cell_id = info$cell_id
sigs = data.table(sigs)
meanSigs  = sigs[, lapply(.SD, mean), by=list(pert_id, cell_id), .SDcols=names(sigs)[1:nGenes]]

# Then normalize each row
M = meanSigs[,3:(nGenes+2)] %>% as.matrix
M = t(apply(M, 1, function(x) x/Norm2(x)))
M %<>% as.data.frame
M$pert_id = meanSigs$pert_id
M$cell_id = meanSigs$cell_id

# Then average across cell lines per drug
M = data.table(M)
finalSigs = M[, lapply(.SD, mean), by=list(pert_id), .SDcols=names(sigs)[1:nGenes]]
finalSigs = finalSigs %>% as.data.frame

# load gmt
load(DataDir('gmt_gene_sets_for_enrichment.RData'))

nDEG = 200
minGenes = 5
pThresh = 1

rownames(finalSigs) = finalSigs$pert_id
finalSigs = finalSigs[,2:ncol(finalSigs)]
all = list()
for(pert_id in rownames(finalSigs)){
  print(pert_id)
  x = finalSigs[pert_id,] %>% as.numeric %>% setNames(colnames(finalSigs))
  up = sort(-x)[1:nDEG] %>% names
  down = sort(x)[1:nDEG] %>% names

  UP = RunGsea(x, dz_genes_up=up, dz_genes_down=NULL, minGeneSetSize=minGenes, nPermute=1000, seed=123,
                          doGSOA=TRUE, doGSEA=FALSE, ListGSC = gmt, appendTerms=FALSE, preprocess=TRUE) %>%
    GetGSEAResults(colname='Adjusted.Pvalue', thresh=pThresh, analysis='GSOA', merge=TRUE) %>% FilterGSOARes %>%
    mutate(dirxn='UP')

  DOWN= RunGsea(x, dz_genes_up=down, dz_genes_down=NULL, minGeneSetSize=minGenes, nPermute=1000, seed=123,
                          doGSOA=TRUE, doGSEA=FALSE, ListGSC = gmt, appendTerms=FALSE, preprocess=TRUE) %>%
    GetGSEAResults(colname='Adjusted.Pvalue', thresh=pThresh, analysis='GSOA', merge=TRUE) %>% FilterGSOARes
  if(!is.null(DOWN)){DOWN$dirxn='DOWN'}

  name = map[which(map$pert_id==pert_id), 'name']
  all[[name]] = UP
  if(!is.null(DOWN)){all[[name]] %<>% rbind(DOWN)}
}

# Order by adjusted.pvalue
all %<>% lapply(function(x) x %>% arrange(Adjusted.Pvalue))
Write2XLS('all', file=PlotDir('CFBE_hits_and_nonhits_hypergeometric_enrichments.xlsx'))
save(all, file=PlotDir('CFBE_hits_and_nonhits_hypergeometric_enrichments.RData'))

# Grab the top K enrichments from each hit compound (hits only!)
K = 10
topResults = lapply(hit_names, function(nm){x=all[[nm]][1:K,]; x$drug=nm; return(x)})
R = Reduce('rbind', topResults) %>% unite('Term_with_dirxn', Term, dirxn, remove=FALSE) %>% as.data.frame
termCounts = R$Term_with_dirxn %>% table
selectedTerms = names(which(termCounts>1))

# Plot all the scores for these selected terms, not just if they are only in the top K for that drug
allResults = lapply(names(all), function(nm){x=all[[nm]]; x$drug=nm; x$hit=(nm %in% hit_names); return(x)})
R2 = Reduce('rbind', allResults) %>% unite('Term_with_dirxn', Term, dirxn, remove=FALSE) %>% as.data.frame
R2$drug %<>% as.factor %>% fct_reorder(R2$hit)
R3 = subset(R2, Term_with_dirxn %in% selectedTerms)
R3$Term_with_dirxn %<>% toupper
R3$drug %<>% fct_recode('15 DELTA PJ2'='15-DELTA PROSTAGLANDIN J2')
R3$logP = -log10(R3$Adjusted.Pvalue)
R3$logP[R3$logP > 20] = 20
print(ggplot(subset(R3, Adjusted.Pvalue<0.05), aes(y=Term_with_dirxn, x=drug)) +
        geom_point(aes(size=logP, colour=FoldEnrichment)) +
        scale_color_gradient2(low = "cyan", high = "red", mid='blue', midpoint=12) +
        scale_size(range = c(0.5, 4)) + theme_bw() + labs(x='', y='') +
        theme(text=element_text(size=10), axis.text.x=element_text(angle=45, hjust=1)))
ggsave(PlotDir('hit_v_nonhit_enrichment_map.svg'), width=7.5, height=6)

# Test whether scores on hits vs. non-hits are significantly different for any of these terms
out = list()
for(term in selectedTerms){
  print(term)
  A = subset(R3, Term_with_dirxn==term)
  a = subset(A, hit==TRUE)$FoldEnrichment
  b = subset(A, hit==FALSE)$FoldEnrichment
  out[[term]] = t.test(a,b)$p.value
}
adjp = p.adjust(unlist(out), method='BH') %>% sort

# For all significant pathways, make boxplot showing fold enrichment scores for hits vs. inactives
sigTerms = names(which(adjp<0.05)) %>% toupper
R4 = subset(R3, Term_with_dirxn %in% sigTerms)
R4 %<>% group_by(Term_with_dirxn) %>% mutate(median_val=median(FoldEnrichment, na.rm=TRUE))
R4$Term_with_dirxn %<>% as.factor %>% fct_reorder(R4$median_val, .desc=TRUE)
ggplot(R4, aes(x=Term_with_dirxn, y=FoldEnrichment, group=interaction(hit, Term_with_dirxn), fill=hit)) + labs(x='') +
  geom_boxplot() + theme_bw() + theme(text=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1))
ggsave(PlotDir('hit_v_nonhit_significant_pathways.svg'), width=5.5, height=7) # use this one for the plot
ggsave(PlotDir('hit_v_nonhit_significant_pathways.pdf'), width=5.5, height=7) # use this one for the labels
