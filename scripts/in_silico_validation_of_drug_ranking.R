set.seed(123)

# This script performs a sanity check on the four overall scores based on their rankings of a
# collection of positive control CFTR correctors.

for(database in c('cmap', 'lincs')){
  print(database)
  nPerm = 100000
  
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
  
  scoreSpecific='top_cell_score'
  scoreMulti='multi_cell_score'
  
  print('KS cell-specific score: ')
  res_ks_spec = ComputeCFEnrichmentScore(scores$ks[,c(compoundName, scoreSpecific)], scoreName=scoreSpecific,
                                         cfDrugs = cfDrugs, compoundName=compoundName)
  
  cat('\n\n')
  print('KS multi-cell score: ')
  res_ks_multi = ComputeCFEnrichmentScore(scores$ks[,c(compoundName, scoreMulti)], scoreName=scoreMulti,
                                          cfDrugs = cfDrugs, compoundName=compoundName)
  
  cat('\n\n')
  print('XSUM cell-specific score:')
  res_xsum_spec = ComputeCFEnrichmentScore(scores$xsum[,c(compoundName, scoreSpecific)], scoreName=scoreSpecific,
                                           cfDrugs = cfDrugs, compoundName=compoundName)
  
  cat('\n\n')
  print('XSUM multi-cell score:')
  res_xsum_multi = ComputeCFEnrichmentScore(scores$xsum[,c(compoundName, scoreMulti)], scoreName=scoreMulti,
                                            cfDrugs = cfDrugs, compoundName=compoundName)
  
  ### Compute empirical p-values of cell specific enrichment scores:
  nTotal = nrow(scores$ks)
  nIntersect = length(intersect(cfDrugs, scores$ks[,compoundName]))
  x = sapply(1:nPerm, function(i) GSEA.EnrichmentScore2(1:nTotal, sample(nTotal,nIntersect,replace=FALSE)))
  pval_ks_spec = length(which(x >= res_ks_spec$out$ES)) / nPerm
  pval_xsum_spec = length(which(x >= res_xsum_spec$out$ES)) / nPerm
  
  ### Compute empirical p-values of multi cell enrichment scores:
  multiCellDrugs = scores$ks[!is.na(scores$ks$multi_cell_score),compoundName]
  nDrugsMultiCell = length(multiCellDrugs)
  nCF = length(intersect(multiCellDrugs, cfDrugs))
  x = sapply(1:nPerm, function(i) GSEA.EnrichmentScore2(1:nDrugsMultiCell, sample(nDrugsMultiCell, nCF, replace=FALSE)))
  pval_ks_multi = length(which(x >= res_ks_multi$out$ES)) / nPerm
  pval_xsum_multi = length(which(x >= res_xsum_multi$out$ES)) / nPerm
  
  ### Plot
  yn = -0.5
  yoff = 0
  p_ks_spec = PlotCFEnrichmentHorizontal(res_ks_spec, scoreType=sprintf('%s KS/CellSpecific', toupper(database)), pval=pval_ks_spec, y_nudge=yn, y_offset=yoff)
  p_ks_multi = PlotCFEnrichmentHorizontal(res_ks_multi, scoreType=sprintf('%s KS/CellConsensus', toupper(database)), pval=pval_ks_multi, y_nudge=yn, y_offset=yoff)
  p_xsum_spec = PlotCFEnrichmentHorizontal(res_xsum_spec, scoreType=sprintf('%s XSum/CellSpecific', toupper(database)), pval=pval_xsum_spec, y_nudge=yn, y_offset=yoff)
  p_xsum_multi = PlotCFEnrichmentHorizontal(res_xsum_multi, scoreType=sprintf('%s XSum/CellConsensus', toupper(database)), pval=pval_xsum_multi, y_nudge=yn, y_offset=yoff)
  
  pdf(PlotDir(sprintf('overall_score_validations_%s.pdf', database)), width=14, height=12)
  MultiPlot(p_ks_spec, p_ks_multi, p_xsum_spec, p_xsum_multi, cols=2)
  dev.off()

}
