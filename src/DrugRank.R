
library(RankProd)
library(WGCNA)
source('src/ConnectivityScores.R')

ComputeAllScores = function(csMetric='ks', database='cmap', debug=FALSE, debugN=5, nCores=1,
                             calcAll=TRUE, calcSigs=calcAll, calcPathways=calcAll, 
                             saveRData=TRUE){

  #SetupParallel(nCores)

  if(debug){
    outDir = OutputDir(sprintf('drug_rankings/TEST_%s/', csMetric))
    rp.perm = 10
    nPerm = 30
  }else{
    outDir = OutputDir(sprintf('drug_rankings/%s/', csMetric))
    rp.perm = 500
    nPerm = 2000
  }

  if(!file.exists(outDir)){
    dir.create(outDir, recursive=T)
  }

  ####################################################################################################
  ### Compute scores for CF signatures
  ####################################################################################################

  # if(debug){
  #   load(OutputDir('sigs/meta3/allSigs_debug.RData'))
  # }else{
  #   load(OutputDir('sigs/meta3/allSigs.RData'))
  # }
  load(DataDir('cf_core_signatures.RData'))
  

  if(calcSigs){

    multiNScores = NULL

    for(sigName in names(allSigs)){
      print(toupper(sigName))
      individualScores = NULL
      for(N in names(allSigs[[sigName]])){
        querySig = list(UP=allSigs[[sigName]][[N]]$up$GeneID,DOWN=allSigs[[sigName]][[N]]$down$GeneID)
        print(system.time(individualScores[[N]] <- computeCsOnCmap(querySig, nPerm=nPerm, csMetric=csMetric,
                                                                   debug=debug, debugN=debugN, database=database)))
      }
      cat('merging individualScores...\n')
      print(system.time(multiNScores[[sigName]] <- MergeDfList(individualScores, allRows=TRUE, allSuffixes=TRUE,
                                                               by='compound')))
    }
    cat('merging overallScores...\n')
    print(system.time(sigScores <- MergeDfList(multiNScores, by='compound', allRows=TRUE, allSuffixes=TRUE)))
    save(sigScores, file=sprintf('%s/sigScores_%s.RData', outDir, database))
  }else{
    load(sprintf('%s/sigScores_%s.RData', outDir, database))
  }

  ####################################################################################################
  ### Compute scores for CF pathways
  ####################################################################################################

  #if(debug){
  #  load(DataDir('cf_genesets/metacore/metaCore_debug.RData')) # not sure where the debug file is. I think I just kept 2 or 3 of the pathways
  #}else{
    load(DataDir('cf_pathways.RData'))
  #}

  if(calcPathways){
    out = NULL
    for(pathwayName in names(metaCore)){
      print(toupper(pathwayName))
      querySig = list(UP=as.character(metaCore[[pathwayName]]))

      print(system.time(out[[pathwayName]] <- computeCsOnCmap(querySig, nPerm=nPerm, csMetric=csMetric,
                                                              debug=debug, debugN=debugN, database=database)))
    }
    pathwayScores = MergeDfList(out, by='compound', allRows=TRUE, allSuffixes=TRUE)
    save(pathwayScores, file=sprintf('%s/pathwayScores_%s.RData', outDir, database))
  }else{
    load(sprintf('%s/pathwayScores_%s.RData', outDir, database))
  }

  ####################################################################################################
  ### Merge all scores and add some finishing touches
  ####################################################################################################

  scoreList = list(sig=sigScores, pathway=pathwayScores)
  allScores = MergeDfList(scoreList, by='compound', allRows=TRUE, allSuffixes=TRUE)

  #minPNames = paste0('adjP.', c(paste0(names(metaCore), '.pathway')))
  allScores = ComputeFinalScores(allScores)
  allScores = PostProcessIds(allScores, database)
  #allScores = AppendAnnotation(allScores, database)

  cat('collapsing by top scoring cell type...\n')
  if(database=='lincs'){
    compoundName = 'pert_id'
  }else{
    compoundName = 'compound'
    if(csMetric=='xsum'){
      allScores = subset(allScores, cell_type %in% c('MCF7', 'PC3', 'HL60'))
    }
  }
  print(dim(allScores))
  print(system.time(allScores <- ComputeMultiCellScoreAndCollapseCellTypes(allScores, compoundName=compoundName, database=database)))
  print(dim(allScores))
  allScores = Factor2Char(allScores)
  rownames(allScores) = as.character(rownames(allScores))

  allScores = allScores[order(allScores[,'top_cell_score'], decreasing=TRUE),]
  simpleScores = SubsetAllScoresColumns(allScores, database=database)

  if(saveRData){
    save(sigScores, pathwayScores, allScores, simpleScores,
         file=sprintf('%s/allScores_%s.RData', outDir, database))
  }
}


AddRank = function(A, colName='meta_rank'){
  A$rank = rank(A[,colName], ties.method='min')
  return(A)
}

AddSign = function(A, colName='mean_rankpt_4'){
  A2 = sign(A[,which(grepl(colName, names(A)))])
  if(is.vector(A2)){
    A$sign = A2
  }else{
    stopifnot(ncol(A2)==3)
    print(colnames(A2))
    A$sign = sapply(1:nrow(A2), function(i) GetSign(A2[i,]))
  }
  A = MoveColumn(A, 'sign', newPosition=1+which(names(A) == 'meta_rank'))
  return(A)
}

GetSign = function(x, l=length(x)){
  if(sum(x) == l){
    out = '+'
  }else if (sum(x) == -l){
    out = '-'
  }else{
    out = NA
  }
  return(out)
}

AppendAnnot = function(A, annot, colName='compound'){
  idxMatch = match(as.character(A[[colName]]), as.character(annot[[colName]]))
  stopifnot(length(which(!is.na(idxMatch))) > 0)
  print(sprintf('%d compounds annotated', length(which(!is.na(idxMatch)))))
  idxCompound = which(names(A) == colName)
  return(cbind(A, annot[idxMatch,-idxCompound]))
}

StripAnnot = function(A, annot){
  idxAnnot = unlist(lapply(names(annot), function(x) which(grepl(x, names(A)))))
  if(length(idxAnnot) == 0){
    warning('did not find any column names matching annotation column names')
    out = A
  }else{
    out = A[,-idxAnnot]
  }
  return(out)
}

MapCMAP2LINCS = function(C=GetMAPAnnot('cmap'), L=GetMapAnnot('lincs'), collapse=',', minColumns=2){
  overlap = intersect(names(C), names(L))
  stopifnot(length(overlap) > 0)
  idx = C[,overlap]
  for(colName in overlap){
    print(colName)
    idx[,colName] = MatchAll(C[,colName], L[,colName], collapse=',') #match(C[,colName], L[,colName], incomparables=NA)
  }
  idx$consensus = unlist(lapply(1:nrow(idx), function(i)
    GetConsensusMapping(idx[i,], collapse=collapse, minColumns=minColumns)))
  cat(sprintf('%d compounds mapped.\n', length(na.omit(idx$consensus))))
  return(idx)
}

GetConsensusMapping = function(row, collapse='|', minColumns=1){
  subrow = row[!is.na(row)]
  if(length(subrow) == 0){
    out = NA
  }else{
    tmp = paste(subrow, collapse=collapse)
    split = paste0('([', collapse, '])')
    outList = unlist(strsplit(tmp, split=split))
    out = GetMultiples(outList, minColumns)
    if(length(out)==0){
      out = NA
    }else{
      out = paste(out, collapse=collapse)
    }
  }
  return(out)
}

GetMapAnnot = function(databaseName, full=FALSE){
  if(tolower(databaseName) == 'cmap'){
    A = GetCMAPAnnot()
    A$unique_id = A$compound
    colNames = GetCMAPColNames()
  }else if(tolower(databaseName) == 'lincs'){
    A = GetLINCSAnnot(full=full)
    load('data/lincs/pertInfo/pertInfo.RData')
    P = pertInfo[,c('inchi_string', 'pubchem_cid', 'pert_id', 'pert_iname')]
    names(P) = c('inchi', 'pubchem', 'unique_id', 'compound')
    rownames(P) = P$unique_id
    P = Factor2Char(P)
    A = merge(A, P, by='unique_id', all=TRUE)
    colNames = GetLINCSColNames()
  }else{
    stop('unexpected databaseName')
  }
  from = colNames$from
  to = colNames$to
  compoundSave = A$compound
  for(i in 1:length(from)){
    print(to[i])
    A = ChangeColumnName(A, from[i], to[i])
    A[to[i]] = unlist(lapply(1:nrow(A), function(x) toupper(A[x, to[i]]) ))
  }
  rownames(A) = A[,'unique_id']
  A$COMPOUND = gsub('[[:punct:]]', '', A$COMPOUND)
  A$COMPOUND = gsub(' ', '', A$COMPOUND)
  A$compound = compoundSave
  return(FixNAStrings(A[, c(to, 'compound')]))
}

GSEA.EnrichmentScore = function(drug.list, drug.set) {
  correl.vector = NULL
  weighted.score.type = 0
  # Computes the weighted GSEA score of drug.set in drug.list.
  # The weighted score type is the exponent of the correlation
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is
  # necessary to input the correlation vector with the values in the same order as in the drug list.
  #
  # Inputs:
  #   drug.list: The ordered drug list (e.g. integers indicating the original position in the input dataset)
  #   drug.set: A drug set (e.g. integers indicating the location of those drugs in the input dataset)
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the drugs in the drug list
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1)
  #   arg.ES: Location in drug.list where the peak running enrichment occurs (peak of the "mountain")
  #   RES: Numerical vector containing the running enrichment score for all locations in the drug list
  #   tag.indicator: Binary vector indicating the location of the drug sets (1's) in the drug list
  #
  # This adapted from GSEA.1.0.R from the Broad Institute.
  tag.indicator = sign(match(drug.list, drug.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator = 1 - tag.indicator
  N = length(drug.list)
  Nh = length(drug.set)
  Nm =  N - Nh
  if (weighted.score.type == 0) {
    correl.vector = rep(1, N)
  }
  alpha = weighted.score.type
  correl.vector = abs(correl.vector**alpha)
  sum.correl.tag    = sum(correl.vector[tag.indicator == 1])
  norm.tag    = 1.0/sum.correl.tag
  norm.no.tag = 1.0/Nm
  RES = cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  max.ES = max(RES)
  min.ES = min(RES)
  if (max.ES > - min.ES) {
    #      ES = max.ES
    ES = signif(max.ES, digits = 5)
    arg.ES = which.max(RES)
  } else {
    #      ES = min.ES
    ES = signif(min.ES, digits=5)
    arg.ES = which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}

GSEA.EnrichmentScore2 = function(gene.list, gene.set, weighted.score.type = 0, correl.vector = NULL) {
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the
  # observed one.
  # The weighted score type is the exponent of the correlation
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset)
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.

  N = length(gene.list)
  Nh = length(gene.set)
  Nm =  N - Nh

  loc.vector = vector(length=N, mode="numeric")
  peak.res.vector = vector(length=Nh, mode="numeric")
  valley.res.vector = vector(length=Nh, mode="numeric")
  tag.correl.vector = vector(length=Nh, mode="numeric")
  tag.diff.vector = vector(length=Nh, mode="numeric")
  tag.loc.vector = vector(length=Nh, mode="numeric")

  loc.vector[gene.list] = seq(1, N)
  tag.loc.vector = loc.vector[gene.set]

  tag.loc.vector = sort(tag.loc.vector, decreasing = F)

  if (weighted.score.type == 0) {
    tag.correl.vector = rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector = correl.vector[tag.loc.vector]
    tag.correl.vector = abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector = correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector = abs(tag.correl.vector)
  } else {
    tag.correl.vector = correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector = abs(tag.correl.vector)
  }

  norm.tag = 1.0/sum(tag.correl.vector)
  tag.correl.vector = tag.correl.vector * norm.tag
  norm.no.tag = 1.0/Nm
  tag.diff.vector[1] = (tag.loc.vector[1] - 1)
  tag.diff.vector[2:Nh] = tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector = tag.diff.vector * norm.no.tag
  peak.res.vector = cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector = peak.res.vector - tag.correl.vector
  max.ES = max(peak.res.vector)
  min.ES = min(valley.res.vector)
  ES = signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)

  return(list(ES = ES))

}

GetKnownCorrectors = function(){
  return(c('chlorzoxazone','chloramphenicol', 'ouabain', 'digoxin', 'valproic acid', 'pizotifen',
           'biperiden', 'MG-132', 'neostigmine bromide', 'midodrine', 'vorinostat',
           'trichostatin A', 'scriptaid', 'MS-275', 'glafenine', 'curcumin', 'ibuprofen',
           'glycerol', 'pyridostigmine', 'vardenafil', 'tadalafil', 'KM11060', 'RDR01752', 'miglustat'))
}

GetKnownCorrectors2 = function(addNewCompounds=TRUE){
  # the integers are pubchem ids
  correctors = c(
    'VALPROIC ACID' = 3121,
    'CHLORZOXAZONE' = 2733,
    'NEOSTIGMINE BROMIDE' = 8246,
    'MIDODRINE' = 4195,
    'TRICHOSTATIN A' = 444732,
    'MS-275' = 4261,
    'GLAFENINE' = 3474,
    'IBUPROFEN' = 3672,
    'GLYCEROL' = 753,
    'PYRIDOSTIGMINE' = 4991,
    'VARDENAFIL' = 110634,
    'TADALAFIL' = 110635,
    'KM11060' = 1241327,
    'RDR01752' =	9566157,
    'MIGLUSTAT' = 51634,
    'CHLORAMPHENICOL' = 5959,
    'OUABAIN' = 439501,
    'DIGOXIN' = 2724385,
    'PIZOTIFEN' = 27400,
    'BIPERIDEN' = 2381,
    'MG-132' = 462382,
    'VORINOSTAT' = 5311,
    'SCRIPTAID' = 5186,
    'CURCUMIN' = 969516)
  if(addNewCompounds){
    correctors = c(correctors,
    '15-DELTA PROSTAGLANDIN J2' = 5311211,
    'AZACITIDINE' = 23760137,
    'BRD-K94991378' = 6731789,
    'CD1530' = 24868309,
    'LDN-193189' = 25195294,
    'MD-II-008-P' = 60194076,
    'STROPHANTHIDIN' = 6185,
    'WITHAFERIN-A' = 73707417)
  }
  return(correctors)
}

GetPreviouslyTested = function(){
  return(c(GetKnownCorrectors(), 'genistein', 'digitoxin','ouabagenin'))
}

ComputeES = function(drugScores, drugSet, scoreName='score', compoundName='compound', higherIsBetter=TRUE){
  if(higherIsBetter){
    drug.list = order(-drugScores[,scoreName])
  }else{
    drug.list = order(drugScores[,scoreName])
  }

  drugSetIdx = which(drugScores[,compoundName] %in% drugSet)
  stopifnot(length(drugSetIdx) == length(drugSet))
  out.all = GSEA.EnrichmentScore2(drug.list, drugSetIdx)
  return(out.all)
}


# This one calls GSEA.EnrichmentScore which computes the score as well as outputs for plotting
ComputeCFEnrichmentScore = function(drugScores, scoreName='score', compoundName='compound',
                                    cfDrugs = GetKnownCorrectors(),
                                    higherIsBetter=TRUE, a=NULL,justScore=TRUE, abs=FALSE, nPerm=100000, seed=123){

  drugScores = na.omit(drugScores[,c(compoundName, scoreName)])

  if(higherIsBetter){
    drug.list = order(-drugScores[,scoreName])
  }else{
    drug.list = order(drugScores[,scoreName])
  }

  drug.set = which(drugScores[,compoundName] %in% cfDrugs)
  if(length(drug.set) != length(cfDrugs)){
    warning(sprintf('found %d out of %d known correctors', length(drug.set), length(cfDrugs)))
  }
  cfDrugs = drugScores[drug.set,compoundName] # restrict to drugs in database
  out = GSEA.EnrichmentScore(drug.list, drug.set)

  dS = drugScores[drug.list,]
  cfRanks = match(cfDrugs, dS[,compoundName])
  b = data.frame(compound=cfDrugs, rank=cfRanks)
  if(!is.null(a)){
    print(sprintf('a=%f',a))
  }
  print(b[order(b$rank),])

  ### Compute empirical p-value
  nTotal = nrow(drugScores)
  nIntersect = length(intersect(cfDrugs, drugScores[,compoundName]))
  set.seed(seed)
  x = sapply(1:nPerm, function(i) GSEA.EnrichmentScore2(1:nTotal, sample(nTotal, nIntersect, replace=FALSE)))
  p = sapply(out$ES, function(es) length(which(x >= es)) / nPerm)

  if(!is.null(na.action(drugScores))){
    result = list(out=out, nDrugs=nrow(drugScores), nCfDrugs=length(cfDrugs), cfRanks=b, pval=p)
  }else{
    result = list(out=out, cfRanks=b, pval=p)
  }
  return(result)
}


# Input should be result$out from above function
PlotCFEnrichment = function(res, scoreType, xlim=c(1, 1309), ylim=c(-0.5, 0.6), print=FALSE, pval){
  library(ggrepel)
  lwd = 1
  n = length(res$out$RES)
  df = data.frame(idx=1:n, RES=res$out$RES)
  df = merge(df, res$cfRanks, by.x='idx', by.y='rank', all=TRUE)
  df$compound_w_rank = ifelse(is.na(df$compound), NA, sprintf('%s-%d', df$compound, df$idx))
  title = sprintf('%s Score\n Max ES=%0.2f (p=%0.1e)', scoreType, res$out$ES, pval)
  p = ggplot(data=df, mapping=aes(x=idx, y=RES)) +
    geom_line(col='black', lwd=lwd) +
    geom_point(aes(x=res$out$arg.ES, y=res$out$ES), shape=8, size=4) +
    geom_hline(col='black', yintercept=0, lwd=lwd/1.5, linetype='dashed') +
    geom_text_repel(data=df, aes(x=idx, y=-0.35, label=compound_w_rank), color='black',
                    size=3, segment.size = 0.25, segment.color='black', box.padding = unit(0, "lines"),
                    point.padding = unit(0, "lines")) +
    xlim(xlim) + ylim(ylim) +
    xlab('Drug Ranking') + ylab('Running Enrichment Score (ES)') +
    ggtitle(title) +
    theme_classic() +
    theme(text=element_text(size=14),  plot.title=element_text(hjust = 0.5)) +
    coord_flip() +
    scale_x_reverse()
  if(print){print(p)}
  return(p)
}

# Input should be result$out from above function
PlotCFEnrichmentHorizontal = function(res, scoreType, ylim=c(-0.6, 0.6), y_offset=0, y_nudge=-0.2, print=FALSE, pval){

  library(ggrepel)

  nDrugs = length(res$out$RES)
  fudge_x = nDrugs / 20
  lwd = 1
  n = length(res$out$RES)
  df = data.frame(idx=1:n, RES=res$out$RES)
  df = merge(df, res$cfRanks, by.x='idx', by.y='rank', all=TRUE)
  df$compound_w_rank = ifelse(is.na(df$compound), NA, sprintf('%s-%d', df$compound, df$idx))
  title = sprintf('%s Score\n Max ES=%0.2f (p=%0.1e)', scoreType, res$out$ES, pval)
  x_max = res$out$arg.ES
  y_max = res$out$ES
  p = ggplot(data=df, mapping=aes(x=idx, y=RES)) +
    geom_line(col='black', lwd=lwd) +

    # max point and corresponding line
    geom_point(aes(x=x_max, y=y_max), shape=8, size=4, color='red') +
    geom_line(col='red', data=data.frame(x=c(-fudge_x, x_max), y=c(y_max, y_max)), aes(x=x, y=y), lwd=lwd/1.5) +

    # baseline: y=0
    geom_hline(col='black', yintercept=0, lwd=lwd/2) + #, linetype='dashed') +

    # show location of known CF compounds
    geom_point(data=subset(df, !is.na(compound)), aes(x=idx, y=y_offset), color = '#00CCFF') +
    geom_text_repel(aes(x=idx, y=y_offset, label=compound_w_rank),
      force        = 0.5,
      nudge_y      = y_nudge,
      direction    = 'x',
      angle        = 90,
      vjust        = 1,
      segment.size = 0.2) +

    xlim(c(-fudge_x, nDrugs)) + ylim(ylim) +
    xlab('Drug Ranking') + ylab('Running Enrichment Score (ES)') +
    ggtitle(title) +
    theme_classic() +
    theme(text=element_text(size=14),  plot.title=element_text(hjust = 0.5))

  if(print){print(p)}
  return(p)
}

ComputeCFEnrichmentScoreRP = function(RP.out){
cfDrugs = dfrrectors()
  drug.list = order(RP.out$pfp[[2]])
  drug.set = which(rownames(RP.out$pfp) %in% cfDrugs)
  stopifnot(length(drug.set) == length(cfDrugs))
  out = GSEA.EnrichmentScore2(drug.list, drug.set)
  return(out$ES)
}

# It is okay for ranks to have NA's
# Either percent or hingePoint should be defined
RankTransformHinge = function(ranks, hingePoint=length(na.omit(ranks)), percentile=NULL){
  n = length(na.omit(ranks))
  if(!is.null(percentile)){
    hingePoint = floor(n*percentile)
  }
  hingePoint = min(hingePoint, length(na.omit(ranks)))
  rr = na.omit(ranks)

  hingeScore = rr
  hingeScore[rr > hingePoint] = 0
  hingeScore[rr <= hingePoint] = (hingePoint+1 - hingeScore[rr <= hingePoint]) / hingePoint

  hingeScoreFull = ranks
  hingeScoreFull[which(!is.na(ranks))] = hingeScore

  return(hingeScoreFull)
}

RankTransformHingeAll = function(drugScores, hingePoints, percentile){
  #hP = hingePoints$cmap
  pStr = paste0('.perc', 100*percentile)
  drugScores[,paste0('hinge.rank.cold.cmap', pStr)] = RankTransformHinge(drugScores$rank.cold.cmap, percentile=percentile)
  drugScores[,paste0('hinge.rank.rnai.cmap', pStr)] = RankTransformHinge(drugScores$rank.rnai.cmap, percentile=percentile)
  drugScores[,paste0('hinge.rank.pighuman.cmap', pStr)] = RankTransformHinge(drugScores$rank.pighuman.cmap, percentile=percentile)

  #hP = hingePoints$lincs
  drugScores[,paste0('hinge.rank.cold.lincs', pStr)] = RankTransformHinge(drugScores$rank.cold.lincs, percentile=percentile)
  drugScores[,paste0('hinge.rank.rnai.lincs', pStr)] = RankTransformHinge(drugScores$rank.rnai.lincs, percentile=percentile)
  drugScores[,paste0('hinge.rank.pighuman.lincs', pStr)] = RankTransformHinge(drugScores$rank.pighuman.lincs, percentile=percentile)

  #drugScores[,paste0('hinge.rank.total.cmap', hingePoints$cmap, '.lincs', hingePoints$lincs)] =
  drugScores$hinge.rank.no.norm =
    apply(drugScores[,grepl('hinge.rank',names(drugScores))], 1, function(x) sum(na.omit(x)))

  return(drugScores)
}

ComputeRank = function(cmap, sigName, takeAbs=TRUE){
  data = cmap[,grepl(paste0('score...',sigName),names(cmap))]
  stopifnot(ncol(data) == 3)

  if(takeAbs){
    data = abs(data)
  }

  for(i in 1:3){
    data[,i] = rank(-data[,i])
  }
  out = apply(data, 1, GeometricMean)
  return(out)
}

TestGSEA.EnrichmentScore = function(){
  drug.list = order(drugScores$score1.rank)
  drug.set = which(drugScores$score1.rank <= 10)
  out = GSEA.EnrichmentScore(drug.list, drug.set)
  stopifnot(out$ES == 1)

  drug.set = which(drugScores$score1.rank <= 20 & drugScores$score1.rank > 10)
  out = GSEA.EnrichmentScore(drug.list, drug.set)
  stopifnot(abs(out$ES - 0.9974) < 0.001)

  drug.set = which(drugScores$score1.rank > (nrow(drugScores)-10))
  out = GSEA.EnrichmentScore(drug.list, drug.set)
  stopifnot(out$ES == -1)
}

TestComputeCFEnrichmentScore = function(){
  cfDrugs = toupper(dfrrectors())
  otherDrugs = setdiff(drugScores$compound, cfDrugs)
  idx = which(drugScores$compound %in% otherDrugs)

  drug.list.A = drugScores[idx,]
  drug.list.B = drugScores[-idx,]
  drug.list = rbind(drug.list.A, drug.list.B)
  drug.list$score = 1:nrow(drug.list)

  stopifnot(ComputeCFEnrichmentScore(drug.list, 'score')==1)
  drug.list$score[(length(otherDrugs)+1):nrow(drug.list)] = -1:-10

  stopifnot(ComputeCFEnrichmentScore(drug.list, 'score')==-1)

  ord = gtools::permute(1:nrow(drug.list))
  stopifnot(ComputeCFEnrichmentScore(drug.list[ord,], 'score')==-1)
}

SampleCMAP = function(n, nMax=1309){
  return(floor(runif(n, 1,nMax)))
}

AppendMetaRank = function(df, colName='score', abs=FALSE, append=TRUE, num.perm=200){
  scoreNames = names(df)[grepl(colName, names(df))]
  print(scoreNames)
  origin = rep(1, length(scoreNames))
  class = rep(1, length(scoreNames))
  data = df[,scoreNames]
  if(is.logical(abs)){
    if(abs){
      data = abs(data)
    }
  }else{
    names = names(data)[grepl(abs, names(data))]
    if(length(names) > 0){
      print(sprintf('taking absolute value of %s', names))
    }
    data[,names] = abs(data[,names])
  }

  RP = RPadvance(data, cl=class, origin=origin, logged=F, gene.names=df$compound,
                  plot=FALSE, rand=42, num.perm=num.perm)

  out = data.frame(pvalue.consensus=RP$pval[[2]],
                    fdr.consensus=RP$pfp[[2]],
                    rank.consensus=RP$RPrank[[2]])
  if(append){
    out = data.frame(df, out)
  }

  return(out)
}

AppendMinP = function(df, colNames, varName, nCol=length(colNames), minus=0){
  P = df[,colNames]
  stopifnot(ncol(P) == nCol)
  df[[varName]] = apply(P, 1, min)
  df[[paste0(varName, '.source')]] = apply(P, 1, function(x)
    {y = unlist(strsplit(names(P)[which(x == min(x))], split='[.]'));
     y[length(y) - minus]})
  return(df)
}

TestDebugQuery = function(database, csMetric, newSigScores, newPathwayScores, newAllScores){
  load(OutputDir(sprintf('drug_rankings/meta3/internal/newPRLs/debugSave/%s/allScores_%s.RData', csMetric, database)))
  stopifnot(identical(sigScores, newSigScores))
  stopifnot(identical(pathwayScores, newPathwayScores))
  stopifnot(identical(allScores, newAllScores))
}

LogP = function(p, replaceZeros=Inf){
  out = log10(p)

  if(is.infinite(replaceZeros)){
    replaceZeros = min(out[!is.infinite(out)])
  }

  if(is.infinite(replaceZeros)){
    warning('no non-infinite values for LogP, hence not doing replacement')
  }else{
    out[is.infinite(out)] = replaceZeros
  }
  return(out)
}

SubsetAllScoresColumns = function(allScores, database, extras=TRUE, orderBy=NULL){

  idxScore = which(grepl('top_cell_score|multi_cell_score', names(allScores)))

  idxTested = which(grepl('previously.tested', names(allScores)))

  if(database == 'lincs'){
    idxCrossRef = which(grepl('cmap', names(allScores)))
    idxId = which(names(allScores) %in% c('pert_id', 'pert_iname', 'unique_id','cell_type'))
  }else if(database == 'cmap'){
    idxCrossRef = which(grepl('lincs', names(allScores)))
    idxId = which(names(allScores) %in% c('compound'))
  }else{
    stop('unexpected database name')
  }

  if(extras){
    idxPmin = which(grepl('pmin|fdr.consensus', names(allScores)))
    idxAdjP = which(grepl('adjPvalue.*pathway', names(allScores)))
    idxReason = which(grepl('reasonForFiltering', names(allScores)))
  }else{
    idxPmin = c()
    idxAdjP = c()
    idxReason = c()
  }

  toKeep = c(idxId, idxScore, idxTested, idxCrossRef, idxPmin, idxAdjP, idxReason)
  allScores = allScores[,toKeep]
  if(!is.null(orderBy)){
    stopifnot(orderBy %in% names(allScores))
    allScores = allScores[order(allScores[[orderBy]]),]
  }
  return(allScores)
}


ComputeFinalScores = function(allScores){
  allScores = ComputeBestNormSigScore(allScores)
  allScores = ComputeBestNormPathwayScore(allScores)
  allScores$top_cell_score = allScores$best.normScore.pathway + allScores$best.normScore.sig
  allScores = MoveColumn(allScores, colName = 'top_cell_score', newPosition = 2)
  return(allScores)
}

ComputeBestNormPathwayScore = function(allScores){
  col = grepl(glob2rx(sprintf('normScore*pathway')), names(allScores))
  pathwayScores = abs(allScores[,col])
  pathwayNames = sapply(names(pathwayScores), function(x) unlist(strsplit(x, split='[.]'))[2])

  allScores$best.normScore.pathway = apply(pathwayScores, 1, max)
  allScores$best.pathway = apply(pathwayScores, 1, function(x) pathwayNames[which.max(x)])
  return(allScores)
}

ComputeBestNormSigScore = function(allScores){
  allScores = ComputeMeanNormSigScores(allScores)
  sigScores = allScores[,c('normScore.cold.sig', 'normScore.rnai.sig', 'normScore.pighuman.sig')]
  sigNames = c('cold', 'rnai', 'pighuman')

  allScores$best.normScore.sig = apply(sigScores, 1, max)
  allScores$best.sig = apply(sigScores, 1, function(x) sigNames[which.max(x)])
  return(allScores)
}

ComputeMeanNormSigScores = function(allScores, colName='normScore'){
  for(sigName in c('cold', 'rnai', 'pighuman')){
    rx = glob2rx(sprintf('%s*%s.sig', colName, sigName))
    col = grepl(rx, names(allScores))
    normScores = allScores[,col]
    if(ncol(normScores) != 3){
      warning(sprintf('ncol(normScores) = %d', ncol(normScores)))
    }
    mn = apply(normScores, 1, mean)
    if(sigName == 'pighuman'){
      mn = abs(mn)
    }
    allScores[[sprintf('normScore.%s.sig', sigName)]] = mn
  }
  return(allScores)
}

PostProcessIds = function(allScores, database){
  allScores = Factor2Char(allScores)
  allScores = ChangeColumnName(allScores, from='compound', to='unique_id')
  allScores$cell_type = sapply(allScores$unique_id, function(x) unlist(strsplit(x, split='_'))[1])
  if(database=='lincs'){
    allScores$pert_id = sapply(allScores$unique_id, function(x) unlist(strsplit(x, split='_'))[2])
  }else if(database=='cmap'){
    nChar = sapply(allScores$cell_type, nchar)
    allScores$compound = sapply(1:nrow(allScores), function(i) substr(allScores$unique_id[i],
                               start = nChar[i]+2, stop=nchar(allScores$unique_id[i])))
  }
  allScores = RemoveDfColumns(allScores, columnNames='unique_id')
  return(allScores)
}

AppendAnnotation = function(allScores, database){
  load(DataDir('cmap_lincs_compound_mapping.RData'))
  if(database == 'lincs'){
    load(DataDir('lincs_compound_metadata.RData'))

    #removing Generic Name because it causes problems sometimes when writing to xls
    lincsAnnot = lincsAnnot[,names(lincsAnnot) != 'Generic_Name.inHouse']
    allScores = merge(allScores, lincsAnnot, by=c('pert_id', 'pert_iname'), all.x=TRUE, all.y=FALSE)

    pMap = pairwiseMap[,c('lincs.id', 'cmap.name', 'matching.fields')]
    pMap = aggregate(pMap, by=list(pMap$lincs.id), paste, collapse='|')
    pMap = RemoveDfColumns(pMap, 'lincs.id')
    names(pMap) = c('pert_id', 'matching.cmap.names', 'cmap.matching.fields')

    allScores = merge(allScores, pMap, by='pert_id', all.x=TRUE, all.y=FALSE)

    rownames(allScores) = allScores$pert_id
  }else if(database == 'cmap'){
    cmapAnnot = RemoveDfColumns(GetCMAPAnnot(), c('SMILES', 'unique_id'))
    allScores = merge(allScores, cmapAnnot, by='compound', all.x=TRUE, all.y=FALSE)

    pMap = pairwiseMap[,c('lincs.id', 'lincs.name', 'cmap.name', 'matching.fields')]
    pMap = aggregate(pMap, by=list(pMap$cmap.name), paste, collapse='|')
    pMap = RemoveDfColumns(pMap, 'cmap.name')
    names(pMap) = c('compound', 'matching.lincs.ids', 'matching.lincs.names', 'lincs.matching.fields')
    allScores = merge(allScores, pMap, by='compound', all.x=TRUE, all.y=FALSE)

    rownames(allScores) = allScores$compound
  }

  dbp = GetDbp(database)
  tested = StandardizeCompoundNames(GetPreviouslyTested())
  drugNames = StandardizeCompoundNames(allScores[[dbp$compoundName]])
  idx = which(drugNames %in% tested)
  allScores$previously.tested = FALSE
  allScores$previously.tested[idx] = TRUE

  return(allScores)
}

RemovePreviouslyTested = function(allScores, colName='compound'){
  tested = GetPreviouslyTested()
  idx = which(allScores[,colName] %in% tested)
  if(length(idx) > 0){
    remainingScores = allScores[-idx,]
    removedScores = allScores[idx,]
    removedScores$reasonForFiltering = 'previously tested'
  }else{
    remainingScores = allScores
    removedScores = allScores[setdiff(1,1),] # empty data frame with appropriate column names
  }

  return(list(remainingScores=remainingScores, removedScores=removedScores))
}

# This one uses collapseRows, which is very slow and awkward
ComputeMultiCellScoreAndCollapseCellTypes_OLD = function(allScores, scoreName='top_cell_score', compoundName='pert_id', database='lincs'){
  # Define rownames to be unique so that collapseRows can work properly
  rownames(allScores) = paste(allScores$cell_type, allScores[[compoundName]], sep='_')

  # Copy column to be collapsed so that collapseRows can work properly (this is a hack in order to use collapseRows to simply select the row with the max 'top_cell_score' per drug)
  data = allScores[,c(scoreName, scoreName)]

  # Run collapseRows, which selects the cell type for each drug with the strongest sig + pathway connectivity score (i.e. top_cell_score)
  collapse.out = collapseRows(data, rowGroup=allScores[[compoundName]], rowID=rownames(allScores), method='MaxMean')

  # Subset allScores to choose the relevant rows based on top_cell_score
  winners = allScores[names(which(collapse.out$selectedRow)),]
  stopifnot(!anyDuplicated(winners[[compoundName]]))

  # Define celltypes_and_rank column
  cellAndRank = sprintf('%s(%d)', allScores$cell_type, rank(-allScores[[scoreName]]))
  ag = aggregate(cellAndRank, by=list(allScores[[compoundName]]), paste, collapse='; ')
  names(ag) = c(compoundName, 'celltypes_and_rank')


  justRank = rank(-allScores[[scoreName]])
  nRow = nrow(allScores)
  rankAg = aggregate(justRank, by=list(allScores[[compoundName]]), GSEA.EnrichmentScore2, gene.list=1:nRow)
  names(rankAg) = c(compoundName, 'multi_cell_score')
  rankAg$multi_cell_score = unlist(rankAg$multi_cell_score)

  out = merge(winners, ag, by=compoundName)
  out = merge(out, rankAg, by=compoundName)
  return(out)
}

# The new version uses dplyr
ComputeMultiCellScoreAndCollapseCellTypes = function(allScores, scoreName='top_cell_score', compoundName='pert_id', database='lincs'){

  # Compute multi-cell score
  justRank = rank(-allScores[[scoreName]])
  nRow = nrow(allScores)
  rankAg = aggregate(justRank, by=list(allScores[[compoundName]]), GSEA.EnrichmentScore2, gene.list=1:nRow)
  names(rankAg) = c(compoundName, 'multi_cell_score')
  rankAg$multi_cell_score = unlist(rankAg$multi_cell_score)

  # Get top cell score per compound
  allScores %<>% ChangeColumnName(from=scoreName, to='top_cell_score')
  allScores %<>% group_by_(compoundName) %>% dplyr::filter(top_cell_score == max(top_cell_score))

  return(merge(allScores, rankAg, by=compoundName))
}

# expand output of MapCMAP2LINCS to list pairwise CMAP/LINCS relationships
ExpandFullMap = function(fullMap, LINCS){

  # the actual matching is done on upper case version of compound, but let's call it 'name' vs.
  # 'COMPOUND' in the output
  fullMap = RemoveDfColumns(fullMap, 'compound')
  fullMap = ChangeColumnName(fullMap, from='COMPOUND', to='name')

  cmap.name = NULL
  lincs.name = NULL
  lincs.id = NULL
  matching.fields = NULL

  consensusChar = Str2Vec(fullMap$consensus, split='[,]')
  consensus = lapply(consensusChar, as.integer)

  for(i in 1:length(consensus)){

    a = Str2Vec(fullMap[i,-ncol(fullMap)])
    b = lapply(a, as.integer)

    for(idx in na.omit(consensus[[i]])){
      cmap.name = c(cmap.name, rownames(fullMap)[i])
      lincs.name = c(lincs.name, lincs.name=LINCS$compound[idx])
      lincs.id = c(lincs.id, lincs.id=LINCS$unique_id[idx])
      matching.fields = c(matching.fields, paste0(names(which(
        unlist(lapply(b, function(x) idx %in% x)))), collapse=','))
    }
  }

  pMap = data.frame(cmap.name=cmap.name, lincs.name=lincs.name, lincs.id=lincs.id, matching.fields=matching.fields)
  pMap = ExpandMatchingFields(Factor2Char(pMap))

  return(pMap)
}

ExpandMatchingFields = function(pMap, fields=c('name','chembl','chebi','drugbank','pubchem','mesh','kegg','smiles','inchi')){
  for(field in fields){
    pMap[[field]] = grepl(field, pMap$matching.fields)
  }
  return(pMap)
}

StandardizeCompoundNames = function(names){
  names = toupper(names)
  names = gsub('[[:punct:]]', '', names)
  names = gsub(' ', '', names)
  return(names)
}

GetATC = function(cmap, colName, annot=GetCMAPAnnot(), topN=1309){
  if(!is.null(colName)){
    idx = which(rank(cmap[,colName]) <= topN)
  }else{
    idx = sample(1:1309, size=topN)
  }

  compounds = cmap$compound[idx]
  ATClong = annot$atc_code[match(compounds, annot$compound)]
  ATClong[which(ATClong == 'NA')] = NA
  ATC = substr(ATClong, 1, 1)
  ATC[which(is.na(ATC))] = 'unknown'
  ATCfactor = as.factor(ATC)
  return(ATCfactor)
}

AppendPertIName = function(allScores){
  load(DataDir('lincs_compound_metadata.RData'))
  allScores = merge(allScores, lincsAnnot[,c('pert_id', 'pert_iname')], by='pert_id', all=FALSE)
  return(allScores)
}

LoadScores = function(baseDir=OutputDir('drug_rankings/')){
  scores = NULL

  load(paste0(baseDir, 'ks/allScores_cmap.RData'))
  allScores$COMPOUND = toupper(allScores$compound)
  rownames(allScores) = allScores$compound
  scores$cmap_ks = allScores

  load(paste0(baseDir, 'ks/allScores_lincs.RData'))
  allScores = AppendPertIName(allScores)
  allScores$PERT_INAME = toupper(allScores$pert_iname)
  rownames(allScores) = allScores$pert_id
  scores$lincs_ks = allScores

  load(paste0(baseDir, 'xsum/allScores_cmap.RData'))
  allScores$COMPOUND = toupper(allScores$compound)
  rownames(allScores) = allScores$compound
  scores$cmap_xsum = allScores

  load(paste0(baseDir, 'xsum/allScores_lincs.RData'))
  allScores = AppendPertIName(allScores)
  allScores$PERT_INAME = toupper(allScores$pert_iname)
  rownames(allScores) = allScores$pert_id
  scores$lincs_xsum = allScores

  return(scores)
}

GetTopDrugsPercentile = function(scores, perc1=0.03, perc2=NULL, scoreName='top_cell_score'){
  out = NULL
  s = scores$cmap_xsum[,scoreName]
  idx = which(s > quantile(s, probs=(1-perc1), na.rm=TRUE))
  out$cmap_xsum = scores$cmap_xsum$compound[idx]

  s = scores$cmap_ks[,scoreName]
  idx = which(s > quantile(s, probs=(1-perc1), na.rm=TRUE))
  out$cmap_ks = scores$cmap_ks$compound[idx]

  s = scores$lincs_xsum[,scoreName]
  idx = which(s > quantile(s, probs=(1-perc1), na.rm=TRUE))
  out$lincs_xsum = scores$lincs_xsum$pert_id[idx]

  s = scores$lincs_ks[,scoreName]
  idx = which(s > quantile(s, probs=(1-perc1), na.rm=TRUE))
  out$lincs_ks = scores$lincs_ks$pert_id[idx]

  if(!is.null(perc2)){
    out2 = NULL

    lincs_ks = setdiff(out$lincs_ks, out$lincs_xsum)
    s = scores$lincs_ks[lincs_ks,scoreName]
    idx = which(s > quantile(s, probs=(1-perc2), na.rm=TRUE))
    print(sprintf('keeping %d out of %d from lincs_ks', length(idx), length(lincs_ks)))
    out2$lincs_ks = union(lincs_ks[idx], intersect(out$lincs_ks, out$lincs_xsum))

    lincs_xsum = setdiff(out$lincs_xsum, out$lincs_ks)
    s = scores$lincs_xsum[lincs_xsum,scoreName]
    idx = which(s > quantile(s, probs=(1-perc2), na.rm=TRUE))
    print(sprintf('keeping %d out of %d from lincs_xsum', length(idx), length(lincs_xsum)))
    out2$lincs_xsum = union(lincs_xsum[idx], intersect(out$lincs_xsum, out$lincs_ks))

    cmap_ks = setdiff(out$cmap_ks, out$cmap_xsum)
    s = scores$cmap_ks[cmap_ks,scoreName]
    idx = which(s > quantile(s, probs=(1-perc2), na.rm=TRUE))
    print(sprintf('keeping %d out of %d from cmap_ks', length(idx), length(cmap_ks)))
    out2$cmap_ks = union(cmap_ks[idx], intersect(out$cmap_ks, out$cmap_xsum))

    cmap_xsum = setdiff(out$cmap_xsum, out$cmap_ks)
    s = scores$cmap_xsum[cmap_xsum,scoreName]
    idx = which(s > quantile(s, probs=(1-perc2), na.rm=TRUE))
    print(sprintf('keeping %d out of %d from cmap_xsum', length(idx), length(cmap_xsum)))
    out2$cmap_xsum = union(cmap_xsum[idx], intersect(out$cmap_ks, out$cmap_xsum))
    out = out2
  }

  return(out)
}


Compare2ScoresOverlap = function(scores, percentile=0.1, scoreName='top_cell_score'){
  out1 = GetTopDrugsPercentile(scores=scores, percentile=percentile, scoreName=scoreName)
  cmap_xsum = out1$cmap_xsum
  cmap_ks = out1$cmap_ks
  lincs_xsum = out1$lincs_xsum
  lincs_ks = out1$lincs_ks

  out = NULL
  out$cmap = NULL
  out$lincs = NULL
  print(sprintf('percentile: %f', percentile))

  out$cmap_p = FisherExact(cmap_xsum, cmap_ks, n_universe = length(na.omit(scores$cmap_ks[[scoreName]])))$p
  out$lincs_p = FisherExact(lincs_xsum, lincs_ks, n_universe = length(na.omit(scores$lincs_ks[[scoreName]])))$p

  out$cmap$both = SubsetScores(cmap_ks, cmap_xsum, whichScore='both', scores)
  out$lincs$both = SubsetScores(lincs_ks, lincs_xsum, whichScore='both', scores)

  out$cmap$ks = SubsetScores(cmap_ks, cmap_xsum, whichScore='ks', scores)
  out$cmap$xsum = SubsetScores(cmap_ks, cmap_xsum, whichScore='xsum', scores)

  out$lincs$ks = SubsetScores(lincs_ks, lincs_xsum, whichScore='ks', scores)
  out$lincs$xsum = SubsetScores(lincs_ks, lincs_xsum, whichScore='xsum', scores)

  return(out)
}

SubsetScores = function(ks_drugs, xsum_drugs, whichScore, scores){
  if(any(grepl('BRD', ks_drugs))){
    database = 'lincs'
  }else{
    database = 'cmap'
  }

  if(whichScore == 'ks'){
    drugs = setdiff(ks_drugs, xsum_drugs)
    scoreSheet = paste0(database, '_ks')
    topScores = SubsetAllScoresColumns(scores[[scoreSheet]], database, extras=FALSE, orderBy='top_cell_score')[drugs,]
    topScores = ChangeColumnName(topScores, from = 'top_cell_score', to=paste0('top_cell_score.', whichScore))
  }else if(whichScore == 'xsum'){
    drugs = setdiff(xsum_drugs, ks_drugs)
    scoreSheet = paste0(database, '_xsum')
    topScores = SubsetAllScoresColumns(scores[[scoreSheet]], database, extras=FALSE, orderBy='top_cell_score')[drugs,]
    topScores = ChangeColumnName(topScores, from = 'top_cell_score', to=paste0('top_cell_score.', whichScore))
  }else if(whichScore == 'both'){
    drugs = intersect(ks_drugs, xsum_drugs )
    scoreSheet1 = paste0(database, '_ks')
    topScores1 = SubsetAllScoresColumns(scores[[scoreSheet1]], database, extras=FALSE)[drugs,]
    scoreSheet2 = paste0(database, '_xsum')
    topScores2 = SubsetAllScoresColumns(scores[[scoreSheet2]], database, extras=FALSE)[drugs,]
    topScores = merge(topScores1, topScores2, by=setdiff(intersect(names(topScores1), names(topScores2)), 'top_cell_score'), suffixes = c('.ks', '.xsum'))
  }else{
    stop('unexpected value for whichScore')
  }
  return(topScores)
}

GetTopDrugs = function(scores, write=FALSE){
  out = Compare2ScoresOverlap(scores, percentile=0.03)
  final = NULL

  final$cmap = merge(out$cmap$both, out$cmap$ks[1:3,], all=TRUE)
  final$cmap = merge(final$cmap, out$cmap$xsum[1:3,], all=TRUE)

  final$lincs = merge(out$lincs$both, out$lincs$ks[1:27,], all=TRUE)
  final$lincs = merge(final$lincs, out$lincs$xsum[1:27,], all=TRUE)

  rownames(final$cmap) = final$cmap$compound
  rownames(final$lincs) = final$lincs$unique_id

  if(write){
    WriteXLS('final', ExcelFileName = 'CF_drugs_150.xlsx',
             row.names = TRUE, col.names=TRUE,
             AdjWidth = FALSE, BoldHeaderRow=TRUE,
             FreezeRow=1, FreezeCol=1)
  }

  return(final)
}

GetDbp = function(database){
  dbParams = NULL
  dbParams$cmap = list(mergeCol='compound', compoundName='compound')
  dbParams$lincs = list(mergeCol='unique_id', #mergeCol=c('pert_id', 'pert_iname', 'nPert', 'cellContributors'),
                         compoundName='pert_iname')
  return(dbParams[[database]])
}

GetPertIds = function(){
  load('data/lincs/pertInfo/pertIds13k.RData')
  return(pertIds)
}

GetVennCategory = function(row){

  if(row$selected_multi_cell_score.ks && row$selected_multi_cell_score.xsum){
    category = 'multi_cell_both'
    score = row$multi_cell_score.ks + row$multi_cell_score.xsum
  }else if(row$selected_top_cell_score.ks && row$selected_top_cell_score.xsum){
    category = 'single_cell_both'
    score = row$top_cell_score.ks + row$top_cell_score.xsum
  }else if(row$selected_multi_cell_score.ks){
    category = 'multi_cell_ks'
    score = row$multi_cell_score.ks
  }else if(row$selected_multi_cell_score.xsum){
    category = 'multi_cell_xsum'
    score = row$multi_cell_score.xsum
  }else if(row$selected_top_cell_score.ks){
    category = 'single_cell_ks'
    score = row$top_cell_score.ks
  }else if(row$selected_top_cell_score.xsum){
    category = 'single_cell_xsum'
    score = row$top_cell_score.xsum
  }else{
    stop('did not find appropriate category')
  }
  return(list(category=category, score=score))
}



