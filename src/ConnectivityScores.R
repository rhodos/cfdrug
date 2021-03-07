
library(foreach)

LINCS_CS_parallel = function(querySig, signatureBlockFiles, sigFilters=NULL, nPerm=50,
                                  show_progress=TRUE, annotFilter=TRUE, debug=FALSE, csMetric='ks'){

  if(debug){signatureBlockFiles = signatureBlockFiles[1]}

  
  ### compute background distribution of scores
  cat('simulating null model...\n')
  RL = loadLincsBlock(signatureBlockFiles[1], sigFilters=sigFilters, debug=debug)$RL
  if(csMetric == 'ks'){
    bgDist = est_emp_Cs(querySig, nPerm, RANKED_LISTS=RL, show_progress=show_progress)
  }else if(csMetric == 'xsum'){
    load(DataDir('lincs_mean_zscores.RData'))
    drugData = cbind(zMean, zMean)
    rownames(drugData) = RL[,1]
    bgDist = est_emp_XSum(querySig, nPerm, drugData)
  }
  
  cat('done!\n')
  
  ###iterate over blocks of 10k signatures
  warning('change this to dorng, and check')
  res_i = foreach(blockFile=signatureBlockFiles) %dopar% {
    CS_LINCS(querySig, blockFile = blockFile, bgDist=bgDist, csMetric=csMetric, debug=debug, 
             sigFilters=sigFilters, nPerm=nPerm, show_progress=show_progress)    
  }
  
  dfList = lapply(res_i, function(x) assemble_df(x))
  combRes = ldply(dfList, rbind)
  resDat = combRes[order(combRes$score,decreasing=FALSE),]
  
  if(annotFilter){
    load(DataDir('lincs_compound_metadata.RData'))
    resDat = resDat[resDat$pert_id %in% lincsAnnot$pert_id,]
  }
  
  resDat_pert_type = split(resDat, f = resDat$pert_type)
  for(type_i in 1:length(resDat_pert_type)){
    resDat_pert_type[[type_i]]$adjPvalue = p.adjust(p = resDat_pert_type[[type_i]]$pvalue, method = "BH")
  }
  
  result = lapply(resDat_pert_type, function(x) collapsePerturbagens(x)) 
  return(result)
}

CS_LINCS = function(querySig, blockFile, bgDist=NULL, csMetric='ks', sigFilters=NULL, debug = FALSE, 
                     nPerm=NULL, show_progress=TRUE){
  
  out = loadLincsBlock(blockFile, sigFilters=sigFilters, debug=debug)

  if(ncol(out$sigMat) == 0){  #check for matrix width (ie some non -666 sigs)
    
    res = vector("list")
    
  }else if(ncol(out$sigMat) > 0){
        
    cat('computing connectivity scores\n')
    
    ns = ncol(out$RL)
    
    if(show_progress){
      pb = txtProgressBar(min=1,max=ns,style=3)
    }
    
    if(csMetric == 'xsum'){
      load(DataDir('lincs_mean_zscores.RData'))
      drugDataVec = rev(zMean)
      mu = mean(bgDist)
      stdev = sd(bgDist)
    }

    
    CS = rep(NA,ns)
    Pvals = rep(NA,ns)

    for (i in 1:ns){
      if(csMetric == 'ks'){
        CS[i]<-cMap_CS(out$RL[,i],querySig)    
        Pvals[i] = pnormmix(CS[i], bgDist)
      }else if(csMetric == 'xsum'){
        names(drugDataVec) = out$RL[,i]
        CS[i] = XSum_CS(drugDataVec, querySig)
      }
      if(show_progress){
        setTxtProgressBar(pb, i)
      }
    }
    
    if(show_progress){
      Sys.sleep(1)
      close(pb)
    }
    
    if(csMetric == 'ks'){
      NCS<-rep(NA,length(CS))
      NCS[CS>=0] = CS[CS>=0]/max(bgDist$mu)
      NCS[CS<0] = -CS[CS<0]/min(bgDist$mu)
    }else if(csMetric == 'xsum'){
      NCS = (CS - mu) / stdev
      Pvals = 2*(1-pnorm(abs(NCS)))
    }

    names(CS) = colnames(out$RL)
    names(Pvals) = colnames(out$RL)
    names(NCS) = colnames(out$RL)
    
    res<-list(CS=CS, Pval=Pvals, NCS=NCS, annotations=out$annot)
  }
  
  return(res)
}

collapsePerturbagens = function(filteredRes){
  
  # Summarize results for same perturbagen
  require(plyr)  
  filteredRes = Factor2Char(filteredRes[!is.na(filteredRes$pert_id),])
  meanScoreTable = ddply(filteredRes, .(pert_id), 
                          plyr::summarize, pert_iname=paste(unique(pert_iname), collapse="|"), 
                          nPert = length(score), 
                          cellContributors = paste(sort(cell), collapse = "|"),
                          signedLogP = -sign(mean(normScore))*log10(pcomb(pvalue)),
                          score = mean(score),
                          normScore = mean(normScore),
                          pvalue = pcomb(pvalue), 
                          .parallel = TRUE)
  
  meanScoreTable = meanScoreTable[order(meanScoreTable$pvalue,decreasing = FALSE),]
  meanScoreTable$adjPvalue = p.adjust(p = meanScoreTable$pvalue, method = "BH")
  
  return(meanScoreTable)
}

### Stouffers method for combining p-values
erf = function(x) 2 * pnorm(2 * x/ sqrt(2)) - 1
erfinv = function(x) qnorm( (x+1)/2 ) / sqrt(2)
pcomb = function(p) (1-erf(sum(sqrt(2) * erfinv(1-2*p))/sqrt(2*length(p))))/2

pertMfc2RegularId = function(x){
  if(is.null(x)){
    out = NA
  }else{
    out = paste(unlist(strsplit(x,split='-'))[1:2], collapse='-')
  }
  return(out)
}

loadLincsBlock = function(blockFile, sigFilters, debug=FALSE){

  cat(sprintf('loading file %s\n', blockFile))
  load(file = blockFile)
  
  #align sig info with matrix ie we have some sig info, with -666 annotations and no data
  allSigs = allSigs[match(colnames(sigMat), names(allSigs))] 
  
  if(ncol(sigMat) > 0){
    
    if(sigFilters$onlyGold){
      gold_status = sapply(allSigs, function(x) x$is_gold)
      allSigs = allSigs[which(gold_status)]
      sigMat =  sigMat[,which(gold_status)]
    }
    
    if(sigFilters$onlyChemPert){
      pert_type = sapply(allSigs, function(x) x$pert_type)
      allSigs = allSigs[which(pert_type == 'trt_cp')]
      sigMat = sigMat[,which(pert_type == 'trt_cp')]
    }
    
    nFilter = sigFilters$debugNumber
    if(debug && ncol(sigMat) > nFilter){
      allSigs = allSigs[1:nFilter]
      sigMat = sigMat[,1:nFilter]
    }
        
    if(ncol(sigMat) > 0){
          
      geneOrder = rownames(sigMat)
      RANKED_LISTS = sapply(1:ncol(sigMat), function(x) 
        rownames(sigMat)[order(sigMat[,x], decreasing = TRUE)])
      colnames(RANKED_LISTS) = colnames(sigMat)
      
      annot = data.frame(t(sapply(allSigs, function(x) c(pertMfc2RegularId(
        x$pert_mfc_id), x$pert_iname, x$pert_idose, x$pert_time, x$cell_id, x$pert_type, x$distil_ss))))
      colnames(annot) = c("pert_id", "pert_iname", "dose", "time", "cell", "pert_type", "distil_ss")

    }else{
      RANKED_LISTS = NULL
      annot = NULL
    }
  }else{
    # yes, this is intentional. There are two ways that sigMat can be empty-
    # either there is no data to start with, or there is no data after filtering
    RANKED_LISTS = NULL
    annot = NULL
  }
    
  return(list(RL=RANKED_LISTS, allSigs=allSigs, sigMat=sigMat, annot=annot))
}



assemble_df = function(resList){
  resDat = data.frame(resList$annotations, signature_id = names(resList$CS), 
                       score = resList$CS, normScore = resList$NCS, pvalue = resList$Pval, row.names = NULL) 
}



############################################################################################################
### X sum, as per Cheng, Agarwal et al 2014
############################################################################################################

XSum_allCmap = function(querySig, nPerm, numDrugGenes=500, data='FC', takeAbs=FALSE, debug=FALSE, debugN=2, database='cmap'){

  if(database=='cmap'){
    warning('should remove other two cell types from FC matrices')
    if(data == 'FC'){
#       load('data/cmap/FC1309.RData')
#       drugData = FC1309
      
      load(DataDir('cmap/cell_specific_FC_cmap.RData'))
      drugData = meanFC
    }else if(data == 'qFC'){
      load(DataDir('cmap/cell_specific_qFC_cmap.RData'))
      drugData = FC
    }
    drugDataVec = drugData[,1]
  }else if(database=='lincs'){
    load(DataDir('lincs_mean_zscores.RData'))
    if(debug && debugN<=5){
      load('/la-forge/rhodos/projects/drug_repurposing/data/lincs/PRL_entrez_lincs_cell_specific_debug.RData')
    }else{
      load('/la-forge/rhodos/projects/drug_repurposing/data/lincs/PRL_entrez_lincs_cell_specific_plusCF.RData')
    }
    drugData = PRL_entrez
    drugDataVec = rev(zMean) # reversing so that the first entry of drugDataVec is the largest z-score, 
                              # since PRL_entrez rank 1 corresponds to most up-regulated gene
    names(drugDataVec) = as.character(PRL_entrez[,1])
  }else{
    stop('unexpected database type')
  }

  if(debug){
     drugData = drugData[,1:debugN]
  }
  
  cat('estimating empirical distribution for XSum...\n')
  print(system.time(CS_EMP <- est_emp_XSum(querySig, nPerm, drugDataVec=drugDataVec, numDrugGenes=numDrugGenes, takeAbs=takeAbs)))
  cat('...done!\n')
  mu = mean(CS_EMP)
  stdev = sd(CS_EMP)
  
  cat('computing connectivity scores...\n')
  warning('change this to dorng, and check')
  CS = foreach(i=1:ncol(drugData)) %dopar%{
    if(database == 'cmap'){
      drugDataVec = drugData[,i]
    }else{
      names(drugDataVec) = as.character(PRL_entrez[,i])
    }
    XSum_CS(drugDataVec, querySig, numDrugGenes=numDrugGenes, takeAbs=takeAbs)
  }
  CS = unlist(CS)
  names(CS) = colnames(drugData)
  cat('...done!\n')
  
  cat('computing pvalues...\n')
  NCS = (CS - mu) / stdev
  Pvals = 2*(1-pnorm(abs(NCS)))
  cat('...done!\n')
  
  return(list(CS=CS,Pval=Pvals,adjP=stats::p.adjust(Pvals,method='fdr'),NCS=NCS))#, mu=mu, stdev=stdev, CS_EMP=CS_EMP))
}

# drugDataVec should be a named vector with names corresponding to entrez gene
# IDs and values either logFC or z-score
# 
# querySig should be a list with at least one element named UP (containing a
# vector of entrez ids) and possibly another list called DOWN, also containing a
# vector of entrez ids
#
# not sure if takeAbs does anything useful.  I made this up..
XSum_CS = function(drugDataVec, querySig, numDrugGenes=500, takeAbs=FALSE){ 

  if(takeAbs){
    drugDataVec = abs(drugDataVec)
  }
  
  minRankUp = length(drugDataVec) - numDrugGenes + 1   
  minRankDown = numDrugGenes
  
  geneRanks = rank(drugDataVec)
  idxUp = which(geneRanks >= minRankUp)
  idxDown = which(geneRanks <= minRankDown)
  changedByCompound = names(drugDataVec)[c(idxUp, idxDown)]
  
  XUP = intersect(querySig$UP, changedByCompound)
  XDOWN = intersect(querySig$DOWN, changedByCompound)

  if(length(XUP) > 0){
    XSum_up = sum(drugDataVec[XUP]) 
  }else{
    XSum_up = 0
  }
  
  if(length(XDOWN) > 0){
    XSum_down = sum(drugDataVec[XDOWN])
  }else{
    XSum_down = 0
  }
  
  if(takeAbs){
    XSum = XSum_up + XSum_down
  }else{
    XSum = XSum_up - XSum_down
  }
  
  return(XSum)
}

est_emp_XSum = function(querySig, nPerm, drugDataVec, numDrugGenes=500, takeAbs=FALSE){
  set.seed(42)
  names = names(drugDataVec)
  ng = length(drugDataVec)

  warning('change this to dorng, and check')
  EMP_CS = foreach(i = 1:nPerm) %dopar% {
   names(drugDataVec) = names[sample(1:ng, ng)]
   XSum_CS(drugDataVec, querySig, numDrugGenes=numDrugGenes, takeAbs=takeAbs)
  }

  return(unlist(EMP_CS))
}


qES = function(RANKEDLIST,REGULON,display=FALSE,returnRS=FALSE){
  
  REGULON<-intersect(as.character(REGULON),RANKEDLIST)
  
  HITS<-is.element(RANKEDLIST,REGULON)+0
  
  hitCases<-cumsum(HITS)
  missCases<-cumsum(1-HITS)
  
  N<-length(RANKEDLIST)
  NR<-length(REGULON)
  
  Phit<-hitCases/NR
  Pmiss<-missCases/(N-NR)
  
  m<-max(abs(Phit-Pmiss))
  t<-which(abs(Phit-Pmiss)==m)
  
  if (length(t)>1){t<-t[1]}
  peak<-t
  ES<-Phit[t]-Pmiss[t]
  RS<-Phit-Pmiss
  
  if (display){
    if (ES>=0){c<-"red"}else{c<-"green"}
    plot(0:N,c(0,Phit-Pmiss),col=c,type="l",xlim=c(0,N),ylim=c(-(abs(ES)+0.5*(abs(ES))),abs(ES)+0.5*(abs(ES))),xaxs="i",bty="l",axes=FALSE,
         xlab="Gene Rank Position",ylab="Running Sum")
    par(new=TRUE)
    plot(0:N,rep(0,N+1),col='gray',type="l",new=FALSE,xlab="",ylab="",ylim=c(-(abs(ES)+0.5*(abs(ES))),abs(ES)+0.5*(abs(ES))))
    axis(side=2)
    
  }
  
  if (returnRS){
    POSITIONS<-which(HITS==1)
    names(POSITIONS)<-RANKEDLIST[which(HITS==1)]
    
    POSITIONS<-POSITIONS[order(names(POSITIONS))]
    names(POSITIONS)<-names(POSITIONS)[order(names(POSITIONS))]
    
    return(list(ES=ES,RS=RS,POSITIONS=POSITIONS,PEAK=t))
  } else {return(ES)}
}

pnormmix = function(x,mixture) {
  lambda = mixture$lambda
  k = length(lambda)
  pnorm.from.mix = function(x,component) {
    if (x>=0){
      lambda[component]*pnorm(-x,mean=-mixture$mu[component],sd=mixture$sigma[component],lower.tail=TRUE)
    }else {
      lambda[component]*pnorm(x,mean=mixture$mu[component],sd=mixture$sigma[component],lower.tail=TRUE)
    }
  }
  pnorms = sapply(1:k, function(i) pnorm.from.mix(x, i))
  return(sum(pnorms))
}

computeCsOnCmap = function(querySig, nPerm, csMetric='ks', numDrugGenes=500, debug=FALSE, debugN=NULL, database='cmap'){
  if(csMetric == 'ks'){
    csList = CS_entrez(querySig, nPerm, debug=debug, debugN=debugN, database=database)
  }else if(csMetric == 'xsum'){
    csList = XSum_allCmap(querySig, nPerm, numDrugGenes=numDrugGenes, data='qFC', 
                        takeAbs=FALSE, debug=debug, debugN=debugN, database=database)
  }else{
    stop('unexpected csMetric, not computing anything')
  }
  
  out = data.frame(compound = names(csList$CS), score = csList$CS,
             normScore = csList$NCS, pvalue = csList$Pval, 
             adjP = csList$adjP)
  
  return(out)
}

CS_entrez = function(querySig, nPerm, debug=FALSE, debugN=NULL, csMetric='ks', numDrugGenes=500, database='cmap'){

  if(database=='cmap'){
    if(debug){
      load(DataDir('PRL_entrez_cmap_cell_specific_debug.RData'))
    }else{
      cat('loading CMAP PRL file...\n')
      load(DataDir('PRL_entrez_cmap_cell_specific.RData'))
      cat('...done!\n')
    }
  }else if(database=='lincs'){
    if(debug){
      load('/la-forge/rhodos/projects/drug_repurposing/data/lincs/PRL_entrez_lincs_cell_specific_debug.RData')
    }else{
      cat('loading LINCS PRL file...\n')
      print(system.time(load('/la-forge/rhodos/projects/drug_repurposing/data/lincs/PRL_entrez_lincs_cell_specific_plusCF.RData')))
    } 
  }

  if(debug){
    PRL_entrez = PRL_entrez[,1:debugN]
  }
  
  RANKED_LISTS = PRL_entrez
    
  cat('simulating null model\n')
  mixmdl<-est_emp_Cs(querySig, nPerm, RANKED_LISTS)
  cat('done!\n')
  
  cat('computing connectivity scores\n')
  
  ns<-ncol(RANKED_LISTS)
  
  CS<-rep(NA,ns)
  Pvals<-rep(NA,ns)
  
  warning('change this to dorng, and check')
  CS = foreach(i=1:ns) %dopar% {
    cMap_CS(RANKED_LISTS[,i],querySig) 
  }
  CS = unlist(CS)
  
  warning('change this to dorng, and check')
  Pvals = foreach(i=1:ns) %dopar%{
    pnormmix(CS[i], mixmdl)
  }
  Pvals = unlist(Pvals)
  
  names(CS)<-colnames(RANKED_LISTS)
  names(Pvals)<-colnames(RANKED_LISTS)
    
  NCS<-rep(NA,length(CS))
  NCS[CS>=0]<-CS[CS>=0]/max(mixmdl$mu)
  NCS[CS<0]<--CS[CS<0]/min(mixmdl$mu)
  
  names(NCS)<-names(CS)
  
  res<-list(CS=CS,Pval=Pvals,adjP=stats::p.adjust(Pvals,method='fdr'),NCS=NCS)
  
  return(res)
}

# This differs from CS_entrez original function in that you input your own drug signatures
CS2 = function(querySig, drugSigs, nPerm){

  # Reformat drugSigs into RANKED_LISTS matrix
  RANKED_LISTS = array(data='', dim=c(length(drugSigs[[1]]), length(drugSigs)))
  for(i in 1:length(drugSigs)){
    RANKED_LISTS[,i] = names(sort(-drugSigs[[i]]))
  }
  
  mixmdl<-est_emp_Cs(querySig, nPerm, RANKED_LISTS)

  ns<-ncol(RANKED_LISTS)
  
  CS<-rep(NA,ns)
  Pvals<-rep(NA,ns)

  CS = foreach(i=1:ns) %dopar% {
    cMap_CS(RANKED_LISTS[,i],querySig) 
  }
  CS = unlist(CS)
  
  Pvals = foreach(i=1:ns) %dopar%{
    pnormmix(CS[i], mixmdl)
  }
  Pvals = unlist(Pvals)
  
  names(CS)<-colnames(RANKED_LISTS)
  names(Pvals)<-colnames(RANKED_LISTS)
  
  NCS<-rep(NA,length(CS))
  NCS[CS>=0]<-CS[CS>=0]/max(mixmdl$mu)
  NCS[CS<0]<--CS[CS<0]/min(mixmdl$mu)
  
  names(NCS)<-names(CS)
  
  res<-list(CS=CS,Pval=Pvals,adjP=stats::p.adjust(Pvals,method='fdr'),NCS=NCS)
  
  return(res)
}


cMap_CS = function(ranked_list, querySig, returnRS=FALSE){
  
  ESUP = qES(ranked_list,querySig$UP,display=FALSE,returnRS=returnRS)
  twoSided = !is.null(querySig$DOWN)
  
  if(twoSided){
    ESDOWN = qES(ranked_list,querySig$DOWN,display=FALSE,returnRS=returnRS)
  }
  
  if (returnRS){
    RSUP<-ESUP$RS
    ESUP<-ESUP$ES
    
    if(twoSided){
      RSDOWN<-ESDOWN$RS
      ESDOWN<-ESDOWN$ES
    }
  }
  
  if(!twoSided){
    TES = ESUP
  }else{
    TES<-(ESUP-ESDOWN)/2
  }
  
  if(returnRS){
    out = list(TES=TES, ESUP=ESUP, ESDOWN=ESDOWN, RSUP=RSUP, RSDOWN=RSDOWN)
  }else{
    out = TES
  }
  
  return(out)
}


est_emp_Cs = function(querySig,nt,RANKED_LISTS,show_progress=TRUE){
  set.seed(42)
  require(mixtools)
  
  EMP_CS = rep(NA,nt)
  
  ng = nrow(RANKED_LISTS)
  
  if (show_progress){
    pb = txtProgressBar(min=1,max=nt,style=3)
  }
  
  warning('change this to dorng, and check')
  EMP_CS = foreach(i = 1:nt) %dopar%{
    cMap_CS(RANKED_LISTS[sample(1:ng,ng),1],querySig)
  }
  EMP_CS = unlist(EMP_CS)
  
  if(show_progress){
    Sys.sleep(1)
    close(pb)
  }
  
  if(is.null(querySig$DOWN)){
    mixmdl = normalmixEM(EMP_CS,k=2,verb=FALSE) # one-sided
  }else{
    mixmdl = normalmixEM(EMP_CS,k=3,verb=FALSE) # two-sided
  }
  
  return(mixmdl)
}
