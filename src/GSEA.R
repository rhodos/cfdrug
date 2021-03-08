library(HTSanalyzeR)
library(GO.db)
library(KEGG.db)
source('R/src/MyHTSAnalyzeR.R')

RunGsea = function(profile, dz_genes_up, dz_genes_down=NULL, minGeneSetSize=15, nPermute=1000, seed=123,
                   doGSOA=FALSE, doGSEA=TRUE, ListGSC = GetListGSC(), appendTerms=TRUE, preprocess=TRUE){

  set.seed(seed)

  gsca = new('GSCA', listOfGeneSetCollections=ListGSC, geneList=profile, hits=as.character(c(dz_genes_up, dz_genes_down)))

  if(preprocess){
    gsca = preprocess(gsca, species='Hs', initialIDs='Entrez.gene', keepMultipleMappings=TRUE,
                      duplicateRemoverMethod='max', orderAbsValue=FALSE)
  }

  gsca = analyze(gsca, para=list(pValueCutoff=0.05, pAdjustMethod='BH', nPermutations=nPermute,
                                 minGeneSetSize=minGeneSetSize, exponent=1), doGSOA=doGSOA, doGSEA=doGSEA)

  if(appendTerms){
    gsca = AppendGSTerms(gsca, go=names(ListGSC)[grepl('GO_', names(ListGSC))],
                         kegg=ifelse('PW_KEGG' %in% names(ListGSC), 'PW_KEGG', NA))
  }

  HTSanalyzeR::summarize(gsca)
  return(gsca)
}

AppendGSTerms = function(gsca, go = c('GO_BP', 'GO_MF', 'GO_CC'), kegg = c('PW_KEGG')){
  return(appendGSTerms(gsca, goGSCs=go, keggGSCs=kegg))
}

GetGSEAResults = function(gsca, colname='Adjusted.Pvalue', thresh=0.05, analysis='GSEA', merge=TRUE){
  if(analysis == 'GSOA'){
    x = gsca@result$HyperGeo.results
  }else if(analysis == 'GSEA'){
    x = gsca@result$GSEA.results
  }else{
    stop('unexpected value for analysis')
  }

  out = list()
  for(n in names(x)){
    y = x[[n]][,colname]
    out[[n]] = as.data.frame(x[[n]][y <= thresh,,drop=FALSE])
  }

  if(merge){
    geneSets = names(which(sapply(out, nrow) > 0))
    out = lapply(geneSets, function(name){out[[name]]$Gene.Set=name; return(out[[name]])})
    out = Reduce('rbind', out)
  }
  return(out)
}

GetListGSC = function(){
  GO_MF = GOGeneSets(species='Hs', ontologies=c('MF'))
  GO_BP = GOGeneSets(species='Hs', ontologies=c('BP'))
  GO_CC = GOGeneSets(species='Hs', ontologies=c('CC'))
  PW_KEGG = KeggGeneSets(species='Hs')
  return(list(GO_MF=GO_MF, GO_BP=GO_BP, GO_CC=GO_CC, PW_KEGG=PW_KEGG))
}

# Takes as input a data.frame such as the output of GetGSEAResults
FilterRes = function(res, fdrThresh=0.1){
  library(dplyr)
  res %<>% dplyr::mutate(Term=rownames(res)) %>% MoveColumn('Term', 1) %>% subset(FDR <= fdrThresh)
  return(res[order(-res$Observed.score),])
}

FilterGSOARes = function(res){
  library(dplyr)
  if(is.null(res)){
    out = NULL
  }else{
    res$FoldEnrichment = res[,'Observed Hits'] / res[,'Expected Hits']
    res %<>% dplyr::mutate(Term=rownames(res)) %>% MoveColumn('Term', 1)
    out = res[order(-res$FoldEnrichment),]
  }
  return(out)
}

GetTFNames = function(chea_names){
  sapply(chea_names, function(x) strsplit(x, split='-')[[1]][1]) %>% as.character
}

# Append top differentially expressed genes from enrichmed pathways to column of enrichment results
AppendTopGeneColumn = function(res, gmt, meta){
  topGenes = vector(mode='character', length=nrow(res))
  for(i in 1:nrow(res)){
    term = res$Term[i]
    library = res$Gene.Set[i]
    score_sign = sign(res$Observed.score[i])
    M = meta[meta$GeneID %in% gmt[[library]][[term]], ]
    if(score_sign > 0){
      topM = subset(M, metaEffect > 0.1) %>% arrange(-metaEffect)
    }else{
      topM = subset(M, metaEffect < -0.1) %>% arrange(metaEffect)
    }
    if(nrow(topM)>5){topM = topM[1:5,]}
    topM %<>% mutate(gene_with_effect=sprintf('%s%s%s(%0.2f)',
                  Symbol, ifelse(metaEffectP.corrected<0.05, '*', ''),
                  ifelse(!is.na(pfp.two.sided.og2b)&pfp.two.sided.og2b<0.05, '#', ''), 2^metaEffect))
    topGenes[i] = paste(topM$gene_with_effect, collapse=', ')
  }
  res$TopGenes = topGenes
  return(res)
}
