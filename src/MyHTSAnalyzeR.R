
# All code that I changed is commented out. These are functions taken from the
# HTSanalyzeR package, just changing the parallel functionality, as the snow
# package wasn't utilizing the cores very well

paraCheck = HTSanalyzeR:::paraCheck

collectionGsea = function(collectionOfGeneSets, geneList, exponent = 1, nPermutations = 1000,
          minGeneSetSize = 15, verbose = TRUE){
  paraCheck("gsc", collectionOfGeneSets)
  paraCheck("genelist", geneList)
  paraCheck("exponent", exponent)
  paraCheck("minGeneSetSize", minGeneSetSize)
  geneList.names <- names(geneList)
  paraCheck("nPermutations", nPermutations)
  nGeneSets <- length(collectionOfGeneSets)
  tagGeneSets <- rep(FALSE, nGeneSets)
  tagGeneSets[which(unlist(lapply(collectionOfGeneSets, length)) < length(geneList))] <- TRUE
  tagGeneSets[which(unlist(lapply(lapply(collectionOfGeneSets, intersect, y = geneList.names), length)) < minGeneSetSize)] <- FALSE
  n.tagGeneSets <- sum(tagGeneSets)
  if (n.tagGeneSets == 0)
    warning('There are no gene sets in your collection that pass the cutoffs on size')
  if (n.tagGeneSets > 0) {
    scoresperm <- matrix(rep(0, (nPermutations * n.tagGeneSets)), nrow = n.tagGeneSets)
    rownames(scoresperm) <- names(collectionOfGeneSets)[which(tagGeneSets)]
    scoresObserved <- rep(0, n.tagGeneSets)
    names(scoresObserved) <- names(collectionOfGeneSets)[which(tagGeneSets)]
    perm.gL <- sapply(1:nPermutations, function(n) names(geneList)[sample(1:length(geneList), length(geneList), replace = FALSE)])
    perm.gL <- cbind(names(geneList), perm.gL)

    if(!is.null(getOption('cores'))){     #if (is(getOption("cluster"), "cluster") && "package:snow" %in% search()) {
      scores <- gseaScoresBatchParallel(geneList, geneNames.perm = perm.gL,
                                        collectionOfGeneSets = collectionOfGeneSets[which(tagGeneSets)],
                                        exponent = exponent, nPermutations = nPermutations)
      sapply(1:n.tagGeneSets, function(i) {
        scoresperm[i, ] <<- unlist(scores[[i]]$scoresperm)
        scoresObserved[i] <<- unlist(scores[[i]]$scoresObserved)
      })
    }else {
      if (verbose)
        pb <- txtProgressBar(style = 3)
      for (i in 1:n.tagGeneSets) {
        scores <- gseaScoresBatch(geneList, geneNames.perm = perm.gL,
                                  geneSet = collectionOfGeneSets[[which(tagGeneSets)[i]]],
                                  exponent = exponent, nPermutations = nPermutations)
        scoresObserved[i] <- scores$scoresObserved
        scoresperm[i, ] <- scores$scoresperm
        if (verbose)
          setTxtProgressBar(pb, i/n.tagGeneSets)
      }
      if (verbose)
        close(pb)
    }
  }
  else {
    scoresObserved <- NULL
    scoresperm <- NULL
  }
  return(list(Observed.scores = scoresObserved, Permutation.scores = scoresperm))
}

gseaScoresBatchParallel = function(geneList, geneNames.perm, collectionOfGeneSets, exponent=1,
                                   nPermutations = 1000){
  paraCheck("genelist", geneList)
  paraCheck("gsc", collectionOfGeneSets)
  paraCheck("exponent", exponent)
  paraCheck("nPermutations", nPermutations)
  if (!is.matrix(geneNames.perm))
    stop("'geneNames.perm' should be a matrix!\n")
  if (ncol(geneNames.perm) != (nPermutations + 1))
    stop("The No of columns of 'geneNames.perm' should be equal to 'nPermutations'!\n")
  gseaScoresBatchLocal <- function(geneList, geneNames.perm,
                                   geneSet, exponent, nPermutations) {
    geneList.names <- names(geneList)
    geneSet <- intersect(geneList.names, geneSet)
    nh <- length(geneSet)
    N <- length(geneList)
    ES <- rep(0, nPermutations + 1)
    Phit <- matrix(0, nrow = N, ncol = nPermutations + 1)
    Pmiss <- Phit
    runningES <- NULL
    if (nh > N) stop("Gene Set is larger than Gene List")
    hits <- matrix(FALSE, nrow = N, ncol = nPermutations +1)
    hits[which(!is.na(match(geneNames.perm, geneSet)))] <- TRUE
    hits <- matrix(hits, ncol = nPermutations + 1, byrow = FALSE)
    if (sum(hits[, 1]) > 0) {
      junk <- sapply(1:(nPermutations + 1), function(i) Phit[which(hits[,i]), i] <<- abs(geneList[which(hits[, i])])^exponent)
      NR <- colSums(Phit)
      Pmiss[which(!hits)] <- 1/(N - nh)
      Pmiss <- sapply(1:(nPermutations + 1), function(i) cumsum(Pmiss[,i]))
      Phit <- sapply(1:(nPermutations + 1), function(i) cumsum(Phit[,i])/NR[i])
      runningES <- Phit - Pmiss
      ESrange <- sapply(1:(nPermutations + 1), function(i) range(runningES[,i]))
      ES <- sapply(1:(nPermutations + 1), function(i) ESrange[which.max(abs(ESrange[,i])), i])
      if (is.list(ES))
        ES <- unlist(ES)
    }
    ES <- list(scoresObserved = ES[1], scoresperm = ES[2:(nPermutations + 1)])
    return(ES)
  }

  scores = foreach(i=1:length(collectionOfGeneSets)) %dopar% {
    gseaScoresBatchLocal(geneList, geneNames.perm = geneNames.perm,
                        geneSet = as.integer(collectionOfGeneSets[[i]]),
                        exponent = exponent, nPermutations = nPermutations)
  }
  return(scores)
}

FDRcollectionGsea = function(permScores, dataScores){

  if (!is.matrix(permScores))
    stop("'permScores' should be a matrix!\\n")
  if (!is.numeric(dataScores) && !is.integer(dataScores))
    stop("'dataScores' should be an integer or numeric vector!\\n")
  if (is.null(names(dataScores)))
    stop("'dataScores' should be named (by gene set identifier)")
  if (nrow(permScores) != length(dataScores))
    stop(paste("The number of rows of the 'permScores' matrix ",
               "should be the same as the length of the 'dataScores' vector",
               sep = ""))

  ldataScores <- length(dataScores)
  FDRgeneset = rep(0, ldataScores)

  meanNegOverall = permScores[permScores<=0] %>% mean
  meanPosOverall = permScores[permScores>=0] %>% mean

  sapply(1:ldataScores, function(i) {
    neg <- which(permScores[i, ] <= 0)
    pos <- which(permScores[i, ] >= 0)

    # The following two statements are a workaround that I added for cases where the statistics are
    # all one sign (happens with very large gene sets for some reason)

    if(length(neg)==0){
      warning(sprintf('length(neg)=0 for i=%d', i))
      NegAvg = meanNegOverall
    }else{
      NegAvg <- abs(mean(permScores[i, neg]))
    }

    if(length(pos)==0){
      warning(sprintf('length(pos)=0 for i=%d', i))
      PosAvg = meanPosOverall
    }else{
      PosAvg <- abs(mean(permScores[i, pos]))
    }

    permScores[i, neg] <<- permScores[i, neg]/NegAvg
    permScores[i, pos] <<- permScores[i, pos]/PosAvg

    dataScores[i] <<- ifelse((dataScores[i] < 0), (dataScores[i]/NegAvg), (dataScores[i]/PosAvg))
  })
  negtot <- length(which(permScores <= 0))
  postot <- length(which(permScores >= 0))
  sapply(1:ldataScores, function(i) {
    if (is.na(dataScores[i])) {
      FDRgeneset[i] <<- 1
    }
    else if (dataScores[i] < 0) {
      FDRgeneset[i] <<- (sum(permScores <= dataScores[i])/negtot)/(sum(dataScores <= dataScores[i])/sum(dataScores <= 0))
    }
    else {
      FDRgeneset[i] <<- (sum(permScores >= dataScores[i])/postot)/(sum(dataScores >= dataScores[i])/sum(dataScores >= 0))
    }
    FDRgeneset[i] <<- ifelse(FDRgeneset[i] > 1, 1, FDRgeneset[i])
  })
  names(FDRgeneset) <- names(dataScores)
  return(FDRgeneset)
}


assignInNamespace('gseaScoresBatchParallel', gseaScoresBatchParallel, ns='HTSanalyzeR')
assignInNamespace('collectionGsea', collectionGsea, ns='HTSanalyzeR')
assignInNamespace('FDRcollectionGsea', FDRcollectionGsea, ns='HTSanalyzeR')