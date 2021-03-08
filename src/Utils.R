merge = base::merge
library(matrixStats)
library(affy)
library(reshape2)
library(devtools)


#### directory and file handling ###################################################################

PlotDir = function(file='', subdir=NULL){
  plotDir = paste(OutputDir(), 'figures', sep='/')
  if(!file.exists(plotDir)){
    dir.create(plotDir, recursive=T)
  }
  return(paste(plotDir, file, sep='/'))
}

BaseDir = function(){
  return(GetConfig('CODEPATH'))
}

GetConfig = function(key){
  a = system("cat config/config.txt;", intern=T)
  b = unlist(strsplit(a,"="))
  return(b[which(b==key)+1])
}

DateStr = function(){
  return(gsub('-', '_', Sys.Date()))
}

CheckDir = function(){
  dirs = unlist(strsplit(getwd(), "/"))
  stopifnot(identical(dirs[length(dirs)], 'code'))
}

DataDir = function(subdir=''){
  dir = GetConfig('DATAPATH')
  return(paste0(dir, subdir))
}

OutputDir = function(subdir=''){
  dir = GetConfig('OUTPATH')
  return(paste0(dir, subdir))
}

#### workspace helper functions ####################################################################

LsVars = function(envInt=0){
  setdiff(ls(env=sys.frame(envInt)), lsf.str(env=sys.frame(envInt)))
}

LsFcns = function(envInt=0){
  lsf.str(env=sys.frame(envInt))
}

GetVarName = function(var){
  return(deparse(substitute(var)))
}

#### set operations and helper functions ###########################################################

"%ni%" = Negate("%in%")

IsPartition = function(partitionList, fullSet){
  isPartition = TRUE

  for(i in 1:length(partitionList)){
    for(j in IncreasingSequence(i+1, length(partitionList))){
      if(IsIntersection(partitionList[[i]], partitionList[[j]])){
        isPartition = FALSE
      }
    }
  }

  all = partitionList[[1]]
  for(i in 2:length(partitionList)){
    all = union(all, partitionList[[i]])
  }

  if(!is.numeric(fullSet)){
    if(!identical(sort(all), sort(fullSet))){
      isPartition = FALSE
    }
  }else{
    if(!identical(sort(as.numeric(all)), sort(as.numeric(fullSet)))){
      isPartition = FALSE
    }
  }

  return(isPartition)
}

IsIntersection = function(setA, setB){
  isIntersection = FALSE
  if(length(intersect(setA, setB))>0){
    isIntersection = TRUE
  }
  return(isIntersection)
}

# find all occurrences of each element of A in B
MatchAll = function(A, B, collapse='|'){
  out = c()
  for(a in A){
    matches = paste(which(B == a), collapse=collapse)
    out = c(out, matches)
  }
  out[out==''] = NA
  return(out)
}

# get elements that appear at least k times
GetMultiples = function(x, k=2){
  y = unique(x)
  out = c()
  for(a in y){
    count = length(which(x == a))
    if(count >= k){
      out = c(out, a)
    }
  }
  return(out)
}

# input is a list of items that may have some repeats
# output is a list of unique items, and the corresponding counts

SymDiff = function(x, y){
  return(setdiff(union(x, y), intersect(x, y)))
}

#### math/stats stuff ##############################################################################

ZScore = function(x){
  return( (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
}

IsBinary = function(x){
  return(all(unique(na.omit(x)) %in% c(0,1)))
}

InInterval = function(a, interval){
  stopifnot(is.numeric(interval))
  stopifnot(length(interval) == 2)
  stopifnot(is.numeric(a))
  stopifnot(interval[1] <= interval[2])
  return(all(a >= interval[1]) & all(a <= interval[2]))
}

FisherExact = function(A, B, n_universe=length(union(A,B)), n_hypotheses=1, printFlag=T, alternative='greater'){
  aa = length(intersect(A, B))
  bb = length(setdiff(A, B))
  cc = length(setdiff(B, A))
  dd = n_universe - (aa + bb + cc)

  x = matrix(c(aa,bb,cc,dd), nrow=2, ncol=2, byrow=T)
  out = fisher.test(x, alternative=alternative)
  p = out$p.value

  frac_observed = aa / n_universe
  frac_expected = (aa + cc)*(aa + bb) / (n_universe^2)
  fold_enrichment = frac_observed / frac_expected

  if(printFlag){
    print(sprintf('  %d in A, %d in B, %d in overlap, p=%.2e, fold enrichment=%.1f', length(A), length(B), aa, p, fold_enrichment))
  }
  return(list(p=p, adjp=min(p*n_hypotheses, 1), overlap=aa, odds.ratio=out$estimate, fold.enrichment=fold_enrichment))
}

# Warning: The output significance levels are not adjusted.
OverlapListOfLists = function(list1, list2=NULL, bg1=NULL, bg2=NULL, printFlag=TRUE, collapse=NULL){
  if(is.null(list2)){
    list2 = list1
  }
  stopifnot(identical(class(list1), class(list2)))
  m = length(list1)
  n = length(list2)

  # If neither background sets are input, default to the union of all inputs
  if(is.null(bg1) && is.null(bg2)){
    all1 = Reduce('union',list1)
    all2 = Reduce('union', list2)
    all = union(all1, all2)
    bg1 = all
    bg2 = all
  }

  if(!is.list(bg1) || length(bg1) != m){
    bg1 = lapply(1:m, function(x) bg1)
  }

  # If bg1 is input, but not bg2, default to match bg1
  if(is.null(bg2)){
    tmp = bg1[[1]]
    bg2 = lapply(1:n, function(x) tmp)
  }

  if(!is.list(bg2) || length(bg2) != n){
    bg2 = lapply(1:n, function(x) bg2)
  }

  numOverlap = matrix(NA, nrow=m, ncol=n, dimnames=list(names(list1), names(list2)))
  sigOverlap = matrix(0, nrow=m, ncol=n, dimnames=list(names(list1), names(list2)))
  overlapList = list()
  if(!is.null(collapse)){
    overlapMatrix = numOverlap
  }

  for(i in 1:m){
    print(names(list1)[[i]])
    overlapList[[i]] = list()
    for(j in 1:n){
      bg = intersect(bg1[[i]], bg2[[j]])
      A = intersect(bg, list1[[i]])
      B = intersect(bg, list2[[j]])
      overlapList[[i]][[j]] = paste(intersect(A,B), collapse=collapse)
      numOverlap[i,j] = length(intersect(A,B))
      sigOverlap[i,j] = FisherExact(A, B, n_universe=length(bg),
                                     printFlag=printFlag)$adjp
      if(!is.null(collapse)){
        overlapMatrix[i,j] = overlapList[[i]][[j]]
      }
    }
    names(overlapList[[i]]) = names(list2)
  }
  names(overlapList) = names(list1)

  out = list(numOverlap=numOverlap, sigOverlap=sigOverlap, overlapList=overlapList)
  if(!is.null(collapse)){
    out$overlapMatrix = overlapMatrix
  }

  return(out)
}

# Warning: Output significance levels are unadjusted.
OverlapListOfDfs = function(listOfDfs, listOfLists, bg1, bg2, colName, printFlag=T){
  #warning('typing of input DF column is hard-coded')
  m = length(listOfDfs)
  n = length(listOfLists)
  numOverlap = matrix(NA, nrow=m, ncol=n, dimnames=list(names(listOfDfs), names(listOfLists)))
  sigOverlap = matrix(0, nrow=m, ncol=n, dimnames=list(names(listOfDfs), names(listOfLists)))
  overlap = list()
  #numHypotheses = m*n

  for(i in 1:m){
    overlap[[i]] = list()
    for(j in 1:n){
      bg = intersect(bg1[[i]], bg2[[j]])
      a = as.integer(listOfDfs[[i]][, colName])
      A = intersect(bg, a)
      B = intersect(bg, listOfLists[[j]])
      numOverlap[i,j] = length(intersect(A,B))
      if(printFlag){
        print('')
        print(names(listOfDfs)[i])
        print(names(listOfLists)[j])
      }
      sigOverlap[i,j] = FisherExact(A, B, n_universe=length(bg),
                                     printFlag=printFlag)$adjp
      overlapIdx = na.omit(match(intersect(A, B), a))
      overlap[[i]][[j]] = listOfDfs[[i]][overlapIdx, ]
    }
    names(overlap[[i]]) = names(listOfLists)
  }
  names(overlap) = names(listOfDfs)
  out = list(numOverlap=numOverlap, sigOverlap=sigOverlap, overlap=overlap)
}

GeometricMean = function(x, na.rm=TRUE){
  if(na.rm){
    l = length(na.omit(x))
  }else{
    l = length(x)
  }
  return(exp(sum(log(x[x > 0]), na.rm=na.rm) / l))
}

CosineDistance = function(x, y, normalize=TRUE, na.rm=TRUE){
  if(na.rm){
    idx1 = which(!is.na(x))
    idx2 = which(!is.na(y))
    idx = intersect(idx1,idx2)
    x = x[idx]
    y = y[idx]
  }
  d = sum(x*y)
  if(normalize){
    d = d/(sqrt(sum(x*x))*sqrt(sum(y*y)))
  }
  return(1-d)
}

# return a matrix with i,j entry corresponding to the cosine distance between rows i and j
CosineDistanceMatrix = function(M, normalize=TRUE, asVector=FALSE){
  D = M %*% t(M)

  if(normalize){
    norm = apply(M, 1, Norm2)
    N = norm %*% t(norm)
    D = D / N
  }
  for(i in 1:nrow(D)){
    for(j in 1:i){
      D[i,j] = NA
    }
  }
  if(asVector){
    D = na.omit(as.vector(D))
  }
  return(1-D)
}

RescaleVec = function(x, a=0, b=1, abs=FALSE, na.rm=TRUE){
  if(abs){
    x = abs(x)
  }
  r = range(x, na.rm = na.rm)
  y = a + (x - r[1])*(b - a) / (r[2] - r[1])
  return(y)
}

CorThresh = function(pred, true, thresh=0){
  idx = which(abs(true) >= thresh)
  true = true[idx]
  pred = pred[idx]
  if(length(idx) > 0){
    out = cor(pred, true, use='pairwise.complete')
  }else{
    out = NA
  }
  return(out)
}


# S1 and S2 should be chemical structure matrices, each with drugs along the
# rows, and structural features along the columns. If S1 has n drugs and S2 has
# m drugs, then the output will be an nxm symmetric matrix with entries in
# [0,1].
JaccardIndexMatrix = function(S1, S2=NULL){
  symmetric = FALSE
  if(is.null(S2)){
    S2 = S1
    symmetric = TRUE
  }

  if(identical(S1, S2)){
    symmetric = TRUE
  }

  n = nrow(S1)
  m = nrow(S2)
  J = array(data=NA, dim=c(n,m), dimnames=list(drug1=rownames(S1), drug2=rownames(S2)))

  if(symmetric){
    for(i in 1:n){
      for(j in IncreasingSequence(i+1,m)){
        J[i,j] = JaccardIndex(S1[i,], S2[j,], type='binary')
      }
    }
    J = SymmetrifyMatrix(J, diag=1)
  }else{
    for(i in 1:n){
      for(j in 1:m){
        J[i,j] = JaccardIndex(S1[i,], S2[j,], type='binary')
      }
    }
  }
  return(J)
}

JaccardIndex = function(arg1, arg2, type='binary'){
  if(type == 'binary'){
    out = JaccardIndexBinary(arg1, arg2)
  }else if(type == 'set'){
    out = JaccardIndexSet(arg1, arg2)
  }else{
    stop('unexpected type in JaccardIndex')
  }
  return(out)
}

JaccardIndexBinary = function(binary_vec1, binary_vec2){
  stopifnot(length(binary_vec1) == length(binary_vec2))

  idx1 = which(binary_vec1 == 1)
  idx2 = which(binary_vec2 == 1)

  A_intersect_B = length(intersect(idx1, idx2))
  A_union_B = length(union(idx1, idx2))

  return(A_intersect_B / A_union_B)
}

JaccardIndexSet = function(A,B){
  S = union(A,B)
  vec1 = S %in% A
  vec2 = S %in% B
  return(JaccardIndexBinary(vec1, vec2))
}

LogNoInf = function(x, base=10){
  stopifnot(all(na.omit(x) >= 0))
  if(!(all(na.omit(x)>0))){
    idx = which(x==0)
    replaceValue = min(x[-idx], na.rm=TRUE)/2
    x[idx] = replaceValue
  }
  return(log(x, base=base))
}

PairwiseAssociation = function(X, Y, print=FALSE, binary=FALSE, adj='matrix'){

  tmp = matrix(data=NA, nrow=ncol(X), ncol=ncol(Y), dimnames=list(colnames(X), colnames(Y)))
  C = tmp
  P = tmp
  A = tmp

  for(j in 1:ncol(Y)){
    if(print)print(j)
    for(i in 1:ncol(X)){
      if(binary && IsBinary(X[,i])){
        list[P[i,j], C[i,j]] = TTest(Y[X[,i]==0,j], Y[X[,i]==1,j])
      }else{
        out = cor.test(X[,i],Y[,j], use='pairwise')
        P[i,j] = out$p.value
        C[i,j] = out$estimate
      }
    }
  }

  if(adj == 'matrix'){
    AdjP = AdjustPMatrix(P)
  }else if(adj == 'col'){
    AdjP = apply(P, 2, p.adjust)
  }else{
    stop('unrecognized argument for adj')
  }


  return(list(C=C, P=P, AdjP=AdjP))
}

TTest = function(x,y, print=TRUE){

  out = tryCatch({
    a = t.test(x,y)
    list(p.value=a$p.value, estimate=a$estimate[2] - a$estimate[1], message='')
  }, error=function(cond){
    if(print){print(cond)}
    return(list(p.value=NA, estimate=NA, message=cond))
  }, warning=function(cond){
    if(print){print(cond)}
    return(list(p.value=NA, estimate=NA, message=cond))
  }
  )
  return(out)
}

IsApprox = function(a, b, thresh=1e-6){
  return(abs(a-b) < thresh)
}

#### data frame helper functions ###################################################################

HeadDf = function(df){
  a = min(nrow(df), 5)
  b = min(ncol(df), 5)
  print(df[1:a, 1:b])
}

FixNAStrings = function(df, printFlag=FALSE){
  for(i in 1:ncol(df)){
    idx = which(df[,i] == 'NA')
    df[idx,i] = NA
    if(printFlag){
      print(sprintf('fix na strings %d',i))
      print(idx)
    }
  }
  return(df)
}

FixBackSlashes = function(df){
  for(i in 1:ncol(df)){
    if(class(df[,i]) == 'character'){
      df[,i] = gsub(df[,i], pattern='\\\\', replacement='\\\\\\\\\\\\\\\\')
    }
  }
  return(df)
}

ChangeColumnName = function(df, from, to){

  stopifnot(length(from) == length(to))
  if(length(from) > 1){
    for(i in 1:length(from)){
      df = ChangeColumnName(df, from[i], to[i])
    }
  }

  idx = which(names(df) == from)
  if(length(idx) != 1){
    warning(sprintf('found %d column name matches for %s', length(idx), from))
  }
  names(df)[idx] = to
  return(df)
}


RemoveDfColumns = function(df, columnNames){
  stopifnot(all(columnNames %in% names(df)))
  idx = which(names(df) %in% columnNames)
  return(df[,-idx, drop=F])
}

RemoveDfRows = function(df, rowNames){
  stopifnot(all(rowNames %in% rownames(df)))
  idx = which(rownames(df) %in% rowNames)
  return(df[-idx,])
}

ChangeDfCase = function(df, case){
  if(any(sapply(df, is.factor))){
    warning('some columns of data frame are factors, may want to run Factor2Char first')
  }
  idx = sapply(df, is.character)
  if(case == 'upper'){
    df[idx] = lapply(df[idx], toupper)
  }else if(case == 'lower'){
    df[idx] = lapply(df[idx], tolower)
  }else{
    warning('unexpected value for case, doing nothing.')
  }
  return(df)
}

Factor2Char = function(df){
  idx = sapply(df, is.factor)
  df[idx] = lapply(df[idx], as.character)
  return(df)
}

MoveColumn = function(df, colName, newPosition){
  currentPosition = which(colName == names(df))
  stopifnot(length(currentPosition) == 1)
  if(currentPosition > newPosition){ # insert column backwards
    for(i in currentPosition:(newPosition+1)){
      df[, c(i-1, i)] = df[,c(i, i-1)]
      names(df)[c(i-1,i)] = names(df)[c(i, i-1)]
    }
  }else if (currentPosition < newPosition){ # insert column forwards
    for(i in currentPosition:(newPosition-1)){
      df[, c(i, i+1)] = df[,c(i+1, i)]
      names(df)[c(i,i+1)] = names(df)[c(i+1, i)]
    }
  }
  return(df)
}

RemoveDuplicateRowsMulti = function(df, colNames){
  for(colName in colNames){
    df = RemoveDuplicateRows(df, colName)
  }
  return(df)
}

RemoveDuplicateRows = function(df, colName){
  x =  df[,colName]
  dups = unique(x[duplicated(x)])
  idxDups = which(x %in% dups)
  if(length(idxDups) > 0){
    out = df[-idxDups,]
  }else{
    out = df
  }
  return(out)
}

CollapseDuplicateRows = function(df, colName){
  x = df[,colName]
  idxDups = which(duplicated(x))
  if(length(idxDups) > 0){
    out = df[-idxDups,]
  }else{
    out = df
  }
  return(out)
}

MergeDfList = function(dfList, by, allRows=TRUE, allSuffixes=TRUE){
  if(allSuffixes){
    for(i in 1:length(dfList)){
      idx = which(names(dfList[[i]]) %in% by)
      names(dfList[[i]])[-idx] = paste(names(dfList[[i]])[-idx], names(dfList)[i], sep='.')
    }
  }
  mergedDf = as.data.frame(Reduce(function(x, y) merge(x, y, by=by, all=allRows), dfList))
  return(FixNAStrings(mergedDf))
}

CompareDfSubset = function(df1, df2, colNames=intersect(names(df1), names(df2)), matchName=NULL,
                            printFlag=T, maxThresh=0, normThresh=0){

  stopifnot(all(colNames %in% names(df1)))
  stopifnot(all(colNames %in% names(df2)))

  namesA = setdiff(names(df1), names(df2))
  namesB = setdiff(names(df2), names(df1))

  df1 = df1[,colNames]
  df2 = df2[,colNames]

  if(!is.null(matchName)){
    idx = match(df1[,matchName], df2[,matchName])
    df2 = df2[idx,]
  }

  same = TRUE
  badNames = c()
  if(!identical(df1, df2)){
    for(name in colNames){
      if(!all(na.omit(df1[,name]) == na.omit(df2[,name]))){
        misMatches = which(na.omit(df1[,name]) != na.omit(df2[,name]))
        if(length(misMatches)>0){
          if(printFlag){print(name); cat('  mismatches: '); cat(misMatches); cat('\n')}
         # if(!VectorCompare(df1[,name], df2[,name], maxThresh=maxThresh, normThresh=normThresh)){
            #warning('Not running VectorCompare')
            same = FALSE
            badNames = c(badNames, name)
        #  }
        }
      }
    }
  }

  goodNames = setdiff(colNames, badNames)

  return(list(same=same, badNames=badNames, goodNames=goodNames, namesA=namesA, namesB=namesB))
}

DataSummary = function(data, varname, groupnames){
  summary_func = function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))}
  data_sum = ddply(data, groupnames, .fun=summary_func, varname)
  data_sum = ChangeColumnName(data_sum, from='mean', to=varname)
  return(data_sum)
}

Write2XLS = function(dfnames, file, AdjWidth=TRUE, rownames=FALSE, colnames=TRUE, FreezeCol=0, FreezeRow=1){
  library(WriteXLS)
  WriteXLS(dfnames, row.names=rownames, col.names=colnames, AdjWidth=AdjWidth,
           FreezeCol=FreezeCol, FreezeRow=FreezeRow, BoldHeaderRow = TRUE,
           ExcelFileName=file)
}

ReadXLS = function(filename, sheetname=NULL) {
  sheets = readxl::excel_sheets(filename)
  if(length(sheets) > 1 && tail(sheets, 1) == 'Sheet1'){
    warning('Removing the last sheet, as it is called Sheet1')
    sheets = head(sheets, -1)
  }
  x = lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  names(x) = sheets
  if(!is.null(sheetname)){
    if(sheetname %in% sheets){
      out = x[[sheetname]]
    }else{
      warning('sheetname not found; returning all sheets.')
    }
  }else{
    out = x
  }

  if(length(out) == 1){
    out = out[[1]]
  }
  return(out)
}

FixMissingRownames = function(df){
  rownames(df) = ResolveDuplicateNAs(rownames(df))
  return(df)
}


#### plotting helper functions #####################################################################

Factor2Color = function(factorVar){
  n = length(unique(factorVar))
  col = rainbow(n)[as.integer(factorVar)]
  names(col) = factorVar
  idx = which(!duplicated(col))
  mapping = list(title=names(col[idx]), col=as.vector(col[idx]))
  return(list(mapping=mapping, col=col))
}

PasteOrCollapse = function(strings, collapse='|', na.omit=FALSE, unique=FALSE){
  if(na.omit){
    strings = na.omit(strings)
  }

  if(length(unique(strings))==0){
    out = NA
  }else if(length(unique(strings))==1){
    out = unique(strings)
  }else{
    if(unique){
      out = paste(unique(strings), collapse=collapse)
    }else{
      out = paste(strings, collapse=collapse)
    }
  }
  return(out)
}

MultiPlot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#### miscellaneous #################################################################################

Eval = function(string){
  return(eval(parse(text=string)))
}

InitRand = function(seed=123){
  set.seed(seed)
}

CheckInputString = function(input, validNames){
  if(!all(input %in% validNames)){
    stop('invalid input')
  }
}

Str2Vec = function(strings, split=',', unique=FALSE){
  if(unique){
    out = lapply(strings, function(x) unique(unlist(strsplit(x, split=split))))
  }else{
    out = lapply(strings, function(x) unlist(strsplit(x, split=split)))
  }
  if(length(strings)==1){
    out = unlist(out)
  }
  if(!is.null(names(strings))){
    names(out) = names(strings)
  }
  return(out)
}

ResolveDuplicateNAs = function(charVec){
  idx = which(is.na(charVec))
  charVec[idx] = paste0('NA.', 1:length(idx))
  return(charVec)
}

# input df should be an empty data frame with appropriate column names and data types to correspond to the input list of lists
List2Df = function(list, df){
  stopifnot(length(which(names(list[[1]]) %in% names(df))) > 0)

  for(i in 1:length(list)){
    for(name in names(df)){
      if(name %in% names(list[[i]])){
        df[i,name] = list[[i]][[name]]
      }else{
        df[i,name] = NA
      }
    }
  }
  return(df)
}

Ind2Sub = function(num.columns, ind){
  c = ((ind-1) %% num.columns) + 1
  r = floor((ind-1) / num.columns) + 1
  return(cbind(r,c))
}

Printf = function(text, level){
  debug = 2 # this should really be set somewhere else
  if(level <= debug){
    cat(text)
    return(text)
  }else{
    return('')
  }
}

IncreasingSequence = function(from, to){
  if (to >= from){
    out = from:to
  }else{
    out = seq(1,1,length=0)
  }
  return(out)
}

GetVarianceQuantile = function(M, quantile){
  v = rowVars(M)
  idx = which(v > quantile(v, probs=(1-quantile)))
  names(idx) = NULL
  return(idx)
}

GetVarianceTopK = function(M, k){
  v = rowVars(M)
  return(which(rank(-v) <= k))
}

GetMedianTopK = function(M, k){
  m = rowMedians(M)
  return(which(rank(-m) <= k))
}

MakeStringsUnique = function(strings){
  unq = unique(strings)
  for(u in unq){
    idx = which(strings == u)
    if(length(idx) > 1){
      strings[idx] = paste(strings[idx], 1:length(idx), sep='.')
    }
  }
  return(strings)
}

GetRowMinAndIdx = function(v){
  m = apply(v, 1, min)
  idx = apply(v, 1, function(x) which(x == min(x))[1])
  return(list(min=m, idx=idx))
}

IsSymmetric = function(A, thresh=1e-12){
  symmetric = TRUE
  if(nrow(A) != ncol(A)){
    symmetric = FALSE
  }
  a = A - t(A)
  if(norm(a) > thresh){
    symmetric = FALSE
  }
  return(symmetric)
}

# SetupParallel = function(nCores=1){
#   library(doMC)
#   library(doRNG)
#   registerDoMC()
#   options(cores=nCores)
# }

StartCluster = function(nCores=10){
  print(sprintf('starting cluster with %d cores', nCores))
  library(doParallel)
  cl = makePSOCKcluster(nCores)
  registerDoParallel(cl)
  options(cores=nCores)
  return(cl)
}

StopCluster = function(cl){
  print('stopping cluster')
  stopCluster(cl)
  options(cores=NULL)
  registerDoSEQ()
}

#### vector manipulation functions #################################

Norm2 = function(x){
  return(sqrt(sum(x^2)))
}

VectorCompare = function(x, y, normThresh=1e-12, maxThresh=NULL){
  isSame = TRUE
  if((length(x) != length(y)) || (is.na(x) != is.na(y)) ||
       class(x) != class(y) || any(is.finite(x) != is.finite(y)) ){
    isSame = FALSE
  }

  x = na.omit(x)
  y = na.omit(y)

  x = x[is.finite(x)]
  y = y[is.finite(y)]

  if(length(x) != length(y)){
    isSame = FALSE
  }

  if(isSame && is.numeric(x) && is.numeric(y) && length(x) > 0){

    if(is.null(normThresh) && is.null(maxThresh)){
      stop('both normThresh and maxThresh can\'t be null')
    }

    if(!is.null(normThresh)){
      if(Norm2(x - y) > normThresh){
        isSame = FALSE
      }
    }

    if(!is.null(maxThresh)){
      if(max(abs(x - y)) > maxThresh){
        isSame = FALSE
      }
    }

  }else if(isSame && !all(x == y)){
    isSame = FALSE
  }

  return(isSame)
}

EmptyVec = function(mode){
  return(vector(mode=mode, length=0))
}

FindAllDuplicates = function(v){
  return(which(duplicated(v) | duplicated(v, fromLast=TRUE)))
}

#### matrix manipulation functions #################################
MatrixCast = function(M, type){
  out = t(apply(M, 1, eval(parse(text=paste0('as.', type)))))
  stopifnot(dim(out) == dim(M))
  rownames(out) = rownames(M)
  colnames(out) = colnames(M)
  return(out)
}

SymmetrifyMatrix = function(M, diag=NA){
  if(nrow(M) != ncol(M)){
    warning('matrix is not symmetric, cannot symmetrify')
    out = M
  }else{
    M[is.na(M)] = 0
    M_t = t(M)
    diag(M_t) = 0
    out = M + M_t
    if(is.numeric(diag)){
      diag(out) = diag
    }
  }
  return(out)
}

GetUpperTriVec = function(M){
  return(as.vector(M[upper.tri(M, diag = FALSE)]))
}

MatrixCor = function(M1, M2, symmetric=TRUE){
  stopifnot(dim(M1) == dim(M2))
  if(symmetric){
    M1 = GetUpperTriVec(M1)
    M2 = GetUpperTriVec(M2)
  }
  return(cor(as.vector(M1), as.vector(M2), use='pairwise'))
}

CorMatrixList = function(list1, list2){
  stopifnot(length(list1) == length(list2))
  out = sapply(1:length(list1), function(i) MatrixCor(list1[[i]], list2[[i]]))
  names(out) = names(list1)
  return(out)
}

AdjustPMatrix = function(P){
  return(matrix(p.adjust(as.numeric(P), method='BH'),
                nrow=nrow(P),
                ncol=ncol(P),
                dimnames=dimnames(P)))
}

ClusterRows = function(M){
  hcr = hclust(dist(M))
  ddr = as.dendrogram(hcr)
  Rowv = rowMeans(M, na.rm=TRUE)
  ddr = reorder(ddr, Rowv)
  rowInd = order.dendrogram(ddr)
  return(rowInd)
}

ClusterCols = function(M){
  return(ClusterRows(t(M)))
}

ClusterRowsAndCols = function(M){
  return(list(rowInd=ClusterRows(M), colInd=ClusterCols(M)))
}

# remove rows and columns that are completely NA
MatrixRemoveNA = function(M, print=TRUE){
  list[nr, nc] = dim(M)
  M = M[rowSums(!is.na(M)) > 0,]
  M = M[,colSums(!is.na(M)) > 0]
  list[nr2, nc2] = dim(M)
  if(print){
    print(sprintf('Matrix dimensions went from (%s, %s) to (%s, %s)', nr, nc, nr2, nc2))
  }
  return(M)
}

