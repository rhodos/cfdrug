
#### directory and file handling ###################################################################

TestPlotDir = function(){
  file1 = PlotDir('test')
  file2 = PlotDir('test', subdir='subdir')
  stopifnot(file.exists(dirname(file1)))
  stopifnot(file.exists(dirname(file2)))
}

TestBaseDir = function(){
  dirs = unlist(strsplit(BaseDir(), "/"))
  stopifnot(identical(dirs[length(dirs)], 'code'))
}

TestGetConfig = function(){
  cfg = GetConfig('DATAPATH')
  stopifnot(is.character(cfg))
  stopifnot(!grepl('=', cfg))
  stopifnot(cfg != '')
}

TestDateStr = function(){
  str = DateStr()
  stopifnot(nchar(str) == 10)
  stopifnot(substr(str, 1, 2) == '20')
}

TestCheckDir = function(){
  setwd(BaseDir())
  CheckDir()
}

TestDataDir = function(){
  dirs = unlist(strsplit(DataDir(), "/"))
  stopifnot(identical(dirs[length(dirs)], 'data'))
  stopifnot(DataDir('/expr') == paste0(DataDir(), '/expr'))
}

TestOutputDir = function(){
  stopifnot(file.exists(OutputDir()))
  stopifnot(OutputDir('test') == paste0(OutputDir(), 'test'))
}

#### workspace helper functions ####################################################################

TestLsVars = function(){
  i = 5
  stopifnot('i' %in% LsVars())
}

TestLsFcns = function(){
  stopifnot('TestLsFcns' %in% LsFcns())
}

TestGetVarName = function(){
  foo = 2
  bar = c(1,2,3)
  joe = list('a','b',foo)
  stopifnot(GetVarName(foo) == 'foo')
  stopifnot(GetVarName(bar) == 'bar')
  stopifnot(GetVarName(joe) == 'joe')
}

#### set operations and helper functions ###########################################################

TestIsPartition = function(){
  stopifnot(IsPartition(list(1:5, 11, 6:10), 1:11))
  stopifnot(!IsPartition(list(1:5, 11, 6:10), 1:10))
  stopifnot(!IsPartition(list(1:5, 11, 6:11), 1:10))
}

TestIsIntersection = function(){
  stopifnot(IsIntersection(1:5, 6:10) == FALSE)
  stopifnot(IsIntersection(1:5, 5:6) == TRUE)
}

TestMatchAll = function(){
  A = 1:4
  B = c(1, 3, 3, 2, 2, 5, 1)
  out = c('1-7', '4-5', '2-3', NA)
  stopifnot(identical(MatchAll(A,B,collapse='-'), out))
}

TestGetMultiples = function(){
  x = c(1, 1, 2, 3, 1, 2)
  stopifnot(identical(GetMultiples(x, 1), unique(x)))
  stopifnot(identical(GetMultiples(x, 2), c(1, 2)))
  stopifnot(identical(GetMultiples(x, 3), c(1)))
  stopifnot(identical(GetMultiples(x, 4), c()))
}

TestSymDiff = function(){
  A = c(1:10)
  B = c(8:20)
  stopifnot(SymDiff(A,B) == c(1:7,11:20))
}

#### math/stats stuff ##############################################################################
TestZScore = function(){
  z = ZScore(rnorm(1000, mean=2, sd=0.1))
  m = mean(z)
  s = sd(z)
  stopifnot(abs(m) < 1e-10)
  stopifnot(abs(s-1) < 1e-10)
}

# TestZScoreMatrix = function(){
#
# }

TestIsBinary = function(){
  a = c(1,0,0,1,0)
  stopifnot(IsBinary(a))
  b = c(1,2,3,4)
  stopifnot(!IsBinary(b))
  c = rnorm(10)
  stopifnot(!IsBinary(c))
  d = 'hello world'
  stopifnot(!IsBinary(d))
}

TestInInterval = function(){
  stopifnot(InInterval(1, 1:2) == TRUE)
  stopifnot(InInterval(1, 2:3) == FALSE)
  stopifnot(InInterval(-1, c(-2, 0.2)) == TRUE)
}

TestFisherExact = function(){
  A = c(1, 2)
  B = c(1, 3)
  out = FisherExact(A, B, n_universe=4, printFlag=F, alternative='two.sided')
  stopifnot(out$p == 1)
  stopifnot(out$adjp == 1)
  stopifnot(out$overlap == 1)

  out = FisherExact(A, B, n_universe=1000, printFlag=F, alternative='two.sided')
  stopifnot(abs(out$p - 0.004) < .0001)
  stopifnot(abs(out$adjp - 0.004) < .0001)
  stopifnot(out$overlap == 1)
}

TestOverlapListOfLists = function(){
  list1 = list()
  list2 = list()
  bgSmall = list()
  bgLarge = list()

  n = 6
  for(i in 1:n){
    list1[[i]] = 1:i
    list2[[i]] = 1:i
    bgSmall[[i]] = 1:n
    bgLarge[[i]] = 1:(2*n)
  }
  out = OverlapListOfLists(list1, list2, bgSmall, bgSmall, printFlag=F)
  stopifnot(IsSymmetric(out$numOverlap))
  stopifnot(IsSymmetric(out$sigOverlap))
  for(i in 1:n){
    stopifnot(!is.unsorted(out$numOverlap[i,]))
    stopifnot(!is.unsorted(-out$sigOverlap[1:i,i]))
    stopifnot(!is.unsorted(out$sigOverlap[i:n,i]))
  }

  out2 = OverlapListOfLists(list1, list2, bgLarge, bgLarge, printFlag=F)
  stopifnot(IsSymmetric(out2$numOverlap))
  stopifnot(IsSymmetric(out2$sigOverlap))
  stopifnot(all(out$numOverlap == out2$numOverlap))
  stopifnot(all(out$sigOverlap >= out2$sigOverlap))

  # test alternative input options
  out3 = OverlapListOfLists(list1, list2, bgSmall, printFlag=F)
  stopifnot(identical(out, out3))

  bgNotAList = 1:10
  bgList = lapply(1:length(list1), function(x) bgNotAList)
  outNotAList = OverlapListOfLists(list1, list2, bg1=bgNotAList, printFlag=F)
  outList = OverlapListOfLists(list1, list2, bg1=bgList, printFlag=F)
  outSmallest = OverlapListOfLists(list1, bg1=bgList, printFlag=F)

  stopifnot(identical(outList, outNotAList))
  stopifnot(identical(outList, outSmallest))

}

TestOverlapListOfDfs = function(){
  DF = list()
  L = list()
  bgSmall = list()
  bgLarge = list()
  n = 20
  for(i in 1:n){
    DF[[i]] = data.frame(a=1:i, b=paste0(1:i, '_bla'))
    L[[i]] = 1:i
    bgSmall[[i]] = 1:n
    bgLarge[[i]] = 1:(2*n)
  }
  out = OverlapListOfDfs(DF, L, bgSmall, bgSmall, colName='a', printFlag=F)
  stopifnot(IsSymmetric(out$numOverlap))
  stopifnot(IsSymmetric(out$sigOverlap))
  for(i in 1:n){
    stopifnot(!is.unsorted(out$numOverlap[i,]))
    stopifnot(!is.unsorted(-out$sigOverlap[1:i,i]))
    stopifnot(!is.unsorted(out$sigOverlap[i:n,i]))
  }

  out2 = OverlapListOfDfs(DF, L, bgLarge, bgLarge, colName='a', printFlag=F)
  stopifnot(IsSymmetric(out2$numOverlap))
  stopifnot(IsSymmetric(out2$sigOverlap))
  stopifnot(all(out$numOverlap == out2$numOverlap))
  stopifnot(all(out$sigOverlap >= out2$sigOverlap))
}

TestGeometricMean = function(){
  a = GeometricMean(c(13, 23, 12, 44, 55))
  stopifnot(abs(a - 24.41932) < 1e-4)
  b = GeometricMean(c(13, 23, 12, 44, 55, NA))
  stopifnot(abs(b - 24.41932) < 1e-4)
}

TestCosineDistance = function(){
  x = c(1,0)
  y = c(0,1)
  a = CosineDistance(x,y,normalize=FALSE)
  b = CosineDistance(x,y,normalize=TRUE)
  stopifnot(a==b)
  stopifnot(a == 1)
}

TestCosineDistanceMatrix = function(){
  M = matrix(data=c(1,0,0,1,-1,0,0,-1), nrow=4, ncol=2, byrow=TRUE)
  d = CosineDistanceMatrix(M, normalize=FALSE, asVector=FALSE)
  v = CosineDistanceMatrix(M, normalize=FALSE, asVector=TRUE)
  D = matrix(data=NA, nrow=4, ncol=4)
  D[1,] = c(NA, 1, 2, 1)
  D[2,] = c(NA, NA, 1, 2)
  D[3,] = c(NA, NA, NA, 1)
  stopifnot(identical(d, D))
  V = c(1, 2, 1, 1, 2, 1)
  stopifnot(all(v == V))
}

TestRescaleVec = function(){
  x = rnorm(100)
  y = RescaleVec(x, a=-5, b=-2, abs=FALSE)
  stopifnot(Norm2(range(y) - c(-5,-2)) < 1e-10)
  stopifnot(Norm2(cor(x,y) - 1) < 1e-10)

  y = RescaleVec(x, a=-2, b=-5, abs=FALSE)
  stopifnot(Norm2(range(y) - c(-5,-2)) < 1e-10)
  stopifnot(Norm2(cor(x,y)+1) < 1e-10)
}

TestCorThresh = function(){
  x = rnorm(10000)
  y = x + rnorm(10000, sd = 1/sqrt(abs(x)))
  ct = c()
  for(i in 0:2){
    ct[i+1] = CorThresh(y, x, thresh=i)
  }
  stopifnot(!is.unsorted(ct))
}

TestJaccardIndexMatrix = function(){
  S1 = matrix(data=c(1,0,0,1,1,1), nrow=3, ncol=2, byrow=TRUE)
  M = JaccardIndexMatrix(S1)
  stopifnot(dim(M) == c(3,3))
  stopifnot(IsSymmetric(M))
  stopifnot(diag(M) == c(1,1,1))
  stopifnot(M[1,2] == 0)
  stopifnot(M[1,3] == 0.5)
  stopifnot(M[2,3] == 0.5)

  M2 = JaccardIndexMatrix(S1, S1)
  stopifnot(identical(M, M2))

  S2 = matrix(data=c(0,0,0,1), nrow=2, ncol=2)
  M3 = JaccardIndexMatrix(S1, S2)
  stopifnot(dim(M3) == c(3,2))
  stopifnot(M3[1,1] == 0)
}

TestJaccardIndex = function(){
  a = as.logical(c(1, 0, 0, 0, 1))
  b = as.logical(c(0, 1, 1, 1, 0))
  stopifnot(JaccardIndex(a,b,'binary') == 0)
  stopifnot(JaccardIndex(which(a), which(b), 'set') == 0)
  b = as.logical(c(1, 0, 0, 0, 0))
  stopifnot(JaccardIndex(a,b,'binary') == 0.5)
  stopifnot(JaccardIndex(which(a), which(b), 'set') == 0.5)
  b = as.logical(c(1, 1, 1, 1, 1))
  stopifnot(JaccardIndex(a,b,'binary') == 0.4)
  stopifnot(JaccardIndex(which(a), which(b), 'set') == 0.4)
}

TestJaccardIndexBinary = function(){} # tested in TestJaccardIndex
TestJaccardIndexSet = function(){} # tested in TestJaccardIndex

TestLogNoInf = function(){
  x = c(0, 1, 0.2, 10, NA)
  out = LogNoInf(x)
  stopifnot(identical(out, c(-1, 0, log10(0.2), 1, NA)))
}

TestPairwiseAssociation = function(){
  InitRand()
  X = matrix(data=rnorm(1000), nrow=100, ncol=10)
  Y = matrix(data=rnorm(500), nrow=100, ncol=5)
  list[C,P,AdjP] = PairwiseAssociation(X,Y)
  stopifnot(dim(C) == c(10,5))
  stopifnot(all(abs(C) < 0.3))

  Y = as.matrix(X[,1]+X[,3]-X[,9])
  list[C,P,AdjP] = PairwiseAssociation(X,Y)
  stopifnot(all(abs(C[c(1,3,9)]) > 0.4))
  stopifnot(all(abs(C[-c(1,3,9)]) < 0.3))
}

TestTTest = function(){
  set.seed(123)
  x = rnorm(100)
  y = rnorm(100, mean=1)
  out = TTest(x,y)
  stopifnot(out$p.value < 1e-5)
  stopifnot(out$estimate > 0.5)
  stopifnot(out$message == '')

  out = TTest(x[1], y, print=FALSE)
  stopifnot(is.na(out$p.value))
  stopifnot(is.na(out$estimate))
  stopifnot(out$message != '')

  # I guess you only need two data points to run a T-test? I thought you needed 3...
  out = TTest(x[1:2], y, print=FALSE)
  stopifnot(out$message == '')
}

TestIsApprox = function(){
  stopifnot(IsApprox(1, 1+1e-10))
  stopifnot(!IsApprox(1, 1.003))
}

#### data frame helper functions ###################################################################

TestHeadDf = function(){} #one-liner, just prints info to screen

TestFixNAStrings = function(){
  A = data.frame(a=c(1, 'NA', 2), b=c(NA, NA, 'NA'))
  B = FixNAStrings(A, printFlag=F)
  stopifnot(which(is.na(A)) == c(4,5))
  stopifnot(which(is.na(B)) == c(2,4,5,6))
}

TestFixBackSlashes = function(){}

TestChangeColumnName = function(){
  A = data.frame(a=c(1, 'NA', 2), b=c(NA, NA, 'NA'))
  A = ChangeColumnName(A, 'b', 'hello')
  stopifnot(names(A) == c('a', 'hello'))
}

TestRemoveDfColumns = function(){
  df = MakeTestDf()
  df1 = RemoveDfColumns(df, c('col3'))
  stopifnot(identical(df1, df[,-3]))

  df2 = RemoveDfColumns(df, c('col1', 'col3'))
  stopifnot(identical(df2, df[,2, drop=F]))

  df3 = RemoveDfColumns(df, c('col1','col2','col3'))
  stopifnot(identical(rownames(df3), c('row1', 'row2', 'row3')))
  stopifnot(all(dim(df3) == c(3,0)))
}

TestRemoveDfRows = function(){
  df = MakeTestDf()
  df1 = RemoveDfRows(df, c('row3'))
  stopifnot(identical(df1, df[-3,]))

  df2 = RemoveDfRows(df, c('row1', 'row3'))
  stopifnot(identical(df2, df[2, ,drop=F]))

  df3 = RemoveDfRows(df, c('row1','row2','row3'))
  stopifnot(identical(colnames(df3), c('col1','col2','col3')))
  stopifnot(all(dim(df3) == c(0,3)))
}

TestChangeDfCase = function(){
  df = Factor2Char(data.frame(a=1:3,b=c('a','b','c')))
  DF = ChangeDfCase(df, 'upper')
  stopifnot(DF$a == 1:3)
  stopifnot(DF$b == c('A','B','C'))
  stopifnot(identical(ChangeDfCase(DF,'lower'), df))
}

TestFactor2Char = function(){
  A = data.frame(a=1:3, b=as.factor(c('a', 'a', NA)))
  B = Factor2Char(A)
  stopifnot(is.integer(B$a))
  stopifnot(is.character(B$b))
  stopifnot(is.na(B[3,'b']))
}

TestMoveColumn = function(){
  A = data.frame(a=c(1, 1, 2), b=c(2, 3, 4), c=c(1, 2, 2), d=c(3, 4, 3))
  B = MoveColumn(A, colName='c', newPosition=1)
  stopifnot(names(B) == c('c', 'a', 'b', 'd'))
  C = MoveColumn(A, colName='c', newPosition=4)
  stopifnot(names(C) == c('a', 'b', 'd', 'c'))
  D = MoveColumn(A, colName='c', newPosition=3)
  stopifnot(names(D) == c('a', 'b', 'c', 'd'))
}

TestRemoveDuplicateRowsMulti = function(){
  A = data.frame(a=c(1, 1, 2, 2, 3, 4), b=c(1, 2, 2, 3, 4, 3))
  B = RemoveDuplicateRowsMulti(A, c('a', 'b'))
  stopifnot(identical(B$a, c(3, 4)))
  stopifnot(identical(B$b, c(4, 3)))
}

TestRemoveDuplicateRows = function(){} # tested by TestRemoveDuplicateRowsMulti()

TestCollapseDuplicateRows = function(){
 A1 = data.frame(a=1:3, b=rep('a',3))
 A2 = data.frame(a=1:3, b=rep('b',3))
 A3 = data.frame(a=c(1,1), b=c('a','b'))
 A = rbind(A1, A2)
 B = CollapseDuplicateRows(A, colName='a')
 CompareDfs(B, A1)
 C = CollapseDuplicateRows(A, colName='b')
 CompareDfs(C, A3)
}

TestMergeDfList = function(){
  dfList = list(one=data.frame(a=1:5, b=c('a','b','c','d','e')),
                two=data.frame(a=5:1, c=6:10),
                three=data.frame(a=1:4, d=1:4))
  out = MergeDfList(dfList, by='a')
  stopifnot(nrow(out)==5)
  stopifnot(identical(names(out),c('a', 'b.one', 'c.two', 'd.three')))
  stopifnot(identical(out$c.two,10:6))
}

TestCompareDfSubset = function(){
  A = data.frame(a=1:3, b=4:6, c=c('a','b','c'))
  B = A[3:1,c('a','c')]
  stopifnot(CompareDfSubset(A,B, printFlag=F)$same==FALSE)
  stopifnot(CompareDfSubset(A,B, printFlag=F, matchName='a')$same)
  stopifnot(CompareDfSubset(A,B, printFlag=F, matchName='c')$same)
}

TestDataSummary = function(){
  n = 10000
  data = melt(data.frame(x=rnorm(n), y=rnorm(n, mean=2)), measure.vars=c('x','y'))
  out = DataSummary(data, varname='value', groupnames='variable')
  rownames(out) = out$variable
  thresh = 0.05
  stopifnot(abs(out['x','value']) < thresh)
  stopifnot(abs(out['y','value']-2) < thresh)
  stopifnot(abs(out['x','sd']-1) < thresh)
  stopifnot(abs(out['y','sd']-1) < thresh)
}

TestWrite2XLS = function(){} # simple wrapper for another function
TestReadXLS = function(){}

TestFixMissingRownames = function(){
  df = MakeTestDf()
  C = as.matrix(df)
  rownames(C) = c('a', NA, NA)
  df = as.data.frame(C)
  df = FixMissingRownames(df)
  stopifnot(!any(is.na(rownames(df))))
}

#### plotting helper functions #####################################################################

TestFactor2Color = function(){
  fac = as.factor(c(1:3,1))
  out = Factor2Color(fac)
  stopifnot(length(unique(out$col)) == 3)
  stopifnot(out$col[1] == out$col[4])
}

TestPasteOrCollapse = function(){
  stopifnot(PasteOrCollapse(c('a','a','a'))=='a')
  stopifnot(PasteOrCollapse(c('a','b','a'))=='a|b|a')
  stopifnot(PasteOrCollapse(c('a','b',NA), na.omit=TRUE) == c('a|b'))
}

#### miscellaneous #################################################################################

TestEval = function(){
  stopifnot(Eval('2+2') == 4)
}

TestInitRand = function(){
  n = 1000
  InitRand()
  x1 = rnorm(n)
  x2 = rnorm(n)
  stopifnot(!identical(x1,x2))
  InitRand()
  x2 = rnorm(n)
  stopifnot(identical(x1,x2))
}

TestCheckInputString = function(){} # one-liner

TestStr2Vec = function(){
  str1 = 'a,b,foo,bar,a'
  strList1 = c('a','b','foo','bar','a')
  str2 = 'apple,peaches,pumpkin pie'
  strList2 = c('apple','peaches','pumpkin pie')
  stopifnot(identical(Str2Vec(str1), strList1))
  stopifnot(identical(Str2Vec(str1, unique=T), unique(strList1)))
  stopifnot(identical(Str2Vec(c(str1,str2)), list(strList1,strList2)))
}

TestResolveDuplicateNAs = function(){
  vec = c('a', NA, NA, 'b')
  vec2 = ResolveDuplicateNAs(vec)
  stopifnot(identical(vec2, c('a', 'NA.1', 'NA.2', 'b')))
}

TestList2Df = function(){
  list = list(a=list(A=1, B=2, C=3), b=list(A=2, C=4))
  df = data.frame(A=EmptyVec('numeric'), B=EmptyVec('numeric'), C=EmptyVec('numeric'))
  df1 = data.frame(A=c(1,2), B=c(2,NA), C=c(3,4))
  stopifnot(identical(List2Df(list, df), df1))
}

TestInd2Sub = function(){
  A = matrix(1:12, nrow=3, ncol=4)
  out = Ind2Sub(ncol(A), which(A<6))
  stopifnot(identical(out[,'r'], c(1,1,1,1,2)))
  stopifnot(identical(out[,'c'], c(1,2,3,4,1)))
}

TestPrintf = function(){
  str = 'Hello world!'
  stopifnot(Printf(str, 0) == str)
  stopifnot(Printf(str, 3) == '')
}

TestIncreasingSequence = function(){
  stopifnot(length(IncreasingSequence(5, 1)) == 0)
  stopifnot(length(IncreasingSequence(1, 5)) == 5)
}

TestGetVarianceQuantile = function(){
  # generate matrix with increasing variance in each row
  InitRand(1234)
  n = 1000
  X = replicate(n, rnorm(n))
  for(i in 1:n){
    X[i,] = X[i,]*i
  }
  idx = GetVarianceQuantile(X, quantile=0.1)
  stopifnot(!(1:50 %in% idx))
  stopifnot(950:1000 %in% idx)
}

TestGetVarianceTopK = function(){
  # generate matrix with increasing variance in each row
  n = 1000
  InitRand(1234)
  X = replicate(n, rnorm(n))
  for(i in 1:n){
    X[i,] = X[i,]*i
  }
  idx = GetVarianceTopK(X, 100)
  stopifnot(!(1:50 %in% idx))
  stopifnot(950:1000 %in% idx)
}

TestGetMedianTopK = function(){
  # generate matrix with increasing median in each row
  n = 1000
  InitRand(1234)
  X = replicate(n, rnorm(n))
  for(i in 1:n){
    X[i,] = X[i,]*i
  }
  idx = GetVarianceTopK(X, 100)
  stopifnot(!(1:50 %in% idx))
  stopifnot(950:1000 %in% idx)
}

TestMakeStringsUnique = function(){
  strings = c('a', 'a', 'b', 'cd')
  strings2 = MakeStringsUnique(strings)
  stopifnot(identical(strings2, c('a.1', 'a.2', 'b', 'cd')))
}

TestGetRowMinAndIdx = function(){
  v = cbind(c(1,2,3), c(0,3,2))
  out = GetRowMinAndIdx(v)
  stopifnot(identical(out$min, c(0, 2, 2)))
  stopifnot(all(out$idx==c(2,1,2)))
}

TestIsSymmetric = function(){
  A = replicate(10, rnorm(10))
  B = A + t(A)
  stopifnot(!IsSymmetric(A))
  stopifnot(IsSymmetric(B))
}

# TestSetupParallel = function(){
#   n=2
#   SetupParallel(nCores=n)
#   set.seed(42)
#   out1 = unlist(foreach(i=1:n) %dorng% {runif(1)})
#   set.seed(42)
#   out2 = unlist(foreach(i=1:n) %dorng% {runif(1)})
#   set.seed(NULL)
#   out3 = unlist(foreach(i=1:n) %dorng% {runif(1)})
#   if(!identical(out1, out2)){
#     warning('Parallel computation is not repeatable')
#   }
#   stopifnot(!identical(out1, out3))
# }

TestStartCluster = function(){}
TestStopCluster = function(){}

#### vector manipulation functions #################################################################
TestNorm2 = function(){
  x = c(2,2,sqrt(41))
  stopifnot(abs(Norm2(x) - 7) < 1e-6)
}

TestVectorCompare = function(){
  x = c(1, 2, 3)
  y = c(1, 2, 3) + 0.01
  stopifnot(!VectorCompare(x,y))
  stopifnot(VectorCompare(x,y,maxThresh=0.011, normThresh=NULL))
  stopifnot(VectorCompare(x, y, normThresh=sqrt(4e-4)))
}

TestEmptyVec = function(){
  a = EmptyVec('character')
  b = EmptyVec('numeric')
  stopifnot(length(a) == 0)
  stopifnot(length(b) == 0)
  stopifnot(class(a) == 'character')
  stopifnot(class(b) == 'numeric')
}

TestFindAllDuplicates = function(){
  a = c(letters[1:6], 'A', letters[1:6])
  stopifnot(a[-FindAllDuplicates(a)] == 'A')
}

#### matrix manipulation functions #################################################################
TestMatrixCast = function(){
  A = matrix(data=rep(c(0,1), 3), nrow=3, dimnames=list(letters[1:3], LETTERS[1:2]))

  B = MatrixCast(A, 'logical')
  C = MatrixCast(A, 'character')

  testB = matrix(data=as.logical(rep(c(0,1), 3)), nrow=3, dimnames=list(letters[1:3], LETTERS[1:2]))
  testC = matrix(data=as.character(rep(c(0,1), 3)), nrow=3, dimnames=list(letters[1:3], LETTERS[1:2]))

  stopifnot(identical(B, testB))
  stopifnot(identical(C, testC))
}

TestSymmetrifyMatrix = function(){
  n = 10
  A = matrix(data=NA, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in IncreasingSequence(i+1,n)){
      A[i,j] = rnorm(1)
    }
  }
  B = SymmetrifyMatrix(A, diag=1)
  stopifnot(all(B == t(B)))
  B[is.na(A)] = NA
  stopifnot(identical(A, B))
}

TestGetUpperTriVec = function(){
  A = matrix(data=1:9, nrow=3, ncol=3)
  a = GetUpperTriVec(A)
  stopifnot(a == c(4,7,8))
}

TestMatrixCor = function(){
  # perfectly correlated matrices
  M1 = matrix(data=rnorm(100), nrow=10, ncol=10)
  M2 = 2*M1
  stopifnot(IsApprox(MatrixCor(M1, M2, symmetric=FALSE),1))
  stopifnot(IsApprox(MatrixCor(M1, -M2, symmetric=FALSE),-1))

    # less correlated
  M2 = M1 + matrix(data=rnorm(100), nrow=10, ncol=10)
  stopifnot(MatrixCor(M1,M2, symmetric=FALSE) > 0.4)

  # no entries line up
  M2[upper.tri(M2, diag=TRUE)] = NA
  stopifnot(is.na(MatrixCor(M1,M2)))
  stopifnot(!is.na(MatrixCor(M1, M2, symmetric=FALSE)))
}

TestCorMatrixList = function(){
  n = 100
  M = matrix(data=rnorm(n*n), nrow=n, ncol=n)
  M1 = list()
  M2 = list()
  for(i in 1:10){
    M1[[i]] = M
    M2[[i]] = M + matrix(data=rnorm(n*n, sd=0.1*i), nrow=n, ncol=n)
  }
  out = CorMatrixList(M1, M2)
  stopifnot(!is.unsorted(-out))
}

TestAdjustPMatrix = function(){
  InitRand()
  n = 10000
  nHits = 100
  trueHits = sample(n,nHits)
  data = runif(n)
  data[trueHits] = data[trueHits]/(n * 20)
  P = matrix(data=data, nrow=sqrt(n), ncol=sqrt(n))
  adjP = AdjustPMatrix(P)

  # check that I get roughly the same hits
  predHits = which(as.numeric(adjP) < 0.05)
  stopifnot(InInterval(length(predHits), c(0.8, 1.2)*nHits))
  stopifnot(JaccardIndex(trueHits, predHits, 'set') > 0.8)

  overlap = intersect(trueHits, predHits)
  stopifnot(cor(adjP[overlap], P[overlap]) > 0.3)
}

# this works but it generates a plot so I don't want to run it every time I test
TestClusterRows = function(){
  # M = MakeTestMatrix()
  # hm = heatmap(M, labRow='', labCol='')
  # stopifnot(identical(hm$rowInd, ClusterRows(M)))
}

# this works but it generates a plot so I don't want to run it every time I test
TestClusterCols = function(){
  # M = MakeTestMatrix()
  # hm = heatmap(M, labRow='', labCol='')
  # stopifnot(identical(hm$colInd, ClusterCols(M)))
}

TestClusterRowsAndCols = function(){}

TestMatrixRemoveNA = function(){
  M = matrix(data=runif(25), nrow=5, ncol=5)
  M[1,1] = NA
  stopifnot(identical(M, MatrixRemoveNA(M, print=FALSE)))
  M[2,] = NA
  M[,c(4,5)] = NA
  M2 = MatrixRemoveNA(M, print=FALSE)
  stopifnot(dim(M2) == c(4,3))
}
