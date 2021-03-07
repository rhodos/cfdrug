

debug = FALSE
nCores = 50 # not used currently, see ComputeAllScores
calcAll = FALSE

ComputeAllScores(csMetric='ks', database='cmap', debug=debug, nCores=nCores, calcAll=calcAll)
ComputeAllScores(csMetric='xsum', database='cmap', debug=debug, nCores=nCores, calcAll=calcAll)
ComputeAllScores(csMetric='ks', database='lincs', debug=debug, nCores=nCores, calcAll=calcAll)
ComputeAllScores(csMetric='xsum', database='lincs', debug=debug, nCores=nCores, calcAll=calcAll)