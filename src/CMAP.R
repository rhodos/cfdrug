LoadCMAP = function(file){
  wb = loadWorkbook(file)
  A = readWorksheet(wb, sheet=1)
  A = FixCMAPColNames(A)
  A = AddRank(A, colName='meta_rank')
  A = AddSign(A, colName='score')
  return(A)
}

FixCMAPColNames = function(A){
  n = names(A)[which(grepl('score.',names(A)))]
  stopifnot(length(n) == 3)
  numbers = sapply(n, function(x) strsplit(x, split='[.]')[[1]][2])
  stopifnot(length(numbers)==3)
  for(i in 1:3){
    pattern = sprintf('[.]%s', numbers[i])
    replacement = sprintf('.%d', i)
    names(A) = gsub(pattern, replacement, names(A))
  }
  A$unique_id = A$compound
  return(A)
}
 
GetCMAPAnnot = function(){
  n = load(DataDir('cmap_compound_metadata.RData'))
  stopifnot(n == 'annot')
  annot$unique_id = annot$compound
  annot = FixNAStrings(annot)
  return(annot)
}

GetCMAPColNames = function(){
  from = c('compound', 'chembl_id', 'chebi_id', 'DrugBank', 'PubChem', 'MeSH_ID', 
            'KEGG_drug', 'SMILES', 'InChI')
  to = c('COMPOUND', 'chembl', 'chebi', 'drugbank', 'pubchem', 'mesh', 
          'kegg', 'smiles', 'inchi')
  return(list(from=from, to=to))
}

GetCMAPAnnot4Map = function(){
  C = GetCMAPAnnot()
  colNames = GetCMAPColNames()
  from = colNames$from
  to = colNames$to
  for(i in 1:length(from)){
    C = ChangeColumnName(C, from[i], to[i])
    print(to[i])
    C[to[i]] = unlist(lapply(1:nrow(C), function(x) toupper(C[x, to[i]]) ))
  }
  rownames(C) = C[,'unique_id']
  return(FixNAStrings(C[, to]))
}
