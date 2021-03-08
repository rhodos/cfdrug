library(WriteXLS)

scores = LoadScores()

cmap.all = merge(scores$cmap_ks, scores$cmap_xsum, by=c('compound', 'COMPOUND'), all=TRUE, 
                  suffixes=c('.ks', '.xsum'))
cmap.all = AppendAnnotation(cmap.all, database='cmap')
rownames(cmap.all) = cmap.all$compound

lincs.all = merge(scores$lincs_ks, scores$lincs_xsum, by=c('pert_id', 'pert_iname', 'PERT_INAME'),
                   all=TRUE, suffixes=c('.ks', '.xsum'))
lincs.all = AppendAnnotation(lincs.all, database='lincs')
rownames(lincs.all) = lincs.all$pert_id

perc1 = 0.015
perc2 = 0.15
topCell = GetTopDrugsPercentile(scores=scores, perc1=perc1, perc2=perc2, scoreName='top_cell_score')
multiCell = GetTopDrugsPercentile(scores=scores, perc1=perc1, perc2=perc2, scoreName='multi_cell_score')

cmapDrugList = list(top_cell_score.ks=topCell$cmap_ks, top_cell_score.xsum=topCell$cmap_xsum,
                     multi_cell_score.ks=multiCell$cmap_ks, multi_cell_score.xsum=multiCell$cmap_xsum)

lincsDrugList = list(top_cell_score.ks=topCell$lincs_ks, top_cell_score.xsum=topCell$lincs_xsum,
                      multi_cell_score.ks=multiCell$lincs_ks, multi_cell_score.xsum=multiCell$lincs_xsum)

cmapDrugsFlat = unique(unlist(cmapDrugList))
lincsDrugsFlat = unique(unlist(lincsDrugList))

colNames = c('top_cell_score.ks', 'multi_cell_score.ks','top_cell_score.xsum', 'multi_cell_score.xsum',
'best.sig.ks', 'best.sig.xsum', 'best.pathway.ks', 'best.pathway.xsum','cell_type.ks', 'cell_type.xsum')

### construct cmap selection df

cmap = cmap.all[cmapDrugsFlat, c('compound', colNames)]

cmap$selected_top_cell_score.ks = FALSE
cmap[cmapDrugList$top_cell_score.ks,'selected_top_cell_score.ks'] = TRUE

cmap$selected_multi_cell_score.ks = FALSE
cmap[cmapDrugList$multi_cell_score.ks, 'selected_multi_cell_score.ks'] = TRUE

cmap$selected_top_cell_score.xsum = FALSE
cmap[cmapDrugList$top_cell_score.xsum, 'selected_top_cell_score.xsum'] = TRUE

cmap$selected_multi_cell_score.xsum = FALSE
cmap[cmapDrugList$multi_cell_score.xsum, 'selected_multi_cell_score.xsum'] = TRUE

cmap$category.venn = sapply(1:nrow(cmap), function(i) GetVennCategory(cmap[i,])$category)
cmap$score.in.category.venn = sapply(1:nrow(cmap), function(i) GetVennCategory(cmap[i,])$score)

cmap = AppendAnnotation(cmap, 'cmap')

cmap = cmap[order(cmap$score.in.category.venn, decreasing=TRUE),]
cmap = cmap[order(cmap$category.venn),]

### construct lincs selection df

lincs=lincs.all[lincsDrugsFlat, c('pert_id', 'pert_iname', colNames)]

lincs$selected_top_cell_score.ks = FALSE
lincs[lincsDrugList$top_cell_score.ks,'selected_top_cell_score.ks'] = TRUE

lincs$selected_multi_cell_score.ks = FALSE
lincs[lincsDrugList$multi_cell_score.ks, 'selected_multi_cell_score.ks'] = TRUE

lincs$selected_top_cell_score.xsum = FALSE
lincs[lincsDrugList$top_cell_score.xsum, 'selected_top_cell_score.xsum'] = TRUE

lincs$selected_multi_cell_score.xsum = FALSE
lincs[lincsDrugList$multi_cell_score.xsum, 'selected_multi_cell_score.xsum'] = TRUE

lincs$category.venn = sapply(1:nrow(lincs), function(i) GetVennCategory(lincs[i,])$category)
lincs$score.in.category.venn = sapply(1:nrow(lincs), function(i) GetVennCategory(lincs[i,])$score)

lincs = AppendAnnotation(lincs, 'lincs')

lincs = lincs[order(lincs$score.in.category.venn, decreasing=TRUE),]
lincs = lincs[order(lincs$category.venn),]

## add cross reference pointers between top selections
ids = Str2Vec(lincs$matching.cmap.names, split='[|]')
topMatchList = lapply(ids, function(x) unlist(lapply(x, function(y) y %in% cmap$compound)))
lincs$matching.cmap.names.in.top = unlist(lapply(topMatchList, function(x) paste0(x, collapse='|')))

ids = Str2Vec(cmap$matching.lincs.ids, split='[|]')
topMatchList = lapply(ids, function(x) unlist(lapply(x, function(y) y %in% lincs$pert_id)))
cmap$matching.lincs.ids.in.top = unlist(lapply(topMatchList, function(x) paste0(x, collapse='|')))

lincsDups = lincs$pert_iname[which(duplicated(lincs$pert_iname))]

## save to file
x = list(cmap=cmap, lincs=FixBackSlashes(lincs))
WriteXLS('x', ExcelFileName=OutputDir('drug_rankings/CF_drugList.xlsx'), 
         row.names=TRUE, col.names=TRUE, AdjWidth=FALSE,
         BoldHeaderRow=TRUE, FreezeRow=1, FreezeCol=1)

y = list(cmap=cmap.all, lincs=FixBackSlashes(lincs.all))
WriteXLS('y', ExcelFileName=OutputDir('drug_rankings/CF_drugList_allScores.xlsx'), 
         row.names=TRUE, col.names=TRUE, AdjWidth=FALSE, 
         BoldHeaderRow=TRUE, FreezeRow=1, FreezeCol=1)


