# Interoperability between SparseEset, Seurat, and anndata
## SparseEset to Seurat
```R
# eset.log2 is the SparseEset with log2 normalized matrix; 
# eset.sel is the SparseEset with raw count matrix

scMINER2Seurat<-function(eset.log2,eset.sel){
  library(Seurat)
  library(SingleCellExperiment)
  library(scMINER)
  dup_inx<-which(duplicated(fData(eset.log2)$geneSymbol)) 
  if(length(dup_inx)>0){
    eset.log2<-eset.log2[-dup_inx,] 
    eset.sel<-eset.sel[-dup_inx,] 
  counts = as.matrix(exprs(eset.sel))
  rownames(counts)<-fData(eset.sel)$geneSymbol
  logcounts = as.matrix(exprs(eset.log2))
  rownames(logcounts)<-fData(eset.log2)$geneSymbol
  sce<-SingleCellExperiment(
    assays = list(counts = counts,data = logcounts), 
    colData = pData(eset.log2)
  )
  SeuratObj <- as.Seurat(sce, counts = "counts", data = "data")
  return(SeuratObj)
}

SeuratObject<-scMINER2Seurat(eset.log2,eset.sel)
```

## Seurat to SparseEset
```R
# Create SparseEset with log2 normalized matrix

library(scMINER)
meta.data<-SeuratObject@meta.data
feature.data<-data.frame(rownames(SeuratObject@assays$RNA@data))
colnames(feature.data)<-"geneSymbol"
rownames(feature.data)<-feature.data$geneSymbol
eset<-CreateSparseEset(data=SeuratObject@assays$RNA@data,meta.data = meta.data,feature.data = feature.data,add.meta = F)

```
