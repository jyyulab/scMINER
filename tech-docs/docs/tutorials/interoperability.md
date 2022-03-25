# Interoperability between SparseEset, Seurat, SingleCellExperiment, and anndata

## SC data type Conversion
|             | SparseEset | Seurat | SingleCellExperiment | anndata |
| ----------- | ----------- | ----------- | ----------- | ----------- |
| **SparseEset** | x | SparseEset<br>Seurat | SparseEset<br>SingleCellExperiment | SparseEset<br>SingleCellExperiment<br>anndata2ri |
| **Seurat** | Seurat<br>SparseEset | x | Seurat<br>SingleCellExperiment | <https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html> |
| **SingleCellExperiment** | SingleCellExperiment<br>SparseEset | SingleCellExperiment<br>Seurat | x | anndata2ri |
| **anndata** |	anndata2ri<br>SingleCellExperiment<br>SparseEset | anndata2ri<br>SingleCellExperiment<br>Seurat | anndata2ri<br>SingleCellExperiment | x |

### Seurat
<https://cran.r-project.org/web/packages/Seurat/index.html>
### SingleCellExperiment
<https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html>
### anndata
<https://cran.r-project.org/web/packages/anndata/index.html>
### anndata2ri
<https://github.com/theislab/anndata2ri>

### And this is an alternative package that helps easy conversion of different single-cell data formats to each other.
<https://github.com/cellgeni/sceasy>

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
