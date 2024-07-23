---
title: "scMINER Guided Analysis on 14K PBMCs from 10x Genomics"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{pbmc14k_tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# About the PBMC14K data set

This dataset was generated from a Peripheral Blood Mononuclear Cells (PBMCs) dataset containing 10 sorted populations of 2,000 cells per population [[Zheng et al., 2017](./data-reference.md#[Zheng et al., 2017])]. We first rectified the gene symbol issues of the original dataset, including the dash-dot conversion (e.g. "RP11-34P13.7" changed to "RP11.34P13.7") and "X" added to those started with numbers (e.g. "7SK" changed to "X7SK"), by referring to the gene annotation file (GRCh37.82) used in the original study. Then we removed 3 cell populations, CD34+ cells, CD4+ helper T cells, and total CD8+ cytotoxic cells, from the dataset because of low sorting purity or a significant overlap with other immune cells based on the sorting strategy, and created a new dataset with seven known cell types and 14k cells in total. The original dataset is freely available under accession number [SRP073767](https://www.ncbi.nlm.nih.gov/sra?term=SRP073767) and [Zenodo](https://zenodo.org/record/3357167#.YhQNF2RKj6V).

<details>
<summary>**How was the PBMC14K dataset generated from the original dataset?**</summary>

```r
## Step 1: rectify the invalid gene symbols 
counts <- read.csv("./Filtered_DownSampled_SortedPBMC_data.csv", row.names = 1) # "Filtered_DownSampled_SortedPBMC_data.csv" is the raw count matrix directly downloaded from Zenodo
d <- t(counts); dim(d) # it includes 21592 genes and 20000 cells

officialGene <- read.table("./genesymbol_from_GTF_GRCh37.txt", header = T, sep = "\t", quote = "", stringsAsFactors = F); head(officialGene) # "genesymbol_from_GTF_GRCh37.txt" contains the official gene ids and symbols extracted from GTF file downloaded from https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/
officialGene$dotted_symbol <- gsub("-", "\\.", officialGene$gene_name); officialGene$dotted_symbol <- make.unique(officialGene$dotted_symbol)
table(row.names(d) %in% officialGene$dotted_symbol); row.names(d)[! row.names(d) %in% officialGene$dotted_symbol] # two genes are not in: X7SK.1 and X7SK.2
row.names(d) <- gsub("X7SK.1", "7SK", row.names(d)); row.names(d) <- gsub("X7SK.2", "7SK.1", row.names(d))
table(row.names(d) %in% officialGene$dotted_symbol) # all true

row.names(officialGene) <- officialGene$dotted_symbol
officialGene <- officialGene[row.names(d),]
row.names(d) <- make.unique(officialGene$gene_name)
write.table(d, file = "./PBMC20K_rawCount.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE) # 21592 genes, 20000 cells

celltype <- read.csv("./Labels.csv"); head(celltype); # "Labels.csv" contains the true label of cell types directly downloaded from Zenodo
table(celltype$x) # 2000 cells for each of 10 cell types: CD14+ Monocyte, CD19+ B, CD34+, CD4+ T Helper2, CD4+/CD25 T Reg, CD4+/CD45RA+/CD25- Naive T, CD4+/CD45RO+ Memory, CD56+ NK, CD8+ Cytotoxic T, CD8+/CD45RA+ Naive Cytotoxic
df <- data.frame(cell_barcode = colnames(d), cell_type = celltype$x); dim(df) ## 
write.table(df, file = "./PBMC20K_trueLabel.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)

## Step 2: extract 7 populations
df.14k <- df[df$cell_type %in% c("CD14+ Monocyte", "CD19+ B", "CD4+/CD25 T Reg", "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory", "CD56+ NK", "CD8+/CD45RA+ Naive Cytotoxic"),]
write.table(df.14k, file = "./PBMC14K_trueLabel.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE)

d.14k <- d[,df.14k$cell_barcode]
d.14k <- d.14k[rowSums(d.14k) > 0,]
write.table(d.14k, file = "./PBMC14K_rawCount.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = FALSE) # 17986 genes, 14000 cells
```
</details>
\



```r
## Load the packages
library(scMINER)
library(ggplot2)

## Create the project space to run the scMINER analysis
createProjectSpace(project_dir = "/Users/qpan", project_name = "PBMC14k", do.unlink = TRUE)
#> A folder of the specified project name was found and has been cleaned up in the specified project directory
#> The project space has been successfully created: /Users/qpan/PBMC14k !
list.dirs("/Users/qpan/PBMC14k")
#> [1] "/Users/qpan/PBMC14k"          "/Users/qpan/PBMC14k/DATA"    
#> [3] "/Users/qpan/PBMC14k/MICA"     "/Users/qpan/PBMC14k/PLOT"    
#> [5] "/Users/qpan/PBMC14k/SJARACNe"
```

# Data Preprocessing

To handle the sparsity in single-cell RNA sequencing data better, scMINER created a new class called **SparseExpressionSet**. The SparseExpressionSet-class is derived from ExpressionSet, and can take a sparse matrix as assayData component. The purpose of data preprocessing is to generate a SparseExpressionSet object containing filtered and normalized single cell gene expression data.

## Create the SparseExpressionSet object

Similar to ExpressioinSet, the SparseExpressionSet-class contains these component parts:
- assayData: required, a `matrix` of expression or activity values, features by cells. scMINER supports these classes: '`dgCMatrix`', '`dgTMatrix`', '`dgeMatrix`', '`matrix`', '`data.frame`'.
- phenoData: optional, a `data frame` of cell meta data
- featureData: optional, a `data frame` of featute meta data
- annotation: optional, a `character` describing the path to the project folder or other user-defined property

scMINER provides a function, `createSparseEset()`, to assemble these component parts and create a SparseExpresionSet object. We will go through all of these component parts one by one below.

### Prepare the assayData

The `assayData` for scMINER's SparseExpressionSet object is usually but not limited to a sparse matrix of raw UMI counts. As for the data format, it's also compatible with the `regular matrix` and `data frame`. Then they were given, scMINER can automatically convert them to a sparse matrix if the `do.sparseConversion` argument of `createSparseEset()` is set `TRUE`. As for the type of quantification measures, it also takes Normalized counts (e.g. **CPM** or **CP10k**), **TPM** (Transcripts Per Million), **FPKM/RPKM** (Fragments/Reads Per Kilobase of transcript per Million) and others.

To make the preparation of assayData easier, scMINER provides four functions for this purpose:
- `readInput_10x.dir()`: to read the gene expression data from a directory containing three files generated by 10x Genomics
- `readInput_10x.h5()`: to read the gene expression data from the HDF5 file generated by CellRanger pipeline of 10x Genomics
- `readInput_h5ad()`: to read the gene expression data from the h5ad file which is commonly used in single cell data sharing
- `readInput_table()`: to read the gene expression data from a table-format file, e.g. tsv, txt and csv.

Here are a few samples to prepare the assayData using the functions above:

```r
## Read the input data generated by 10x Genomics from a directory
input_dir <- "path-to-directory"
list.files(input_dir, full.names = FALSE)  # you should see three files: matrix.mtx, barcodes.tsv and features.tsv (or genes.tsv)
assayData_matrix <- readInput_10x.dir(input_dir, featureType = "gene_symbol", removeSuffix = TRUE,
    addPrefix = "demoSample")

## Read the input data generated by 10x Genomics from the HDF5 file
hdf5_file <- "path-to-hdf5_file"
assayData_matrix <- readInput_10x.h5(hdf5_file, featureType = "gene_symbol", removeSuffix = TRUE,
    addPrefix = "demoSample")

## Read the h5ad file
h5ad_file <- "path-to-h5ad_file"
assayData_matrix <- readInput_h5ad(h5ad_file, removeSuffix = FALSE, addPrefix = "demoSample")

## Read the table format file
table_file <- "path-to-table_file"
assayData_matrix <- readInput_table(table_file, sep = "\t", is.geneBYcell = TRUE, removeSuffix = FALSE,
    addPrefix = "demoSample")
```

In this case, the raw count matrix of PBMC14k data has been added to the package and can be easily loaded by:


```r
data(pbmc14k_rawCount)
dim(pbmc14k_rawCount)
#> [1] 17986 14000
```

### Prepare the phenoData

The `phenoData` for scMINER's SparseExpressionSet object is a data frame of cell annotation information (a.k.a. meta data of cells). This is optional and can be set `NULL` in `createSparseEset()`. If so, scMINER will automatically generate a column nameed "`CellID`" using the column names of the `assayData` matrix. scMINER can also automatically add the project ID if the `projectID` argument is specified, and add four quality control metrics, "`nUMI`", "`nFeature`", "`pctMito`" and "`pctSpikeIn`", if the `addMetadata` argument is set `TRUE` in `createSparseEset()`. If there is any other cell annotations to add, the data frame of `phenoData` is the way to go.

In this case, we would like to add the true labels of cell types into the SparseExpressionSet object.

```r
## Generally, the data frame can be generated by: pbmc14k_trueLabel <-
## read.table('PBMC14K_trueLabel.txt', header = T, sep = '\t', quote = '', stringsAsFactors =
## F, row.names = 1)

## In this case, the true labels of cells has been added to the package and can be loaded by:
data(pbmc14k_trueLabel)

table(pbmc14k_trueLabel$true_label)  # check the true label
#> 
#>               CD14+ Monocyte                      CD19+ B 
#>                         2000                         2000 
#>              CD4+/CD25 T Reg   CD4+/CD45RA+/CD25- Naive T 
#>                         2000                         2000 
#>          CD4+/CD45RO+ Memory                     CD56+ NK 
#>                         2000                         2000 
#> CD8+/CD45RA+ Naive Cytotoxic 
#>                         2000
table(row.names(pbmc14k_trueLabel) == colnames(pbmc14k_rawCount))  # the row names of phenoData must match he column names of assayData
#> 
#>  TRUE 
#> 14000
```

### Prepare the featureData

The `featureData` for scMINER's SparseExpressionSet object is a data frame of feature/gene annotation information (a.k.a. meta data of features/genes). Similar to the `phenoData`, it is also optional and can be set `NULL` in `createSparseEset()`. If so, scMINER will automatically generate a column nameed "`GeneSymbol`" using the row names of the `assayData` matrix. scMINER can also automatically add one quality control metrics, "`nCell`", if the `addMetadata` argument is set `TRUE` in `createSparseEset()`. If there is any other feature annotations to add, like type of genes or genome coordinates, the data frame of `featureData` is the way to go.

In this case, we don't have more feature annotation to add and will go with `NULL`.

```r
featureData_info <- NULL
```

### Prepare the annotation

The `annotation` for scMINER's SparseExpressionSet object is a portal slot allowing the users to save some "key words" about the project. It's also optional and can be leave as `NULL`. However, it's **highly recommended** to set it as **the path to the project space**. scMINER always saves the output files in the directories specified by the users. However, if it's not specified, scMINER will automatically retrieve the path to the project space from the `annotation` slot. If a path is found and it does exists, scMINER will use it as the project space.


```r
annotation_info <- ".//PBMC14k"
```

### Generate the sparse eset object.
Now we are ready to generate the SparseExpressionSet object:

```r
pbmc14k_raw.eset <- createSparseEset(input_matrix = pbmc14k_rawCount, cellData = pbmc14k_trueLabel,
    featureData = featureData_info, annotation = annotation_info, addMetaData = TRUE)
#> Creating sparse eset from the input input_matrix ...
#> 	Adding meta data based on input_matrix ...
#> Done! The sparse eset has been generated: 17986 genes, 14000 cells.
pbmc14k_raw.eset
#> SparseExpressionSet (storageMode: environment)
#> assayData: 17986 features, 14000 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: CACTTTGACGCAAT GTTACGGAAACGAA ... ACGTGCCTTAAAGG (14000
#>     total)
#>   varLabels: cell_barcode true_label ... CellID (7 total)
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: AL627309.1 AP006222.2 ... SRSF10.1 (17986 total)
#>   fvarLabels: GeneSymbol nCell
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation: .//PBMC14k
saveRDS(pbmc14k_raw.eset, file = "/Users/qpan/PBMC14k/DATA/pbmc14k_raw.rds")
```

## QC and filter the SparseExpressionSet object

### QC metrics available in scMINER
scMINER provides 4 quality control metrics that [commonly used by the community](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/) to identifiy the high-quality cells, including:
- **"nUMI"**: number of total UMIs in each cell. Cells with abnormally high nUMI usually indicate doublets, while those with abnormally low nUMI usually indicate poorly sequenced cells or empty droplets.
- **"nFeature"**: number of expressed features/genes in each cell. Similar to nUMI.
- **"pctMito"**: percentage of UMIs of mitochondrial genes (defined by "^mt-|^MT-") in each cell. Cells with aberrantly high pctMito usually indicate dying cells.
- **"pctSpikeIn"**: percentage of UMIs of spike-in RNAs (defined by "^ERCC-|^Ercc-")) in each cell. This is used to estimate the normalization factor. Cells with extremely high or low pctSpikeIn need to be removed.

These 4 metrics can work together to identify the cells of poor quality. They are stored in `phenoData` and can be accessed by:

```r
head(pData(pbmc14k_raw.eset))
#>                  cell_barcode     true_label nUMI nFeature    pctMito
#> CACTTTGACGCAAT CACTTTGACGCAAT CD14+ Monocyte  764      354 0.01832461
#> GTTACGGAAACGAA GTTACGGAAACGAA CD14+ Monocyte  956      442 0.01569038
#> AGTCACGACAGGAG AGTCACGACAGGAG CD14+ Monocyte 7940     2163 0.01977330
#> TTCGAGGACCAGTA TTCGAGGACCAGTA CD14+ Monocyte 4177     1277 0.01149150
#> CACTTATGAGTCGT CACTTATGAGTCGT CD14+ Monocyte  629      323 0.02066773
#> GCATGTGATTCTGT GCATGTGATTCTGT CD14+ Monocyte  875      427 0.02628571
#>                pctSpikeIn         CellID
#> CACTTTGACGCAAT          0 CACTTTGACGCAAT
#> GTTACGGAAACGAA          0 GTTACGGAAACGAA
#> AGTCACGACAGGAG          0 AGTCACGACAGGAG
#> TTCGAGGACCAGTA          0 TTCGAGGACCAGTA
#> CACTTATGAGTCGT          0 CACTTATGAGTCGT
#> GCATGTGATTCTGT          0 GCATGTGATTCTGT
```

scMINER also provides one quality control metrix to identify the high-quality genes:
- **"nCell"**: number of cells that each feature/gene was identified. Genes with extremely low nCell indicate genes poorly sequenced and are usually of low variance.

The nCell metrics is stored in `featureData` and can be accessed by:

```r
head(fData(pbmc14k_raw.eset))
#>                  GeneSymbol nCell
#> AL627309.1       AL627309.1    50
#> AP006222.2       AP006222.2     2
#> RP11-206L10.3 RP11-206L10.3     1
#> RP11-206L10.2 RP11-206L10.2    33
#> RP11-206L10.9 RP11-206L10.9    17
#> LINC00115         LINC00115   115
```

### Generate the QC report
To help assess the data quality and determine the cutoffs for data filtration, scMINER provides a function, `drawSparseEsetQC()`, to generate a html-format quality control report:

```r
## scMINER supports group-specific QC highlights
drawSparseEsetQC(input_eset = pbmc14k_raw.eset, output_html_file = "/Users/qpan/PBMC14k/PLOT/pbmc14k_rawCount.html",
    overwrite = FALSE, group_by = "true_label")
```

The qaulity control report consists of 4 parts:
- **Key Statistics**: it highlights 5 key statistics of given eset object, including `number of cells`, `number of genes`, mean of `genes per cell`, mean of `UMIs per cell` and mean of `cells per gene`.
- **Detailed statistics of key metrics**: it summarizes and visualizes the detailed statistics of 5 key metrics that scMINER uses for filtration: `nUMI`, `nFeature`, `pctMito`, `pctSpikeIn`, `nCell`.
- **Detailed statistics per cell and gene**: it lists the detailed statistics of each gene and cell.
- **Filtration cutoffs by scMINER**: it provides the cutoffs estimated automatically by scMINER based on Median ± 3 * MAD (maximum absolute deviance), and the pseudo-filtration statistics on both genes and cells with these cutoffs.

### Filter the sparse eset object

From the quality control report generated above, we have got a better sense about the data quality and the cutoffs to use for filtration. scMINER provides two modes to do the filtration:
- **"auto"**: in this mode, scMINER will use the cutoffs estimated by Median ± 3*MAD (maximum absolute deviation). Based on our tests, in most cases, this mode works well with the matrix of both raw UMI counts and TPM (Transcripts Per Million) values.
- **"manual"**: in this mode, the users can manually specify the cutoffs, both low and high, of all 5 metrics: **nUMI**, **nFeature**, **pctMito**, **pctSpikeIn** for cells, and **nCell** for genes. No cells or features would be removed under the default cutoffs of each metrics.

Here are the codes to do the filtration:

```r
## to use the auto mode
pbmc14k_filtered.eset <- filterSparseEset(pbmc14k_raw.eset, mode = "auto")
#> Checking the availability of the 5 metrics ('nCell', 'nUMI', 'nFeature', 'pctMito', 'pctSpikeIn') used for eset filtration ...
#> Checking passed! All 5 metrics are available.
#> Filtration is done!
#> Filtration Summary:
#> 	8846/17986 genes passed!
#> 	13605/14000 cells passed!
#> 
#> For more details:
#> 	Gene filtration statistics:
#> 		Metrics		nCell
#> 		Cutoff_Low	70
#> 		Cutoff_High	Inf
#> 		Gene_total	17986
#> 		Gene_passed	8846(49.18%)
#> 		Gene_failed	9140(50.82%)
#> 
#> 	Cell filtration statistics:
#> 		Metrics		nUMI		nFeature	pctMito		pctSpikeIn	Combined
#> 		Cutoff_Low	458		221		0		0		NA
#> 		Cutoff_High	3694		Inf		0.0408		0.0000		NA
#> 		Cell_total	14000		14000		14000		14000		14000
#> 		Cell_passed	13826(98.76%)	14000(100.00%)	13778(98.41%)	14000(100.00%)	13605(97.18%)
#> 		Cell_failed	174(1.24%)	0(0.00%)	222(1.59%)	0(0.00%)	395(2.82%)
```


```r
## to use the manual mode
pbmc14k_filtered.eset <- filterSparseEset(pbmc14k_raw.eset, mode = "manual", gene.nCell_min = 10,
    cell.nUMI_min = 500, cell.nUMI_max = 6500, cell.nFeature_min = 200, cell.nFeature_max = 2500,
    cell.pctMito_max = 0.1)
```

The function, `filterSparseEset()`, generates a summary table and print it to the console. The user can refer to it and adjust the cutoffs accordingly.


```r
pbmc14k_filtered.eset
#> SparseExpressionSet (storageMode: environment)
#> assayData: 8846 features, 13605 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: CACTTTGACGCAAT GTTACGGAAACGAA ... ACGTGCCTTAAAGG (13605
#>     total)
#>   varLabels: cell_barcode true_label ... CellID (7 total)
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: LINC00115 NOC2L ... SRSF10.1 (8846 total)
#>   fvarLabels: GeneSymbol nCell
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation: .//PBMC14k
saveRDS(pbmc14k_filtered.eset, file = "/Users/qpan/PBMC14k/DATA/pbmc14k_filtered.rds")
```

## Normalize the SparseExpressionSet object

After the filtration of low-quality genes and cells, the eset object is ready for normalization. We recommend to use "log2CPM" method for normalization: the raw counts in each cell are normalized to a library size of 1 million, followed by log2 transformation.

```r
pbmc14k_log2cpm.eset <- normalizeSparseEset(pbmc14k_filtered.eset, scale_factor = 1e+06, log_base = 2,
    log_pseudoCount = 1)
#> Done! The data matrix of eset has been normalized and log-transformed!
#> The returned eset contains: 8846 genes, 13605 cells.
dim(pbmc14k_log2cpm.eset)
#> Features  Samples 
#>     8846    13605
saveRDS(pbmc14k_log2cpm.eset, file = "/Users/qpan/PBMC14k/DATA/pbmc14k_log2cpm.rds")
```

This normalized and log-transformed SparseExpresionSet object can be directly used for Mutual Information-based clustering, network inference and other downstream analysis. The preprocessing section is done!


# Mutual Information-based Clustering Analysis (MICA)

MICA (Mutual Information-based Clustering Analysis) is a clustering tool designed for single cell genomics data. Most existing single-cell clustering algorithms select the highly variable genes first and then perform principal component analysis (PCA) dimension reduction followed by graph-based or consensus k-means clustering. The selection of top variable features improves the clustering speed but is arbitrary and may lose the information that can distinguish close cell states. And the linear-transformation of PCA and co-expression analysis using linear Pearson or Spearman correlations may not capture the nonlinear cell-cell distance and gene-gene correlations. To address these concerns, MICA uses all high-quality features and employs the mutual information, which can catch both linear and non-linear correlations, to measure cell-cell similarity for unsupervised clustering analysis. A typical MICA run consists of 4 steps:
- Mutual information estimation for cell-cell distance quantification
- Dimension reduction on the non-linear mutual informaion-based distance space
- Clustering on dimension-reduced spaces
- Clustering visualization using UMAP or t-SNE method

As a component of scMINER framework, MICA works seamlessly with the SparseExpressionSet object. The input for MICA can be easily generated from the SparseExpressionSet object by `generateMICAinput()`, and the output of MICA, the clustering results, can be effortlessly visualized by `MICAplot()` and integrated into SparseExpressionSet object by `addMICAoutput()`.


## Generate the MICA input file

The standard input of MICA is **a normalized and log-transformed gene expression matrix**. Given the fact that the mutual information estimation of big data set could be compute-intensive, we developed the MICA using Python framework to take it's strengths in speed and memory consumption. For a better user experience in the interoperability between Python and R, scMINER saves the normalized and log-transformed gene expression matrix into a file that can be directly read by MICA. MICA accepts `.h5ad` or `.txt` format as the input file, which can be easily generated by embedded function `generateMICAinput()`:

```r
generateMICAinput(input_eset = pbmc14k_log2cpm.eset, output_file = "/Users/qpan/PBMC14k/MICA/micaInput.txt",
    overwrite = FALSE)
#> Writing MICA input to: /Users/qpan/PBMC14k/MICA/micaInput.txt 
#> 
#> For dataset with more than 5k cells, MICA GE mode is recommended.
#> Suggested command line is:
#> 
#>  mica ge -i /Users/qpan/PBMC14k/MICA/micaInput.txt -o project_space/MICA/mica_output_dir -ar 3.0 -ss 0.2 -nw 1 -nnm 80 
#> 
#>     Where options represent:
#>         -res:   the number of communities (default: 1.822)
#>         -minr:  the minimum number of communities to sweep (default: 1.822)
#>         -maxr:  the maximum number of communities to sweep (default: 1.822)
#>         -ss:    the step size to sweep number of communities from the mimimum to the maximum (default: 1)
#>         -nw:    the number of workers to run in parallel (default: 1, suggested: 25)
#>         -annef: the ef value of hnsw
#>         -annm:  the M value of hnsw
#>     Use mica ge -h to find more options of MICA GE mode.
```

To use the .h5ad format, run the codes below.

```r
generateMICAinput(input_eset = pbmc14k_log2cpm.eset, output_file = "/Users/qpan/PBMC14k/MICA/micaInput.h5ad",
    overwrite = FALSE)
```

In addition to generate the standard MICA input file, this function, `generateMICAinput()`, also prints the recommended commands of running MICA to the console. 

## Run the MICA

MICA features two different modes named by their different dimension reduction techniques:
- Multi-Dimensional Scaling (MDS) mode: this mode is more accurate and robust for small datasets (less than 5,000 cells, be default) due to its global dimension reduction nature;
- Graph Embedding (GE) mode: this mode works better with large datasets (more than 5,000 cells, by default) using a graph embedding approach to explore distant neighbor cells.

In this case, since there are 13,605 cells, we will use the `MICA GE` mode for the clustering:

```shell
mica ge -i /Users/qpan/PBMC14k/MICA/micaInput.txt -o /Users/qpan/PBMC14k/MICA/micaOutput -minr 0.1, -maxr 9.0 -ss 0.05 -nw 40
```

## Visualize and add the MICA output

MICA generates several files and save all of them in the output directory specified by the user with `-o` argument. The core, and only, output file we need for subsequent analysis is the clustering label file named in the format of `ProjectName_clustering_VisualizeMethod_euclidean_NumberOfDimensions_Resolution.txt`. In this case, since we used a range of resolutions, there are several clustering label files generated, one for each resolution. Based on the knowledge about PBMC14k dataset, we compared the results of different resolutions and picked `clustering_UMAP_euclidean_20_2.05.txt` for subsequent analysis.


```r
data("pbmc14k_MICAoutput")
head(pbmc14k_MICAoutput)
#>               ID        X        Y label
#> 1 CACTTTGACGCAAT 14.91650 13.04096     6
#> 2 GTTACGGAAACGAA 14.57031 10.27093     6
#> 3 CACTTATGAGTCGT 14.28869 13.61674     6
#> 4 GCATGTGATTCTGT 14.12546 13.36319     6
#> 5 TAGAATACGTATCG 14.91227 11.19407     6
#> 6 CAAGAAGACCCTCA 15.34154 12.25821     6
```
As shown above, the clustering label file contains four columns:
- `ID`: cell barcodes;
- `X`: coordinates of UMAP_1 or tSNE_1;
- `Y`: coordinates of UMAP_2 or tSNE_2;
- `label`: labels of predicted clusters.

The clustering result can be easily easily added to the SparseExpressionSet object by `addMICAoutput()` and visualized by `MICAplot()` function:

```r
pbmc14k_log2cpm.eset <- addMICAoutput(pbmc14k_log2cpm.eset, mica_output_file = system.file("extdata",
    "pbmc14k_projectSpace/MICA/clustering_UMAP_euclidean_20_2.05.txt", package = "scMINER"), visual_method = "umap")
```


# Cell type annotation

The clustering analysis of single cell RNA-seq data is challenging, while the interpretation of the clusters could be more difficult, since the clustering results could be noisy. Different clusters could be of different cell types, same cell type in different conditions, or same cell type in same condition but with different QC statistics (e.g. doublets, different gene coverage). This makes the cell type annotation a common bottleneck of the whole scRNA-seq data analysis.

Currently, there are two types of strategies to annotate the clusters: **supervised** and **unsupervised**. The **supervised** methods use a list of known markers of potential cell types curated from somoe existing studies of the same/similar contexts. While in contrast, the **unsupervised** methods are usually based on the differentially expressed genes. scMINER provide several useful functions to support both types of strategies.

## When a mark list is available

In most cases, some existing studies on the same/similar tissue type are available, and from these existing studies, we can figure out a list of cell type to expect to see and curate a list of markers for each of them. In this case, we know the 7 cell types involved in the dataset, and based on some other PBMCs studies, we manually curated a marker list and generated a signature table as shown below:

```r
data("pbmc14k_signatureTable")
head(pbmc14k_signatureTable)
#>   signature_name signature_feature weight
#> 1       Monocyte              CD14      1
#> 2       Monocyte               LYZ      1
#> 3       Monocyte            S100A8      1
#> 4       Monocyte            S100A9      1
#> 5       Monocyte           S100A12      1
#> 6             NK            FCGR3A      1
```
The signature table contains three columns:
- `signature_name`: name of cell types/signatures;
- `signature_feature`: markers genes/features of corresponding cell type/signature;
- `weight`: weight of corresponding maker/feature in corresponding cell type/signature. It ranges from -1 to 1, so both positive and negtive markers are supoorted.

With this signature table, scMINER can estimate a signature score, which is mathematically `the weighted mean of the expression of marker genes involved`, for each cell type from each of the clusters:


```r
## Violin plot of marker genes across clusters
draw_bubbleplot(input_eset = pbmc14k_log2cpm.eset, signature_table = pbmc14k_signatureTable, group_by = "clusterID")
#> 31 features of 7 signatures were found in the input eset and will be used in calculation.
```

<img src="/private/var/folders/61/192x8zbn7slbdn2hvw_11lpw0000gp/T/RtmpXR975a/preview-10879f530be1.dir/pbmc14k_tutorial_files/figure-html/signature-bubble-1.png" style="display: block; margin: auto;" />

In the bubble plot above, the color of the bubbles is proportional to the mean of signature score, and the size of the bubbles is proportional to the percentage of cells with higher signature score than mean. The cell type of each cluster is clear, except the cluster 7, which shows equally-high signature score of both CD4+ TCM and CD4+ Reg and higher percentage of CD4+ TCM cells.

scMINER also provides a variety of functions to visualize the selected features:

```r
## For the demonstration purposes, we picked two well known markers for each of the 7 known
## cell types, plus 'CD3D' and 'CD4'.
genes_of_interest <- c("CD14", "LYZ", "GZMB", "NKG7", "CD19", "MS4A1", "CD8A", "CD8B", "SELL", "CCR7",
    "IL2RA", "FOXP3", "IL7R", "S100A4", "CD3D", "CD4")
```


```r
## Violin plot of marker genes across clusters
feature_vlnplot(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, group_by = "clusterID",
    ncol = 4)
```

<img src="/private/var/folders/61/192x8zbn7slbdn2hvw_11lpw0000gp/T/RtmpXR975a/preview-10879f530be1.dir/pbmc14k_tutorial_files/figure-html/featurePlot-vln-1.png" style="display: block; margin: auto;" />


```r
## Box plot of marker genes across clusters
feature_boxplot(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, group_by = "clusterID",
    ncol = 4)
```

<img src="/private/var/folders/61/192x8zbn7slbdn2hvw_11lpw0000gp/T/RtmpXR975a/preview-10879f530be1.dir/pbmc14k_tutorial_files/figure-html/featurePlot-box-1.png" style="display: block; margin: auto;" />


```r
## UMAP scatter plot of marker genes
feature_scatterplot(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, ncol = 4, location_x = "UMAP_1",
    location_y = "UMAP_2", point.size = 0.5)
```

<img src="/private/var/folders/61/192x8zbn7slbdn2hvw_11lpw0000gp/T/RtmpXR975a/preview-10879f530be1.dir/pbmc14k_tutorial_files/figure-html/featurePlot-scatter-1.png" style="display: block; margin: auto;" />

```
#> TableGrob (4 x 4) "arrange": 16 grobs
#>     z     cells    name           grob
#> 1   1 (1-1,1-1) arrange gtable[layout]
#> 2   2 (1-1,2-2) arrange gtable[layout]
#> 3   3 (1-1,3-3) arrange gtable[layout]
#> 4   4 (1-1,4-4) arrange gtable[layout]
#> 5   5 (2-2,1-1) arrange gtable[layout]
#> 6   6 (2-2,2-2) arrange gtable[layout]
#> 7   7 (2-2,3-3) arrange gtable[layout]
#> 8   8 (2-2,4-4) arrange gtable[layout]
#> 9   9 (3-3,1-1) arrange gtable[layout]
#> 10 10 (3-3,2-2) arrange gtable[layout]
#> 11 11 (3-3,3-3) arrange gtable[layout]
#> 12 12 (3-3,4-4) arrange gtable[layout]
#> 13 13 (4-4,1-1) arrange gtable[layout]
#> 14 14 (4-4,2-2) arrange gtable[layout]
#> 15 15 (4-4,3-3) arrange gtable[layout]
#> 16 16 (4-4,4-4) arrange gtable[layout]
```


```r
## Bubble plot of marker genes across clusters
feature_bubbleplot(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, group_by = "clusterID",
    xlabel.angel = 45)
```

<img src="/private/var/folders/61/192x8zbn7slbdn2hvw_11lpw0000gp/T/RtmpXR975a/preview-10879f530be1.dir/pbmc14k_tutorial_files/figure-html/featurePlot-bubble-1.png" style="display: block; margin: auto;" />


```r
## Heatmap of marker genes across clusters
feature_heatmap(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, group_by = "clusterID",
    scale_method = "none", annotation_columns = c("true_label"))
```

<img src="/private/var/folders/61/192x8zbn7slbdn2hvw_11lpw0000gp/T/RtmpXR975a/preview-10879f530be1.dir/pbmc14k_tutorial_files/figure-html/featurePlot-heatmap-1.png" style="display: block; margin: auto;" />


## When a mark list is not available

The existing studies on same/similar context are not always available. And, there is a major concern about the reliability of the reference studies which largely depends on the expertise of the original authors who defined the markers and assigned the cell types. So, we always encourage the users to try the unsupervised methods as well, which can serve as a cross-validation method.

scMINER provides a function, `getDE()`, to perform the differential expression analysis and identify the markers of each cluster. The `getDE()` function supports three different methods to perform the differential expression analysis, `limma`, `wilcoxon` and `t.test`. And it allows the user to define the groups to compare in a flexible way:


```r
## 1. To perform differential expression analysis in a 1-vs-rest manner for all groups
de_res1 <- getDE(input_eset = pbmc14k_log2cpm.eset, group_by = "clusterID", use_method = "limma")
#> 7 groups were found in group_by column [ clusterID ].
#> Since no group was specified, the differential analysis will be conducted among all groups in the group_by column [ clusterID ] in the 1-vs-rest manner.
#> 	 1 / 7 : group 1 ( 1 ) vs the rest...
#> 	 2505 cells were found for g1.
#> 	 11100 cells were found for g0.
#> 	 2 / 7 : group 1 ( 2 ) vs the rest...
#> 	 2022 cells were found for g1.
#> 	 11583 cells were found for g0.
#> 	 3 / 7 : group 1 ( 3 ) vs the rest...
#> 	 2014 cells were found for g1.
#> 	 11591 cells were found for g0.
#> 	 4 / 7 : group 1 ( 4 ) vs the rest...
#> 	 1918 cells were found for g1.
#> 	 11687 cells were found for g0.
#> 	 5 / 7 : group 1 ( 5 ) vs the rest...
#> 	 1912 cells were found for g1.
#> 	 11693 cells were found for g0.
#> 	 6 / 7 : group 1 ( 6 ) vs the rest...
#> 	 1786 cells were found for g1.
#> 	 11819 cells were found for g0.
#> 	 7 / 7 : group 1 ( 7 ) vs the rest...
#> 	 1448 cells were found for g1.
#> 	 12157 cells were found for g0.
head(de_res1)
#>      feature g1_tag      g0_tag    g1_avg   g0_avg    g1_pct    g0_pct   log2FC
#> 1251    CD3E      1 2,3,4,5,6,7  8.354660 3.874230 0.7920160 0.3819820 4.480430
#> 3820    LDHB      1 2,3,4,5,6,7  9.555670 5.614992 0.8806387 0.5458559 3.940678
#> 7765  TMEM66      1 2,3,4,5,6,7  8.604421 5.041570 0.8103792 0.5051351 3.562851
#> 1250    CD3D      1 2,3,4,5,6,7  7.281998 4.097082 0.6990020 0.3965766 3.184916
#> 1235    CD27      1 2,3,4,5,6,7  5.566280 2.428199 0.5481038 0.2482883 3.138081
#> 3992     LTB      1 2,3,4,5,6,7 10.436707 7.315803 0.9141717 0.6430631 3.120905
#>               Pval           FDR   Zscore
#> 1251 2.225074e-308  0.000000e+00 37.53784
#> 3820 4.919041e-274 6.216262e-271 35.37012
#> 7765 1.509154e-228 9.535698e-226 32.27621
#> 1250 3.193659e-174 1.130044e-171 28.13981
#> 1235 7.997857e-219 4.716603e-216 31.57555
#> 3992 5.621459e-159 1.841756e-156 26.86509
```


```r
## 2. To perform differential expression analysis in a 1-vs-rest manner for one specific group
de_res2 <- getDE(input_eset = pbmc14k_log2cpm.eset, group_by = "clusterID", g1 = c("1"), use_method = "limma")

## 3. To perform differential expression analysis in a rest-vs-1 manner for one specific group
de_res3 <- getDE(input_eset = pbmc14k_log2cpm.eset, group_by = "clusterID", g0 = c("1"), use_method = "limma")

## 4. To perform differential expression analysis in a 1-vs-1 manner for any two groups
de_res4 <- getDE(input_eset = pbmc14k_log2cpm.eset, group_by = "clusterID", g1 = c("1"), g0 = c("3"),
    use_method = "limma")
```

scMINER also provides another function, `getTopFeatures()`, to easily extract the group-specific markers from the differential expression result:

```r
cluster_markers <- getTopFeatures(input_table = de_res1, number = 10, group_by = "g1_tag", sort_by = "log2FC",
    sort_decreasing = TRUE)
dim(cluster_markers)
#> [1] 16 11
head(cluster_markers)
#>      feature g1_tag      g0_tag    g1_avg   g0_avg    g1_pct    g0_pct   log2FC
#> 1251    CD3E      1 2,3,4,5,6,7  8.354660 3.874230 0.7920160 0.3819820 4.480430
#> 3820    LDHB      1 2,3,4,5,6,7  9.555670 5.614992 0.8806387 0.5458559 3.940678
#> 7765  TMEM66      1 2,3,4,5,6,7  8.604421 5.041570 0.8103792 0.5051351 3.562851
#> 1250    CD3D      1 2,3,4,5,6,7  7.281998 4.097082 0.6990020 0.3965766 3.184916
#> 1235    CD27      1 2,3,4,5,6,7  5.566280 2.428199 0.5481038 0.2482883 3.138081
#> 3992     LTB      1 2,3,4,5,6,7 10.436707 7.315803 0.9141717 0.6430631 3.120905
#>               Pval           FDR   Zscore
#> 1251 2.225074e-308  0.000000e+00 37.53784
#> 3820 4.919041e-274 6.216262e-271 35.37012
#> 7765 1.509154e-228 9.535698e-226 32.27621
#> 1250 3.193659e-174 1.130044e-171 28.13981
#> 1235 7.997857e-219 4.716603e-216 31.57555
#> 3992 5.621459e-159 1.841756e-156 26.86509
```
## Add cell type annotations to SparseExpressionSet object

Based on the supervised and unsupervised methods, we have annotated the cell types for each cluster. To add the cell type annotation information into the sparse eset object:

```r
celltype_map <- c(`1` = "CD4TN", `2` = "CD4TCM", `3` = "CD8TN", `4` = "NK", `5` = "B", `6` = "Monocyte",
    `7` = "CD4Treg")
pbmc14k_log2cpm.eset$cell_type <- as.character(celltype_map[pbmc14k_log2cpm.eset$clusterID])
```

The `draw_barplot()` function can visualize the cell composition of self-defined groups. We can use it to show the purity of MICA clusters:

```r
## Violin plot of marker genes across clusters
draw_barplot(input_eset = pbmc14k_log2cpm.eset, group_by = "cell_type", color_by = "true_label",
    xlabel.angel = 45)
```

<img src="/private/var/folders/61/192x8zbn7slbdn2hvw_11lpw0000gp/T/RtmpXR975a/preview-10879f530be1.dir/pbmc14k_tutorial_files/figure-html/group-barplot-1.png" style="display: block; margin: auto;" />

Don't forget to save the SparseExpressionSet object after the cell type annotation added.

```r
saveRDS(pbmc14k_log2cpm.eset, file = "/Users/qpan/PBMC14k/DATA/pbmc14k_log2CPM_annotated.rds")
```


# Mutual Information-based Network Inference by SJARACNe

[SJARACNe](https://academic.oup.com/bioinformatics/article/35/12/2165/5156064) is a scalable software tool for gene network reverse engineering from big data. As an improved implementation of the ARACNe, SJARACNe achieves a dramatic improvement in computational performance in both time and memory usage and implements new features while preserving the network inference accuracy of the original algorithm.

Similar to MICA, SJARACNe is also a component of scMINER framework, and can work seamlessly with the SparseExpressionSet object. The input for MICA can be easily generated from the SparseExpressionSet object by `generateSJARACNeInput()`, and the output of SJARACNe, the networks, can be effortlessly assessed by `drawNetworkQC()` and directly taken for driver activity estimation by `getActivity_individual()` and `getActivity_inBatch()`.


## Generate the SJARACNe input files

The network inference should be performed in a cluster- or cell type-specific basis. For each group (one cluster or cell type), the `generateSJARACNeInput()` function will create a folder named by the label of this group, and generated the standard input files in it:
- a "**`.exp.txt`**" file: a tab-separated genes/transcripts/proteins by cells/samples expression matrix with the first two columns being ID and symbol.
- a "**`TF`**" folder containing a "**`.tf.txt`**" file: a list of significant gene/transcript/protein IDs of TF drivers.
- a "**`SIG`**" folder containing a "**`.sig.txt`**" file: a list of significant gene/transcript/protein IDs of SIG drivers.
- a bash script (**`runSJARACNe.sh`**) to run SJARACNe. Further modification is needed to run it.
- a json file (**`config_cwlexec.json`**) containing parameters to run SJARACNe.

Usually, the ground truth of cell types is not available. Then the cluster labels, or cell type annotations of the clusters, can be used for grouping in network rewiring, since it's expected that cells with same cluster label/annotated cell type are of similar gene expression profiles. To generate from annotated cell types, you can run:

```r
generateSJARACNeInput(input_eset = pbmc14k_log2cpm.eset, group_name = "cell_type", sjaracne_dir = "/Users/qpan/PBMC14k/SJARACNe/bycelltype",
    species_type = "hg", driver_type = "TF_SIG")
```

In this case, we know the true label of cell types, so we use them to define the groups for network rewiring. **NOTE: scMINER creates a folder for each of groups using the group labels. Any illegal characters in the group labels may cause issues in subsequenct analysis. To avoid this, scMINER only accept letters(A-Za-z), numbers(0-9), underscores('_') and periods('.').**

```r
truelabel_map <- c(`CD14+ Monocyte` = "Monocyte", `CD19+ B` = "B", `CD4+/CD25 T Reg` = "CD4Treg",
    `CD4+/CD45RA+/CD25- Naive T` = "CD4TN", `CD4+/CD45RO+ Memory` = "CD4TCM", `CD56+ NK` = "NK",
    `CD8+/CD45RA+ Naive Cytotoxic` = "CD8TN")
pbmc14k_log2cpm.eset$true_label.new <- as.character(truelabel_map[pbmc14k_log2cpm.eset$true_label])
generateSJARACNeInput(input_eset = pbmc14k_log2cpm.eset, group_name = "true_label.new", sjaracne_dir = "/Users/qpan/PBMC14k/SJARACNe",
    species_type = "hg", driver_type = "TF_SIG", downSample_N = NULL)
#> No non-word characters (any character other than [a-zA-Z0-9_] ) were found in the group_name elements. group_name check passed...
#> Extracting TF drivers from the reference list...
#> A TF driver list of 2008 genes was generated!
#> Extracting SIG drivers from the reference list...
#> A SIG driver list of 9704 genes was generated!
#> 7 groups were found from the input eset.
#> 1/7: Generating SJARACNe inputs for group: Monocyte
#> 	Creating the output dirctory...
#> 	1837 cells were found in this group...
#> 	The metacell analysis by was skipped, since the superCell_N is null.
#> 	The downsampling was skipped, since the downSample_N is null.
#> 	Writing gene expression matrix into /Users/qpan/PBMC14k/SJARACNe/Monocyte/Monocyte.8573_1837.exp.txt...
#> 	Writing TF driver list into /Users/qpan/PBMC14k/SJARACNe/Monocyte/TF/Monocyte.832_1837.tf.txt...
#> 	Writing SIG driver list into /Users/qpan/PBMC14k/SJARACNe/Monocyte/SIG/Monocyte.4191_1837.sig.txt...
#> 	Writing command lines to run SJARACNe to /Users/qpan/PBMC14k/SJARACNe/Monocyte/runSJARACNe.sh...
#> 	Done.
#> 2/7: Generating SJARACNe inputs for group: B
#> 	Creating the output dirctory...
#> 	1902 cells were found in this group...
#> 	The metacell analysis by was skipped, since the superCell_N is null.
#> 	The downsampling was skipped, since the downSample_N is null.
#> 	Writing gene expression matrix into /Users/qpan/PBMC14k/SJARACNe/B/B.8572_1902.exp.txt...
#> 	Writing TF driver list into /Users/qpan/PBMC14k/SJARACNe/B/TF/B.835_1902.tf.txt...
#> 	Writing SIG driver list into /Users/qpan/PBMC14k/SJARACNe/B/SIG/B.4148_1902.sig.txt...
#> 	Writing command lines to run SJARACNe to /Users/qpan/PBMC14k/SJARACNe/B/runSJARACNe.sh...
#> 	Done.
#> 3/7: Generating SJARACNe inputs for group: CD4Treg
#> 	Creating the output dirctory...
#> 	1985 cells were found in this group...
#> 	The metacell analysis by was skipped, since the superCell_N is null.
#> 	The downsampling was skipped, since the downSample_N is null.
#> 	Writing gene expression matrix into /Users/qpan/PBMC14k/SJARACNe/CD4Treg/CD4Treg.8659_1985.exp.txt...
#> 	Writing TF driver list into /Users/qpan/PBMC14k/SJARACNe/CD4Treg/TF/CD4Treg.838_1985.tf.txt...
#> 	Writing SIG driver list into /Users/qpan/PBMC14k/SJARACNe/CD4Treg/SIG/CD4Treg.4214_1985.sig.txt...
#> 	Writing command lines to run SJARACNe to /Users/qpan/PBMC14k/SJARACNe/CD4Treg/runSJARACNe.sh...
#> 	Done.
#> 4/7: Generating SJARACNe inputs for group: CD4TN
#> 	Creating the output dirctory...
#> 	1994 cells were found in this group...
#> 	The metacell analysis by was skipped, since the superCell_N is null.
#> 	The downsampling was skipped, since the downSample_N is null.
#> 	Writing gene expression matrix into /Users/qpan/PBMC14k/SJARACNe/CD4TN/CD4TN.8612_1994.exp.txt...
#> 	Writing TF driver list into /Users/qpan/PBMC14k/SJARACNe/CD4TN/TF/CD4TN.831_1994.tf.txt...
#> 	Writing SIG driver list into /Users/qpan/PBMC14k/SJARACNe/CD4TN/SIG/CD4TN.4180_1994.sig.txt...
#> 	Writing command lines to run SJARACNe to /Users/qpan/PBMC14k/SJARACNe/CD4TN/runSJARACNe.sh...
#> 	Done.
#> 5/7: Generating SJARACNe inputs for group: CD4TCM
#> 	Creating the output dirctory...
#> 	1967 cells were found in this group...
#> 	The metacell analysis by was skipped, since the superCell_N is null.
#> 	The downsampling was skipped, since the downSample_N is null.
#> 	Writing gene expression matrix into /Users/qpan/PBMC14k/SJARACNe/CD4TCM/CD4TCM.8660_1967.exp.txt...
#> 	Writing TF driver list into /Users/qpan/PBMC14k/SJARACNe/CD4TCM/TF/CD4TCM.838_1967.tf.txt...
#> 	Writing SIG driver list into /Users/qpan/PBMC14k/SJARACNe/CD4TCM/SIG/CD4TCM.4209_1967.sig.txt...
#> 	Writing command lines to run SJARACNe to /Users/qpan/PBMC14k/SJARACNe/CD4TCM/runSJARACNe.sh...
#> 	Done.
#> 6/7: Generating SJARACNe inputs for group: NK
#> 	Creating the output dirctory...
#> 	1936 cells were found in this group...
#> 	The metacell analysis by was skipped, since the superCell_N is null.
#> 	The downsampling was skipped, since the downSample_N is null.
#> 	Writing gene expression matrix into /Users/qpan/PBMC14k/SJARACNe/NK/NK.8654_1936.exp.txt...
#> 	Writing TF driver list into /Users/qpan/PBMC14k/SJARACNe/NK/TF/NK.841_1936.tf.txt...
#> 	Writing SIG driver list into /Users/qpan/PBMC14k/SJARACNe/NK/SIG/NK.4202_1936.sig.txt...
#> 	Writing command lines to run SJARACNe to /Users/qpan/PBMC14k/SJARACNe/NK/runSJARACNe.sh...
#> 	Done.
#> 7/7: Generating SJARACNe inputs for group: CD8TN
#> 	Creating the output dirctory...
#> 	1984 cells were found in this group...
#> 	The metacell analysis by was skipped, since the superCell_N is null.
#> 	The downsampling was skipped, since the downSample_N is null.
#> 	Writing gene expression matrix into /Users/qpan/PBMC14k/SJARACNe/CD8TN/CD8TN.8568_1984.exp.txt...
#> 	Writing TF driver list into /Users/qpan/PBMC14k/SJARACNe/CD8TN/TF/CD8TN.825_1984.tf.txt...
#> 	Writing SIG driver list into /Users/qpan/PBMC14k/SJARACNe/CD8TN/SIG/CD8TN.4145_1984.sig.txt...
#> 	Writing command lines to run SJARACNe to /Users/qpan/PBMC14k/SJARACNe/CD8TN/runSJARACNe.sh...
#> 	Done.
#> SparseExpressionSet (storageMode: environment)
#> assayData: 8846 features, 13605 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: CACTTTGACGCAAT GTTACGGAAACGAA ... ACGTGCCTTAAAGG (13605
#>     total)
#>   varLabels: cell_barcode true_label ... true_label.new (12 total)
#>   varMetadata: labelDescription
#> featureData
#>   featureNames: LINC00115 NOC2L ... SRSF10.1 (8846 total)
#>   fvarLabels: GeneSymbol nCell
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> Annotation:
```

## Run the SJARACNe

By default, the `generateSJARACNeInput()` function also generates a `runSJARACNe.sh` file in the folder of each group. This file much be modified before you can run it:
- **removed unneeded lines**: There are usually 4 lines in this file: the lines starting with "sjaracne lsf" are the command lines to run on IBM LSF cluster, while the lines starting with "sjaracne local" are the command lines runing on a single machine (Linux/OSX). Please select the lines based on your situation and remove the others.
- **-n**: number of bootstrap networks to generate. Default: 100.
- **-pc**: p value threshold to select edges in building consensus network. Default: e-2 for single-cell data, e-3 for meta-cell data, and e-5 for bulk sample data.
Please use "sjaracne lsf -h" or "sjaracne local -h" to check more details of arguments available in SJARACNe.

There is another file, `config_cwlexec.json`, available in the folder. It contains the information (e.g. memory request for each step of SJARACNe run) used for LSF job submission. This file is only needed for LSF runs and the default values works well in most cases. If you are running SJARACNe on a big dataset, you may need to request more memory from it.

In this case, we use LSF to run the SJARACNe:

```shell
## let's use B cell as an example
# for TF
sjaracne lsf -e /Users/qpan/PBMC14k/SJARACNe/B/B.8572_1902.exp.txt -g /Users/qpan/PBMC14k/SJARACNe/B/TF/B.835_1902.tf.txt -o /Users/qpan/PBMC14k/SJARACNe/B/TF/bt100_pc001 -n 100 -pc 0.01 -j /Users/qpan/PBMC14k/SJARACNe/B/config_cwlexec.json
# for SIG
sjaracne lsf -e /Users/qpan/PBMC14k/SJARACNe/B/B.8572_1902.exp.txt -g /Users/qpan/PBMC14k/SJARACNe/B/SIG/B.4148_1902.sig.txt -o /Users/qpan/PBMC14k/SJARACNe/B/SIG/bt100_pc001 -n 100 -pc 0.01 -j /Users/qpan/PBMC14k/SJARACNe/B/config_cwlexec.json
```

We created a folder named "bt100_pc001" in both TF and SIG folders of each group, to save the networks generated under 100 bootstraps and 0.01 consensus p value.

To run SJARACNe on a local machine:

```shell
## let's use B cell as an example
# for TF
sjaracne local -e /Users/qpan/PBMC14k/SJARACNe/B/B.8572_1902.exp.txt -g /Users/qpan/PBMC14k/SJARACNe/B/TF/B.835_1902.tf.txt -o /Users/qpan/PBMC14k/SJARACNe/B/TF/bt100_pc001 -n 100 -pc 0.01
# for SIG
sjaracne lsf -e /Users/qpan/PBMC14k/SJARACNe/B/B.8572_1902.exp.txt -g /Users/qpan/PBMC14k/SJARACNe/B/SIG/B.4148_1902.sig.txt -o /Users/qpan/PBMC14k/SJARACNe/B/SIG/bt100_pc001 -n 100 -pc 0.01
```

## Assess the quality of networks

scMINER provides a function, `drawNetworkQC()`, to quickly assess the quality of networks in batch.


```r
stat_table <- drawNetworkQC(sjaracne_dir = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe",
    generate_html = FALSE)
#> The sjaracne_dir is specified for network QC: /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe .
#> 18 network files were found in the directory: /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe .
#> 	Processing file 1 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/SIG/bt100_pc001/sjaracne_workflow-df798096-8dee-4baf-8f70-891c689dc769/consensus_network_ncol_.txt ...
#> 	Processing file 2 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/TF/bt100_pc001/sjaracne_workflow-fb2a69b9-f98e-47ff-87a0-6d538822fc6e/consensus_network_ncol_.txt ...
#> 	Processing file 3 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4TCM/SIG/bt100_pc001/sjaracne_workflow-424f1068-13d1-4f0e-9c26-56acd9a2027c/consensus_network_ncol_.txt ...
#> 	Processing file 4 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4TCM/TF/bt100_pc001/sjaracne_workflow-52b3cdf5-5914-4c8c-a77a-05f17c755d83/consensus_network_ncol_.txt ...
#> 	Processing file 5 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4TN/SIG/bt100_pc001/sjaracne_workflow-7b5bb68e-1de5-4d0e-80ec-8d8aa037866f/consensus_network_ncol_.txt ...
#> 	Processing file 6 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4TN/TF/bt100_pc001/sjaracne_workflow-89716541-eb53-435c-8a45-bab63d6b5198/consensus_network_ncol_.txt ...
#> 	Processing file 7 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/SIG/bt100_pc001/sjaracne_workflow-2703654d-4235-4082-ae27-b76d6b124007/consensus_network_ncol_.txt ...
#> 	Processing file 8 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/TF/bt100_pc001/sjaracne_workflow-d6744c02-79df-4060-9178-50748b5bdda0/consensus_network_ncol_.txt ...
#> 	Processing file 9 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/SIG/bt100_pc001/sjaracne_workflow-afa3c623-36ad-455b-8f48-2bedeeccb09a/consensus_network_ncol_.txt ...
#> 	Processing file 10 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/TF/bt100_pc001/sjaracne_workflow-1264efea-a69b-41c0-b1ea-d2b0542f29cc/consensus_network_ncol_.txt ...
#> 	Processing file 11 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/SIG/bt100_pc001/sjaracne_workflow-09d281a6-e0f7-4e7b-a253-5724de4d3dde/consensus_network_ncol_.txt ...
#> 	Processing file 12 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/TF/bt100_pc001/sjaracne_workflow-f41fdae8-399a-480a-9e13-b03195cdf857/consensus_network_ncol_.txt ...
#> 	Processing file 13 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Network_superCell/gamma_50/PBMC20K/SIG/sjaracne_workflow-bc050977-276a-42f4-8ff7-fd00aa8cf213/consensus_network_ncol_.txt ...
#> 	Processing file 14 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Network_superCell/gamma_50/PBMC20K/TF/sjaracne_workflow-ff40c564-c158-4a4c-8551-fca9d793c464/consensus_network_ncol_.txt ...
#> 	Processing file 15 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Network_superCell/supercell_default/PBMC20K/SIG/sjaracne_workflow-396debe4-6df3-46e0-bc8b-d0331cf321d7/consensus_network_ncol_.txt ...
#> 	Processing file 16 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Network_superCell/supercell_default/PBMC20K/TF/sjaracne_workflow-039a2e97-286d-402c-bbf9-3be8d36aa14d/consensus_network_ncol_.txt ...
#> 	Processing file 17 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/SIG/bt100_pc001/sjaracne_workflow-195700c8-57dd-458c-af1b-fd9a1e75a652/consensus_network_ncol_.txt ...
#> 	Processing file 18 / 18 : /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/TF/bt100_pc001/sjaracne_workflow-831d81a2-802e-4744-9695-0ed97230935a/consensus_network_ncol_.txt ...
stat_table
#>                                        network_tag network_node network_edge
#> 1                                B.SIG.bt100_pc001         8572       391889
#> 2                                 B.TF.bt100_pc001         8572        95341
#> 3                           CD4TCM.SIG.bt100_pc001         8660       382153
#> 4                            CD4TCM.TF.bt100_pc001         8660        94319
#> 5                            CD4TN.SIG.bt100_pc001         8612       401658
#> 6                             CD4TN.TF.bt100_pc001         8612        95152
#> 7                          CD4Treg.SIG.bt100_pc001         8659       381306
#> 8                           CD4Treg.TF.bt100_pc001         8659        97596
#> 9                            CD8TN.SIG.bt100_pc001         8568       395750
#> 10                            CD8TN.TF.bt100_pc001         8568        94308
#> 11                        Monocyte.SIG.bt100_pc001         8573       400924
#> 12                         Monocyte.TF.bt100_pc001         8573        90193
#> 13          Network_superCell.gamma_50.PBMC20K.SIG         8838       358570
#> 14           Network_superCell.gamma_50.PBMC20K.TF         8842       260344
#> 15 Network_superCell.supercell_default.PBMC20K.SIG         8842       298867
#> 16  Network_superCell.supercell_default.PBMC20K.TF         8844       222732
#> 17                              NK.SIG.bt100_pc001         8654       342785
#> 18                               NK.TF.bt100_pc001         8653        99335
#>    driver_count targetSize_mean targetSize_median targetSize_minimum
#> 1          4148        94.47662              94.0                 33
#> 2           835       114.18084              96.0                 64
#> 3          4209        90.79425              91.0                 31
#> 4           838       112.55251              95.5                 60
#> 5          4180        96.09043              95.0                 43
#> 6           831       114.50301              99.0                 64
#> 7          4214        90.48552              90.0                 29
#> 8           838       116.46301              96.0                 63
#> 9          4145        95.47648              95.0                 38
#> 10          825       114.31273             101.0                 68
#> 11         4191        95.66309              95.0                 41
#> 12          832       108.40505              98.0                 70
#> 13         4309        83.21420              73.0                  1
#> 14          850       306.28706             260.5                  8
#> 15         4309        69.35878              64.0                  1
#> 16          850       262.03765             229.0                 21
#> 17         4202        81.57663              82.0                 15
#> 18          841       118.11534              92.0                 64
#>    targetSize_maximum
#> 1                 396
#> 2                 913
#> 3                 281
#> 4                 689
#> 5                 303
#> 6                 743
#> 7                 401
#> 8                 846
#> 9                 243
#> 10                691
#> 11                255
#> 12                580
#> 13                511
#> 14               1241
#> 15                346
#> 16                955
#> 17                343
#> 18               1001
#>                                                                                                                                                                                                                         network_path
#> 1                                /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/SIG/bt100_pc001/sjaracne_workflow-df798096-8dee-4baf-8f70-891c689dc769/consensus_network_ncol_.txt
#> 2                                 /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/TF/bt100_pc001/sjaracne_workflow-fb2a69b9-f98e-47ff-87a0-6d538822fc6e/consensus_network_ncol_.txt
#> 3                           /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4TCM/SIG/bt100_pc001/sjaracne_workflow-424f1068-13d1-4f0e-9c26-56acd9a2027c/consensus_network_ncol_.txt
#> 4                            /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4TCM/TF/bt100_pc001/sjaracne_workflow-52b3cdf5-5914-4c8c-a77a-05f17c755d83/consensus_network_ncol_.txt
#> 5                            /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4TN/SIG/bt100_pc001/sjaracne_workflow-7b5bb68e-1de5-4d0e-80ec-8d8aa037866f/consensus_network_ncol_.txt
#> 6                             /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4TN/TF/bt100_pc001/sjaracne_workflow-89716541-eb53-435c-8a45-bab63d6b5198/consensus_network_ncol_.txt
#> 7                          /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/SIG/bt100_pc001/sjaracne_workflow-2703654d-4235-4082-ae27-b76d6b124007/consensus_network_ncol_.txt
#> 8                           /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD4Treg/TF/bt100_pc001/sjaracne_workflow-d6744c02-79df-4060-9178-50748b5bdda0/consensus_network_ncol_.txt
#> 9                            /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/SIG/bt100_pc001/sjaracne_workflow-afa3c623-36ad-455b-8f48-2bedeeccb09a/consensus_network_ncol_.txt
#> 10                            /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/CD8TN/TF/bt100_pc001/sjaracne_workflow-1264efea-a69b-41c0-b1ea-d2b0542f29cc/consensus_network_ncol_.txt
#> 11                        /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/SIG/bt100_pc001/sjaracne_workflow-09d281a6-e0f7-4e7b-a253-5724de4d3dde/consensus_network_ncol_.txt
#> 12                         /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Monocyte/TF/bt100_pc001/sjaracne_workflow-f41fdae8-399a-480a-9e13-b03195cdf857/consensus_network_ncol_.txt
#> 13          /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Network_superCell/gamma_50/PBMC20K/SIG/sjaracne_workflow-bc050977-276a-42f4-8ff7-fd00aa8cf213/consensus_network_ncol_.txt
#> 14           /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Network_superCell/gamma_50/PBMC20K/TF/sjaracne_workflow-ff40c564-c158-4a4c-8551-fca9d793c464/consensus_network_ncol_.txt
#> 15 /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Network_superCell/supercell_default/PBMC20K/SIG/sjaracne_workflow-396debe4-6df3-46e0-bc8b-d0331cf321d7/consensus_network_ncol_.txt
#> 16  /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/Network_superCell/supercell_default/PBMC20K/TF/sjaracne_workflow-039a2e97-286d-402c-bbf9-3be8d36aa14d/consensus_network_ncol_.txt
#> 17                              /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/SIG/bt100_pc001/sjaracne_workflow-195700c8-57dd-458c-af1b-fd9a1e75a652/consensus_network_ncol_.txt
#> 18                               /Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/NK/TF/bt100_pc001/sjaracne_workflow-831d81a2-802e-4744-9695-0ed97230935a/consensus_network_ncol_.txt
```

There is no simple standards to tell the reliability of networks. Empirically, a network with 50-300 target size is good.

# Activity-based analysis

How to understand the driver activity? Mathematically, the activity of one driver is a type of mean of the expressions of its targets. And biologically, the activity can be interpreted as a measurement of how actively the driver functions, like the enzymes in digesting their subtracts, kinase in activating their downstream genes.

## Calculate the activities

With the gene expression profiles and networks, the activity can be effortlessly calculated by `getActivity_individual()`:

```r
## let's use B cell as an example
activity_B.eset <- getActivity_individual(input_eset = pbmc14k_log2cpm.eset[, pData(pbmc14k_log2cpm.eset)$true_label.new ==
    "B"], network_file.tf = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/TF/bt100_pc001/sjaracne_workflow-fb2a69b9-f98e-47ff-87a0-6d538822fc6e/consensus_network_ncol_.txt",
    network_file.sig = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe/B/SIG/bt100_pc001/sjaracne_workflow-df798096-8dee-4baf-8f70-891c689dc769/consensus_network_ncol_.txt",
    driver_type = "TF_SIG")
```

If you need to calculate the activity for multiple groups, this is usually the case, you can do it using `getActivity_individual()` as shown above one by one and merge the esets after that. Or, scMINER privides another function, `getActivity_inBatch()`, to calculate the activity in batch:

```r
## let's use B cell as an example
activity.eset <- getActivity_inBatch(input_eset = pbmc14k_log2cpm.eset, sjaracne_dir = "/Volumes/projects/scRNASeq/yu3grp/scMINER/NG_Revision/QPan/scminer_R/Datasets/PBMC14K/SJARACNe",
    group_name = "true_label.new", driver_type = "TF_SIG", activity_method = "mean", do.z_normalization = TRUE)
#> 7 groups were found in true_label.new ...
#> Checking network files for each group ...
#> 	Group 1 / 7 : Monocyte ...
#> 		TF network check passed!
#> 		SIG network check passed!
#> 	Group 2 / 7 : B ...
#> 		TF network check passed!
#> 		SIG network check passed!
#> 	Group 3 / 7 : CD4Treg ...
#> 		TF network check passed!
#> 		SIG network check passed!
#> 	Group 4 / 7 : CD4TN ...
#> 		TF network check passed!
#> 		SIG network check passed!
#> 	Group 5 / 7 : CD4TCM ...
#> 		TF network check passed!
#> 		SIG network check passed!
#> 	Group 6 / 7 : NK ...
#> 		TF network check passed!
#> 		SIG network check passed!
#> 	Group 7 / 7 : CD8TN ...
#> 		TF network check passed!
#> 		SIG network check passed!
#> Calculating activity for each group ...
#> 	Group 1 / 7 : Monocyte ...
#> 	Activity calculation is completed successfully!
#> 	Group 2 / 7 : B ...
#> 	Activity calculation is completed successfully!
#> 	Group 3 / 7 : CD4Treg ...
#> 	Activity calculation is completed successfully!
#> 	Group 4 / 7 : CD4TN ...
#> 	Activity calculation is completed successfully!
#> 	Group 5 / 7 : CD4TCM ...
#> 	Activity calculation is completed successfully!
#> 	Group 6 / 7 : NK ...
#> 	Activity calculation is completed successfully!
#> 	Group 7 / 7 : CD8TN ...
#> 	Activity calculation is completed successfully!
#> NAs were found in the activity matrix and have been replaced by the minimum value:  -0.4401971 .
```


```r
saveRDS(activity.eset, file = "/Users/qpan/PBMC14k/DATA/activity.eset")
```


## Differential activity analysis

Similar to `getDE()`, scMINER provides a function, `getDA()`, to perform the differential activity analysis and identify the cell-type-specific drivers.


```r
## 1. To perform differential expression analysis in a 1-vs-rest manner for all groups
da_res1 <- getDA(input_eset = activity.eset, group_by = "cell_type", use_method = "t.test")
#> 7 groups were found in group_by column [ cell_type ].
#> Since no group was specified, the differential analysis will be conducted among all groups in the group_by column [ cell_type ] in the 1-vs-rest manner.
#> 	 1 / 7 : group 1 ( B ) vs the rest...
#> 	 1912 cells were found for g1.
#> 	 11693 cells were found for g0.
#> 	 2 / 7 : group 1 ( CD4TCM ) vs the rest...
#> 	 2022 cells were found for g1.
#> 	 11583 cells were found for g0.
#> 	 3 / 7 : group 1 ( CD4TN ) vs the rest...
#> 	 2505 cells were found for g1.
#> 	 11100 cells were found for g0.
#> 	 4 / 7 : group 1 ( CD4Treg ) vs the rest...
#> 	 1448 cells were found for g1.
#> 	 12157 cells were found for g0.
#> 	 5 / 7 : group 1 ( CD8TN ) vs the rest...
#> 	 2014 cells were found for g1.
#> 	 11591 cells were found for g0.
#> 	 6 / 7 : group 1 ( Monocyte ) vs the rest...
#> 	 1786 cells were found for g1.
#> 	 11819 cells were found for g0.
#> 	 7 / 7 : group 1 ( NK ) vs the rest...
#> 	 1918 cells were found for g1.
#> 	 11687 cells were found for g0.
head(de_res1)
#>      feature g1_tag      g0_tag    g1_avg   g0_avg    g1_pct    g0_pct   log2FC
#> 1251    CD3E      1 2,3,4,5,6,7  8.354660 3.874230 0.7920160 0.3819820 4.480430
#> 3820    LDHB      1 2,3,4,5,6,7  9.555670 5.614992 0.8806387 0.5458559 3.940678
#> 7765  TMEM66      1 2,3,4,5,6,7  8.604421 5.041570 0.8103792 0.5051351 3.562851
#> 1250    CD3D      1 2,3,4,5,6,7  7.281998 4.097082 0.6990020 0.3965766 3.184916
#> 1235    CD27      1 2,3,4,5,6,7  5.566280 2.428199 0.5481038 0.2482883 3.138081
#> 3992     LTB      1 2,3,4,5,6,7 10.436707 7.315803 0.9141717 0.6430631 3.120905
#>               Pval           FDR   Zscore
#> 1251 2.225074e-308  0.000000e+00 37.53784
#> 3820 4.919041e-274 6.216262e-271 35.37012
#> 7765 1.509154e-228 9.535698e-226 32.27621
#> 1250 3.193659e-174 1.130044e-171 28.13981
#> 1235 7.997857e-219 4.716603e-216 31.57555
#> 3992 5.621459e-159 1.841756e-156 26.86509
```


```r
## 2. To perform differential expression analysis in a 1-vs-rest manner for one specific group
da_res2 <- getDA(input_eset = activity.eset, group_by = "cell_type", g1 = c("B"), use_method = "t.test")

## 3. To perform differential expression analysis in a rest-vs-1 manner for one specific group
da_res3 <- getDA(input_eset = activity.eset, group_by = "cell_type", g0 = c("B"), use_method = "t.test")

## 4. To perform differential expression analysis in a 1-vs-1 manner for any two groups
da_res4 <- getDA(input_eset = activity.eset, group_by = "cell_type", g1 = c("CD4Treg"), g0 = c("CD4TCM"),
    use_method = "t.test")
```

scMINER also provides another function, `getTopFeatures()`, to easily extract the group-specific markers from the differential expression result:

```r
top_drivers <- getTopFeatures(input_table = da_res1, number = 10, group_by = "g1_tag", sort_by = "log2FC",
    sort_decreasing = TRUE)
dim(top_drivers)
#> [1] 16 11
head(top_drivers)
#>           feature g1_tag                                 g0_tag   g1_avg
#> 2472    MS4A1_SIG      B CD4TCM,CD4TN,CD4Treg,CD8TN,Monocyte,NK 2.050847
#> 1755 HLA-DPA1_SIG      B CD4TCM,CD4TN,CD4Treg,CD8TN,Monocyte,NK 2.750036
#> 619     CD79B_SIG      B CD4TCM,CD4TN,CD4Treg,CD8TN,Monocyte,NK 2.361906
#> 1390    FCER2_SIG      B CD4TCM,CD4TN,CD4Treg,CD8TN,Monocyte,NK 2.037607
#> 618     CD79A_SIG      B CD4TCM,CD4TN,CD4Treg,CD8TN,Monocyte,NK 2.091783
#> 1757 HLA-DQA1_SIG      B CD4TCM,CD4TN,CD4Treg,CD8TN,Monocyte,NK 2.054284
#>           g0_avg    g1_pct     g0_pct   log2FC          Pval FDR   Zscore
#> 2472 -0.33975624 0.9937238 0.00821004 2.390604 2.225074e-308   0 37.53784
#> 1755  0.37945455 0.9979079 0.53553408 2.370582 2.225074e-308   0 37.53784
#> 619   0.01006943 0.9973849 0.44197383 2.351836 2.225074e-308   0 37.53784
#> 1390 -0.28234219 0.9937238 0.03095869 2.319949 2.225074e-308   0 37.53784
#> 618  -0.10332717 0.9937238 0.13786026 2.195110 2.225074e-308   0 37.53784
#> 1757  0.13361395 0.9979079 0.21080989 1.920670 2.225074e-308   0 37.53784
```







