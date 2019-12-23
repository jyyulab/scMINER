---
layout: default
title: Get started
nav_order: 2

---
# Get started
{:.no_toc}
--- 

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}
---  

## Overview
**scMINER** is a **systems biology** toolkit for QC and analysis of high-throughput single cell RNA-seq data. 
For the better understanding of the sources of heterogeneity from single-cell RNA-seq data, scMINER enable users 
to identify and interpret **cell type-specific hidden drivers**, where the drivers can be either transcription 
factors (TF) or signaling factors (SIG). 

**scMINER** also includes new mutual information based approach [MICA](https://github.com/jyyulab/MICA) for classifying 
cells and improved scalable version of [SJARACNe](https://github.com/jyyulab/SJARACNe) for cell type-specific 
gene network reverse engineering.
{: .mb-lg-10 }
<img src="../docs/plots/scMINER_Overview.png" alt="Overview"/> 


## Installation

### Dependencies
The entry point of scMINER is a R package, which includes some of the major functions such as draw.scRNAseq.QC, generateMICAcmd, etc. It has been
throughly tested in R 3.6.1. So it is highly recommended to have **R version 3.6.1 or higher** in your environment.  
scMINER also includes two independent python packages [MICA](https://github.com/jyyulab/MICA) and 
[SJARACNe](https://github.com/jyyulab/SJARACNe), for which **python 3.6.1 or higher** is required.

### Installing scMINER

* Install [scMINER](https://github.com/jyyulab/scMINER) R package from github

```R
#install dev_tool first install.packages(devtools)
devtools::install_github("jyyulab/scMINER") 
#or
devtools::install_local('scMINER_0.0.1.tar.gz')
```

* Install [MICA](https://github.com/jyyulab/MICA)

 Installing the official package from PyPi:

```bash
pip install MICA
```

 Or you can install from source: 

```bash
git clone https://github.com/jyyulab/MICA
cd MICA
python setup.py install
```

* Install [SJARACNe](https://github.com/jyyulab/SJARACNe)

 Installing the official package from PyPi:
```bash
pip install SJARACNe
```

 Or you can install from source: 
```bash
git clone https://github.com/jyyulab/SJARACNe
cd SJARACNe
python setup.py build
python setup.py install
```

## Interoperability
### SparseExpressionSet

Majority of data processing/visualizaiton in scMINER are completed under R environment, utilizing SparseExpressionSet object. It is a customized data object adapted from `ExpressionSet` class. The only thing that is different between these two classes is users are able to store sparsematrix (class dgTmatrix, etc.) in the object, for scRNA-seq data, it will boost up computational performances.


To import scRNA-seq data to R, we provide an easy function `readscRNAseqData `:  

For standard 10x genomics output, you can set `is.10x=TRUE`, it will take matrix.mtx, feature.tsv, barcodes.tsv as input.  
You may want to set `CreateSparseEset = F` at the very beginning in order to manually check if your data was correctly imported. Imported data will be stored in a list.

```R
d.12k <- readscRNAseqData(file="../PBMC12k_input/", is.10x = T, CreateSparseEset = F, add.meta = F)
```

If your data is not directly from 10x genomics standard output, but random public dataset or even in other object format -- just make sure you have your expression matrix as a matrix or sparsematrix, feature data and phenotype data as data frame, then use `CreaseSparseEset` function to create scMINER object. By doing this step, **no** filtering, transformation, or preprocessing is performed.

```R
eset.12k<-CreateSparseEset(data=d.12k$raw.data,feature.data = d.12k$feature.data, add.meta = T)
```

