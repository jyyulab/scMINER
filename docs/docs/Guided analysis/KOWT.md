---
layout: default
title: KO vs WT analysis on Th17
parent: Guided tutorials
nav_order: 6
---

# KO vs WT analysis in Th17 scRNA-seq data
{:.no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---
## Demo data
{: .d-inline-block :}

Here we demonstrate a KO vs WT type of analysis using public dataset from a nature paper published in late 2018. Unprocessed data could be downloaded under accession number GSE121599. 

## Step 1: Data preprocessing

In this project, we have 6 standard scRNA-seq outputs, you can merge them first then perform quality control and data filtering. 

### Read 10x genomics data
```R
d<-readscRNAseqData(file = "01_Quality_Control/",add.meta = T,is.10x = TRUE,CreateSparseEset = T)

#Manually checking phenotype and feature data to see if data was read correctly
head(pData(d))
head(fData(d))

d$group<-sapply(strsplit(colnames(d),"_",fixed = T),"[",1) # store sample info in group 
```

### Perform Quality control

`draw.scRNAseq.QC` will generate a html file with essential visualization to assist your quality control and data filtering.

```R
cfs <- draw.scRNAseq.QC(SparseEset = d,project.name = "Th17",plot.dir = "./01_Quality_Control/")

# make adjustment to computed threshold
# no additional filtering added
cfs$ERCC_cf<-Inf
cfs$mito_cf<-0.1
cfs$umi_cf_hi<-40000

eset.sel<-preMICA.filtering(SparseEset = d,cutoffs = cfs)

draw.scRNAseq.QC(SparseEset = eset.sel,project.name = "Th17_afterQC",plot.dir = "./01_Quality_Control/",output.cutoff = FALSE)
```

### Normalization and log2 transformation
After re-run quality control function, you can go ahead and normalize your data. Here we recommend to perform CPM normalization and log2 transformation.

> **Note:** Please clean up your working environment if needed, to reduce memory usage.

```
norm = 1e6 
exp.norm <- sweep(exprs(eset.sel), 2, norm/unname(Matrix::colSums(exprs(eset.sel))), '*')
exp.log2 <- log(exp.norm+1,base=2)

eset.log2 <- CreateSparseEset(data=exp.log2, 
                              meta.data = pData(eset.sel), 
                              feature.data = fData(eset.sel), 
                              add.meta = F)

rownames(eset.log2)<-NULL

#Don't forget to save your working image
save.image(file = "02_MICA/Th17_beforeMICA.RData")
```

## Step 2: Run Clustering via MICA

After data preprocessing, you can use `generateMICAinput` to generate standard MICA input, this can be outputed in two types of format tab-delimited text file or h5 file. 

```R
#generate MICA input as h5 file
generateMICAinput(d= exprs(eset.log2) ,filename="02_MICA/Th17_MICA_input.h5")
```

In order to generate runnable command to perform MICA, we could use `generateMICAcmd` function. With data at this size, we recommend you to run MICA on High Performance Computing Facility(HPCF).

```R
generateMICAcmd(save_sh_at = "02_MICA/",
                input_file = "02_MICA/Th17_MICA_input.txt",
                config_file = "02_MICA/config_cwlexec.json",
                project_name = "Th17",num_cluster = c(6,7,8),
                host = "lsf",output_path = "./",queue = "standard")
```

## Step 3: Composition analysis 

###Read MICA output
After loading your pre-saved preprocessing result, you can load MICA output into your eset via function `readMICAoutput`. 

```R
load("02_MICA/Th17_beforeMICA.RData")

eset.log2<-readMICAoutput(eset=eset.log2,output_file = "Th17_k8_tsne_ClusterMem.txt",load_ClusterRes = T)
head(pData(eset.log2))

eset.log2$condition<-gsub("[0-9]","",eset.log2$group)
eset.log2$replicates<-gsub("[A-Z]|[a-z]","",eset.log2$group)

eset.log2$tSNE_1<-eset.log2$X 
eset.log2$tSNE_2<-eset.log2$Y

```

## Step 4: Feature visualization



## Step 5: KO/WT specific network generation
### Network validation
### Infer Activity


## Step 6: Context specific hidden driver identification


## Step 7: Functional enrichment study  





