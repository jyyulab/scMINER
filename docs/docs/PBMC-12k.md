---
layout: default
title: PBMC(12k) scRNA-seq data
nav_order: 4
---

# Analysis on PBMC(12k) scRNA-seq data via scMINER
{:.no_toc}
Here we demonstrate our pipeline using PBMC (10x genmomics) scRNA-seq data [link to data matrix]. Full data contains 68k cells(link to 10x website), in order to provide a quicker guidance, we've down sampled this data to 12k cells.
Original data website can be downloaded here: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/fresh_68k_pbmc_donor_a


## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---
## Data preprocessing
This can be done by any scRNA-seq preprocessing pipeline. We encourage user to feed in all genes from your data for MICA, instead of highly variable genes only. Here, in order to stick to the focus, we only demonstrate a quick function to conduct gene/cell filtering, without any data exploratory visualization.

```R
eset.demo<-pre.MICA(data.input=,
				)
```




## Run MICA clustering

MICA is implemented in Python now. 

This will give you one visualization for each choice of k.



![]()


## Install MINIE package in R
{: .d-inline-block :}

**MINIE** analysis was wrapped up as a R package to help bridge unsupervised clustering and gene regulatory network analysis.  

> **_Note:_** Detailed information about individual functions are documented in package manual.


## Cell type analysis from MICA output
{: .d-inline-block :}

First, after clustering via MICA(see [MICA] ({{site.baseurl }}{% link docs/MICA.md %}), you can load MICA output (in .txt) as well as input expression matrix in R under an expressionSet. This is going to be the major data structure we used for MINIE analysis in R.

### Reading MICA output
{: no_toc }

```R
eset.demo <- readMICAoutput( 
	input_file="PBMC12k_MICA_input.txt",
	load_clust_label=TRUE, 
	output_file="PBMC12k.ggplot.txt")
```


### Marker gene highlighting
{: no_toc }

Picked marker genes could be visualized on t-SNE scatterplot, heatmap or violinplot. This will help pick up a reasonable number of cluster.

```R
gn.sel <- c("GZMK","GZMH","GZMA","CCR7","CD8A","SELL")
gene_highlighting(input_eset=eset.demo, target = gn.sel, title.size = 8)
```

![]()


```R
gene_vlnplot(eset.demo,target=gn.sel,group_tag = "label")
```


![]()


```R
gene_heatmap(eset = eset.demo,target = gn.sel,group_tag = "label",
			 save_plot = TRUE,width = 6,height = 6,
             name = "log2_expression",plot_name="./GeneHeatmap.png")
```

![]()


### Assign cell type to cluster
{: no_toc }

Here we curated a reference signature list of 8 immune cell types(link) for cell type annotation. In `AssignCellTypes.hmp` function, we calculated cell type scores for each clusters, and visualize scores using heatmap. 

```R
ref<-read.xlsx("Immune_signatures.xlsx")
hmp<-AssignCellTypes.Hmp(ref=ref,eset=eset.demo,save_plot = TRUE)

# Manually assign your cell type label
celltype<-c("MemoryT","NaiveT","CD8em","CD8eff")
eset.demo$celltype <- celltype[eset.demo$label]
```

## Network generation via SJARACNe
{: no_toc }

In order to generate cell type/group/cluster specific network, group information should be stored under `pData([your_expressionSet])`. And R function `generateSJAracneInput` will help to partition your expression matrix and conduct a loose filtering of your scRNA-seq data(filter about 0 expressed genes in cluster). Besides, a reference TF list should be provided as `tf.ref` to guide hub gene selection. Each group will create one directory which contains filtered expression matrix in .exp format, as long as the filtered TF list in .txt. 

```R
generateSJAracneInput(eset=eset.demo,tf.ref=tf.ref,wd.src="Sjaracne/", group_tag="celltype")
```
_Warning:_ SJARACNe has not been integrated into MINIE yet, please consult [here](https://github.com/jyyulab/SJARACNe) to run SJARACNe for network generation on scRNA-seq data.


## Find cell type specific master regulator 
{: no_toc }

Identify master regulator from content based network is the key step in MINIE to help understanding your scRNA-seq data.  


### Calculate Inferred activity
{: no_toc }
TF acitivities are calculated by integrating expression profile of their targets. Targets identified from SJARACNe of perticular TF was normalized and averaged to infer TF activity.

```R
acs.demo <- GetActivityFromSJARANCE(
			 SJaracne_output_path="Sjaracne/",
			 SJaracne_input_eset=eset.demo,
			 activity.method="unweighted", 
			 activity.norm=TRUE, 
			 group_tag = "celltype",
			 save_network_file=FALSE, save_path=NA)
```


### Find Differential activity TF
{: no_toc }

The function `FindDAG` was designed for identify highly differentiated TF from SJARACNe inferred activity matrix. In order to do so, we did two sided student's t-test to compare mean acitivty from one cell type V.S. the others. 

```R
res <- FindDAG(eset = acs.demo,group_tag = "celltype")

```
This function will output a full matrix that contians all TF occurred in original dataset, statistics such as t.statistics, p-value, 95%CI, etc. are outputed to help idenify master regulators. 

```R
print(head(res))




```

You can also visualize top master regulator candidates in heatmap or violinplots. Only png as plotting device is supported.

```R
gn.sel <- TopMasterRegulator(res)

gene_heatmap(eset = eset.demo,target = gn.sel,group_tag = "label",
			 width = 6,height = 6, save_plot=TRUE,
             name = "log2_expression",plot_name="./TopTFHeatmap.png")
```

![](https://guides.github.com/activities/hello-world/branching.png)

```

---
## R session info


R version 3.5.0 (2018-04-23)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] bindrcpp_0.2.2      Biobase_2.42.0      BiocGenerics_0.28.0 dplyr_0.7.8        
[5] RColorBrewer_1.1-2  ggplot2_3.1.0       reshape2_1.4.3     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.0       rstudioapi_0.8   bindr_0.1.1      magrittr_1.5     tidyselect_0.2.5
 [6] munsell_0.5.0    colorspace_1.3-2 R6_2.3.0         rlang_0.3.0.1    stringr_1.3.1   
[11] plyr_1.8.4       tools_3.5.0      grid_3.5.0       gtable_0.2.0     withr_2.1.2     
[16] digest_0.6.18    yaml_2.2.0       lazyeval_0.2.1   assertthat_0.2.0 tibble_1.4.2    
[21] crayon_1.3.4     zip_1.0.0        purrr_0.2.5      glue_1.3.0       labeling_0.3    
[26] openxlsx_4.1.0   stringi_1.2.4    compiler_3.5.0   pillar_1.3.1     scales_1.0.0    
[31] pkgconfig_2.0.2 
```
---
