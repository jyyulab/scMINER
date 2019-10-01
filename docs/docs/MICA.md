---
layout: default
title: Clustering by MICA
nav_order: 3
---

# Clustering Analysis with MICA
{:.no_toc}
MICA(Mutual Information based Clustering Analysis) is a nonlinear clustering analysis tool designed for scRNA-seq data. To install MIE (Mutual inforamtion estimator for distance matrix generation, module required for MICA) and MICA, please refer our [MIE](https://github.com/jyyulab/MIE) and [MICA](https://github.com/jyyulab/MICA) github page.  
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }
 
1. TOC
{:toc}

---
## Preprocssing
Preprocessing is fairly simple for MICA analysis. You can use our customized script in R, with function `draw.scRNAseq.QC` and `preMICA.filtering` in scMINER R package. For detailed information, please see complementary tutorial in tab `Sample Analysis with PBMC(12k) scRNA-seq data`.

## Basic usage
MICA is implemented in python. For those who are not familiar with python a well-designed function in scMINER R package could help generate command for   you could use function `generate_MICA_rmd` in R package `scMINER` to generate essential command for running MICA locally:

In R console: 
```R
scMINER::generate_MICA_cmd(save_sh_at, #path to save shell script 
                  			input_file, #your MICA input file
                  			project_name, 
                  			num_cluster, #a vector of numerical number
                  			output_path, #path to MICA output
                  			host="local", 
		                    visualization="tsne" #or "umap")
```


or, you can create your own shell script to run MICA like below: 

```SHELL
#!/usr/bin/env bash
mica local \
-i ./test_data/inputs/PBMC_Demo_MICA_input_mini.txt \
-p "test_local" \
-k 3 4 \
-o ./test_data/outputs/test_local/ \
```

## MICA Outputs

Each assigned number of k will output one folder containing following files.

1. `[Project_name]_k[number]_tsne.png`  --visualization of clustering result (default as tSNE)

  <img src="./plots/2_0_cwl_local_k3_tsne.png" width="600"/> 

2. `[Project_name]_dist.h5`  -- h5 file containing distance matrix calculated.
3. `[Project_name]_mds.pdf`  -- pdf file of t-SNE visualization of mds transformed distance matrix, with perplexity set to 30
4. `[Project_name]_tsne_ClusterMem.txt`  -- txt file containing visualization coordinates and clustering labels


## Useful parameters

### Visualize with U-map or t-SNE
MICA incorporate [UMAP](https://umap-learn.readthedocs.io/en/latest/parameters.html) as optional clustering visualization, with `min_dist` parameter set to `0.25`, this controls how points packed together. Low values of min_dist will result in clumpier embeddings. You can tune this parameter with :

```SHELL
--min_dist 0.1 (or other number ranging from 0-1) 
```

tSNE visualization is our default visualization method in the pipeline, if you want to use t-SNE, just set :

```SHELL
--visualization tsne (all lower cap, no "-")
```
and you can also set parameter (perplexity) for tsne using

```SHELL
-pp 20 (or any other integers larger than 5)
```

### Try other dimension reduction methods
MICA also incorporated other dimension reduction methods such as pca or lpl, 
you can use them via adding parameter:

```SHELL
-dr PCA  (or: MDS | PCA | LPL | LPCA) 
```

### Try other distance matrix calculation methods
MICA also incorporated other dimension reduction methods such as pca or lpl, 
you can use them via adding parameter:

```SHELL
--dist MI  (or: euclidean | spearman | pearson)
```

## Post-clustering analysis
We offer a handful of useful functions in scMINER ranging from visualization to driver estimation to help you explore your scRNA-seq data in a system biology way after clustering. 