
---
layout: default
title: Online Example
nav_order: 2
---


# Example(PBMC)
_______


## Run MICA
* Run MICA pipeline (MIE and MICA) using filtered PBMC data after working environment was set up [In LSF]
* Sample input(Cell by Gene Square matrix): "PBMC_Demo_MICA_input.txt"


**First run MIE :**   
`[PYTHON_PATH] [scMINER.py] MIE Pipeline PBMC_demo [input.txt] [output_dir] PBMC_demo`

**Then MICA:** `[PYTHON_PATH] [scMINER.py] MICA Clust PBMC_demo [path_to_PBMC_demo.whole.h5] [path_to_PBMC_demo_mi.h5] [output_dir] PBMC_demo --k 3 4 5 --perplexity 30 --retransformation FALSE`



## Define clusters and Prepare for MINIE(SJARACNE)

*  First of all, install R package **'MINIE'**
*  MICA outputs a txt file with t-SNE coordination and cluster membership of all cells, User can load the input and output together in a expressionset in **R** with Function: 

```
eset.demo <- readMICAoutput(input_file=[Your_input_txt], 
								 load_clust_label = TRUE, #whether to load your clustering result
			   					 output_file=[Your_scMINER_output_txt])		 
```

* To visualize, you can use either

```   					
gene_highlighting(input_eset = eset.demo, target = [target_gn_selected], title.size = 8)
#or
gene_vlnplot(eset.demo, target=gn.sel, group_tag = "label")
#only if "load_cluster_label=TRUE in previous function or manually add your label to eset"
```

* Celltype annotation (manually or assigned by pre-selected marker genes)








* To generate SJARACNE input, you need a list of knowledge based TF genes

```
generateSJAracneInput( eset=eset.demo, # input expressionset
							tf.ref=tf.ref, #knowledge base TF list in geneSymbol
							wd.src="Sjaracne/", # output path
							group_tag="celltype" ) # For cluster/subtype based network
							
```



## Run SJARACNE
*  Run SJARACNE   

```

```




## Get network inferred Activity and find master regulators
*  TF Activities are calculated by allocating expressions from their downstream targets

```
acs.demo <- GetActivityFromSJARANCE(
						SJaracne_output_path="Sjaracne/", # path to sjaracne output folder
						SJaracne_input_eset = eset.demo, # your expressionset
						activity.method = "unweighted", # default for scRNAseq data
						activity.norm = TRUE, #default for scRNAseq data
						save_network_file = FALSE, 
						save_path = NA)

```						  
*  Find master regulators by identifying highly variable TF

```
res.demo <- aov.test( eset = acs.demo, # activity set 
						  group_tag = "celltype", # group case
						  trouble_shooting = FALSE) # default 

```






###R session info

```
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












