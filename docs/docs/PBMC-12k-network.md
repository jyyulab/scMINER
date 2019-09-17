---
layout: default
title: Cell-type specific network analysis
nav_order: 6
---

# Analysis on cell-type specific network
{:.no_toc}
Here we demonstrate our advanced downstream analysis pipeline using PBMC (10x genmomics) scRNA-seq data [link to data matrix]. Full data contains 68k cells(link to 10x website), in order to provide a quicker guidance, we've down sampled this data to 12k cells. This tutorial provide a guidance for study master regulators utilzing cell-type specific networks. 

Original data website can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/fresh_68k_pbmc_donor_a).



##Network visualization
###


































---
```
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
