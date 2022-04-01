# Installation
## Dependencies
scMINER requires R version >= 4.0.5 and Python version >= 3.8.10 to run. We recommend using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) to create and manage the environments.


## Installing scMINER
###### Install [scMINER](https://github.com/jyyulab/scMINER) R package

To install from GitHub `master` branch

```R
install.packages(devtools)  # install devtools package first
devtools::install_github("jyyulab/scMINER", ref='master', dependencies='Depends') 
```
or to clone from `https://github.com/jyyulab/scMINER/` and install from source files:

```R
devtools::install(pkg='.',dependencies=TRUE)            # Install the package with dependencies.
devtools::install_deps(pkg = ".", dependencies = TRUE)  # Install package dependencies if needed.
```

or to download the released tar.gz file is under project space, you can install via:

```R
devtools::install(pkg='scMINER_0.1.0.tar.gz', dependencies=TRUE)
```

###### Install [MICA](https://github.com/jyyulab/MICA) python package

To install the package from PyPI:

```bash
pip install MICA
```

or to install from source: 

```bash
git clone https://github.com/jyyulab/MICA
cd MICA
python setup.py install
```

###### Install [SJARACNe](https://github.com/jyyulab/SJARACNe) python package

As the core SJARACNe functions were implemented in C++ which requires compiling, we recommend using installing from source: 
```bash
git clone https://github.com/jyyulab/SJARACNe
cd SJARACNe
python setup.py build
python setup.py install
```

###### R session info
```R
> sessionInfo()
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] stats4    grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scMINER_0.1.0               NetBID2_2.0.1               openxlsx_4.1.0.1            msigdbr_7.0.1               GSVA_1.32.0                 biomaRt_2.40.4             
 [7] reshape_0.8.8               arm_1.10-1                  MASS_7.3-51.4               MCMCglmm_2.29               ape_5.3                     coda_0.19-3                
[13] ordinal_2019.4-25           umap_0.2.3.1                plotrix_3.7-6               plot3D_1.1.1                igraph_1.2.4.1              aricode_0.1.2              
[19] ConsensusClusterPlus_1.48.0 DESeq2_1.24.0               tximport_1.12.3             impute_1.58.0               limma_3.40.6                GEOquery_2.52.0            
[25] lme4_1.1-21                 SummarizedExperiment_1.14.1 DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.55.0          GenomicRanges_1.36.1       
[31] GenomeInfoDb_1.20.0         IRanges_2.18.3              S4Vectors_0.22.1            kableExtra_1.1.0            knitr_1.25                  Matrix_1.2-17              
[37] rmarkdown_1.16              cowplot_1.0.0               dplyr_0.8.3                 pheatmap_1.0.12             ComplexHeatmap_2.0.0        scales_1.0.0               
[43] RColorBrewer_1.1-2          ggplot2_3.2.1               reshape2_1.4.3              Biobase_2.44.0              BiocGenerics_0.30.0         BiocManager_1.30.4         

loaded via a namespace (and not attached):
  [1] backports_1.1.4        circlize_0.4.8         Hmisc_4.2-0            plyr_1.8.4             lazyeval_0.2.2         GSEABase_1.46.0        splines_3.6.1          digest_0.6.21         
  [9] htmltools_0.3.6        magrittr_1.5           checkmate_1.9.4        memoise_1.1.0          cluster_2.1.0          readr_1.3.1            annotate_1.62.0        askpass_1.1           
 [17] prettyunits_1.0.2      colorspace_1.4-1       blob_1.2.0             rvest_0.3.4            xfun_0.10              jsonlite_1.6           crayon_1.3.4           RCurl_1.95-4.12       
 [25] graph_1.62.0           genefilter_1.66.0      zeallot_0.1.0          survival_2.44-1.1      glue_1.3.1             gtable_0.3.0           zlibbioc_1.30.0        XVector_0.24.0        
 [33] webshot_0.5.1          GetoptLong_0.1.7       Rhdf5lib_1.6.1         shape_1.4.4            abind_1.4-5            DBI_1.0.0              Rcpp_1.0.2             viridisLite_0.3.0     
 [41] xtable_1.8-4           progress_1.2.2         htmlTable_1.13.2       clue_0.3-57            reticulate_1.13        foreign_0.8-72         bit_1.1-14             Formula_1.2-3         
 [49] htmlwidgets_1.3        httr_1.4.1             acepack_1.4.1          pkgconfig_2.0.3        XML_3.98-1.20          nnet_7.3-12            locfit_1.5-9.1         tidyselect_0.2.5      
 [57] rlang_0.4.0            later_0.8.0            AnnotationDbi_1.46.1   munsell_0.5.0          tools_3.6.1            RSQLite_2.1.2          evaluate_0.14          stringr_1.4.0         
 [65] bit64_0.9-7            zip_2.0.4              purrr_0.3.2            nlme_3.1-141           mime_0.7               xml2_1.2.2             compiler_3.6.1         shinythemes_1.1.2     
 [73] rstudioapi_0.10        png_0.1-7              tibble_2.1.3           geneplotter_1.62.0     stringi_1.4.3          cubature_2.0.3         lattice_0.20-38        nloptr_1.2.1          
 [81] tensorA_0.36.1         vctrs_0.2.0            pillar_1.4.2           lifecycle_0.1.0        GlobalOptions_0.1.1    ucminf_1.1-4           data.table_1.12.2      bitops_1.0-6          
 [89] corpcor_1.6.9          httpuv_1.5.2           R6_2.4.0               latticeExtra_0.6-28    promises_1.0.1         gridExtra_2.3          boot_1.3-23            assertthat_0.2.1      
 [97] rhdf5_2.28.0           openssl_1.4.1          rjson_0.2.20           withr_2.1.2            GenomeInfoDbData_1.2.1 hms_0.5.1              rpart_4.1-15           tidyr_1.0.0           
[105] minqa_1.2.4            misc3d_0.8-4           numDeriv_2016.8-1.1    shiny_1.3.2            base64enc_0.1-3 
```

###### Python session info
```Python
>>> import MICA
>>> import session_info
>>> session_info.show()
-----
MICA                NA
session_info        1.0.0
-----
Python 3.8.10 | packaged by conda-forge | (default, Sep 13 2021, 21:46:58) [GCC 9.4.0]
Linux-3.10.0-1160.15.2.el7.x86_64-x86_64-with-glibc2.10
-----
Session information updated at 2022-02-15 11:56
```
