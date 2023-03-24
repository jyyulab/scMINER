
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMINER

<!-- badges: start -->
<!-- badges: end -->

*scMINER* is an R package for preprocessing, QC, clustering, and *hidden
driver analysis* of single-cell RNA-seq data. scMINER enables mutual
information-based cell clustering, cell-type-specific gene regulatory
network (GRN) reverse engineering and protein activity inference, to
identify hidden transcriptional factors (TFs) and signaling factors
(SIGs) driving cellular lineage differentiation and tissue specific
specification. scMINER software consists of three components:
[MICA](https://github.com/jyyulab/MICA) (Mutual Information based
Clustering analysis), [SJARACNe](https://github.com/jyyulab/SJARACNe)
and [MINIE](https://github.com/jyyulab/scMINER/).

## Installation

You can install the development version of scMINER from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jyyulab/scMINER")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scMINER)
#> 载入需要的程辑包：Biobase
#> 载入需要的程辑包：BiocGenerics
#> 
#> 载入程辑包：'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which.max, which.min
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 载入需要的程辑包：ggplot2
#> Warning: 程辑包'ggplot2'是用R版本4.1.3 来建造的
## basic example code
```

## Tutorial

Read the [documentation](https://jyyulab.github.io/scMINER/site/) for
detailed installation instruction and guided analysis. If you’d like to
contribute, please open an issue or a pull request in the [github
repository](https://github.com/jyyulab/scMINER/issues).
