# scMINER (Single-cell Mutual Information based Network Engineering Ranger)

=======

*scMINER* is a **systems biology** analysis framework for high-throughput single cell RNA-seq data implemented primarily in R and python. This package offers a combination of several individual tools including but not limit to [MICA](https://github.com/jyyulab/MICA) (Mutual Information based Clustering analysis) and [SJARACNe] (https://github.com/jyyulab/SJARACNe). Installation instructions for each individual tools are available on github through links above.


## Installation
To install scMINER R package, in R environment with `R >= 3.6.0` and above, run

```R
devtools::install_github("jyyulab/scMINER",ref='master',dependencies='Depends') 
```
or clone from `https://github.com/jyyulab/scMINER/` and install from source files like below:

```R

devtools::install(pkg='.',dependencies=TRUE) ## Install the package with dependencies.
devtools::install_deps(pkg = ".", dependencies = TRUE) ## Install package dependencies if needed.

```
If you are a Yu Lab member, tar.gz file is under project space, you can install via:

```R
devtools::install(pkg='/Volumes/yu3grp/scRNASeq/yu3grp/scMINER/scMINER_0.1.0.tar.gz', dependencies=TRUE)

```


## Documentation and guided anlaysis
Comprehensive documentations could be found at:  
- https://jyyulab.github.io/scMINER/


