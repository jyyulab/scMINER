
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
specification.

scMINER software consists of three components:

-   [MICA](https://github.com/jyyulab/MICA) (Mutual Information based
    Clustering analysis) : perform clustering analysis

-   [SJARACNe](https://github.com/jyyulab/SJARACNe) : generate cluster
    specific networks

-   [MINIE](https://github.com/jyyulab/scMINER/)(Mutual
    Information-based Network Inference Engine) : identify
    cell-type-specific hidden drivers

## Installation

Users need to install MICA, SJARACNe and scMINER(MINIE) to run the
scMINER analysis. Using conda to create a virtual environment is
strongly recommended.

``` bash
conda create -n scminer python=3.7.6        # Create a python virtual environment
source activate scminer                     # Activate the virtual environment
## install MICA
pip install setuptools==57.5.0              # install setuptools
pip install MICA                            # Install MICA and its dependencies
## install SJARACNE
# pip install SJARACNe                        # Install SJARACNe and its dependencies
git clone https://github.com/jyyulab/SJARACNe.git
cd SJARACNe/
python setup.py build     # build SJARACNe binary
python setup.py install
cd ../
## install scMINER(MINIE)
conda install -c -r r-base=4.0.3            # Install R (Bioconductor Version 3.12)
conda install -c conda-forge r-devtools     # Install R devtools
conda install -c conda-forge r-biocmanager  # Install R BiocManager
```

In R, you can install the current version of scMINER from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("jyyulab/scMINER")
```

## Tutorial

Read the [documentation](https://jyyulab.github.io/scMINER/site/) for
detailed installation instruction and guided analysis.

### Example1: scMINER Guided Analysis on 14k PBMCs from 10x Genomics

This tutorial introduce you scMINER’s basic analysis using a PBMC
dataset with 10 sorted populations of 2k cells per population [Zheng et
al., 2017](https://www.nature.com/articles/ncomms14049). Check [tutorial
for example](https://jyyulab.github.io/scMINER/site/tutorials/PBMC-14k/)

### Example2: scMINER Guided Analysis on WT and KO CD8+ T cell in chronic infection model

TOX is a master transcription factor for CD8+ T cell exhaustion during
chronic infection. This tutorial introduce you scMINER’s basic analysis
using a WT and TOX KO CD8+ T dataset (GSE119940) [Yao et al., Nat
Immunol 2019](https://www.nature.com/articles/s41590-019-0403-4). Check
[tutorial for
example](https://jyyulab.github.io/scMINER/site/tutorials/CD8T/)

If you’d like to contribute, please open an issue or a pull request in
the [github repository](https://github.com/jyyulab/scMINER/issues).
