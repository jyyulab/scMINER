---
layout: default
title: Welcome to scMINER
nav_order: 1
permalink: /
---


# scMINER tutorial
{: .fs-9 }

scMINER is a toolbox for single-cell analysis based on mutual information. This package is a combination of several individual tools including but not limit to [MICA](https://github.com/jyyulab/MICA) and [MINIE](). Separate installation of each individual tools is also available on github. 

{: .fs-6 .fw-300 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/jyyulab/scMINER){: .btn .fs-5 }

---

## Getting started
### Dependencies
scMINER is depend on [python3](https://www.python.org/downloads/) and [R](https://www.r-project.org/). 

### Local installation: 
1. Install scMINER:

2. Install scMINER-MINIE:

```R
#install dev_tool first install.packages(devtools)
devtools::install_github("jyyulab/MINIE") 

```

3. Install R dependencies:

```R
install.packages(c("dplyr","plyr","ggplot2","Rcolorbrewer","reshape2","BiocGenerics"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase", version = "3.8") # For R version >= 3.5

```


---

## About the project


### License

Just the Docs is distributed by an [MIT license](https://github.com/jyyulab/scMINER/blob/master/LICENSE).

### Contributing

When contributing to this repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change. Read more about becoming a contributor in our GitHub repo.