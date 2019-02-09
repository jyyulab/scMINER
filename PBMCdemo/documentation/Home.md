#HOME

##Installation

###Install Python 3

###Install R 3.5.0
If you don't have R preinstalled, please infer below official website of R:  
https://www.r-project.org/

###Install R pacakges 
In order to smoothly run scMINER, in R console, please install following packages:

```
install.packages(c("dplyr","plyr","ggplot2","Rcolorbrewer","reshape2","BiocGenerics"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase", version = "3.8") # For R version >= 3.5

```


###Online graphic interface



























