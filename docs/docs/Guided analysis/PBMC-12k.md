---
layout: default
title: Guided analysis on PBMC(12k)
parent: Guided tutorials
nav_order: 5
---

# Analysis on PBMC(12k) scRNA-seq data via scMINER
{:.no_toc}

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---
## Demo data
{: .d-inline-block :}

Here we demonstrate our pipeline using PBMC (10x genmomics) scRNA-seq data. Full data contains 68k cells, in order to provide a faster guidance, we've down sampled this data to 12k cells. Original data can be downloaded [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/fresh_68k_pbmc_donor_a).

## Step 1: Data preprocessing
{: .d-inline-block :}

### Read 10x genomics data
Read 10x genomics data with function `readscRNAseqData `in scMINER package. This function could help read data from either 10x genomics standard output(usually contains three individual files: matrix.mtx, barcodes.tsv, features.tsv), or other text files by passing arguments to `read.delim()`. 

This function helps create a Sparse matrix object using Expressionset prototype, If set `CreateSparseEset=T`. Otherwise, it will create a list object that stores expression data, feature data and sample data under three separate slots. If `add.meta=T`, then additional sample info such as total number of UMI will be calcualated and outputed in sample data. Here, we defined `is.10x=T`,`CreateSparseEset = F`, and `add.meta=F`.

```R
#do not create SpearseEset or add meta this time
#please confirm that data was sucessfully read
d.12k <- readscRNAseqData(file="../PBMC12k_input/",is.10x = T, CreateSparseEset = F, add.meta=F)
```

After loading data to environment properly, you can now create Sparse Matrix expression by using `CreateSparseEset`function:

```R
# please set add.meta=T if you want to run our quality control pipeline
eset.12k<-CreateSparseEset(data=d.12k$raw.data,feature.data = d.12k$feature.data, add.meta = T)

# Define your group information:
# If there is no additional sample info available now, you can just put a character here
eset.12k$group<-sapply(strsplit(eset.12k$cellName,"-"),"[",2) #user define, optional

```

### Quality control and data filtering
Quality control assessments could be done using `draw.scRNAseq.QC` function, which outputs a html report generated through Rmarkdown([PBMC12K_QC_report.html](./PBMC12k_scRNAseq_QC.html)). The report includes three essential quality control figures at both gene and cell level. Suggested cutoff will be returned as a list if set `output.cutoff=TRUE`.

```R
cutoffs <- draw.scRNAseq.QC(SparseEset=eset.12k, 
                          project.name = "PBMC12k",
                          plot.dir = "./QC/",
                          group = "group", # this indicate which meta data information will be use in x axis to group violin plots
                          output.cutoff = TRUE) #whether or not to output suggested cutoffs
```
The first plot is a histogram which helps visualize distribution of expressed genes among each cells. Blue veritcal line shows the recommended cutoff, which is calculated by 0.5% * number of cells. Genes expressed lower number of cells than threshold should be filtered.

<center><img src="./plots/1_1_Gene_QCmetrics_before_filtering.png" alt="drawing" width="700"></center>

The second plot helps visualize total UMI count, and total number of gene expressed in violin plots. Horizontal blue line indicates suggested high/low cutoffs. Suggested thresholds were computed based on Median ± 3 * MAD (Maximum absolute deviance). Suggested threshold numbers are also printed right above blue lines in labels.
![](../plots/1_2_Cell_QC_1.png)

The third plot visualizes mitochondrial gene expression percentage, and spike-in genes expression percentage for each cell V.S. total number of UMI in scatter plots. Cells with high percentage of mitochondrial gene expression but low total number of UMI count are often considered as low quality.
![](../plots/1_3_Cell_QC_2.png)


Then you could perform filtering with function `preMICA.filtering`. We recommend to input `cutoffs` using thresholds list which directly outputted from `draw.scRNAseq.QC` functon. You could also manually change cutoffs by re-assign thresholds in `cutoffs` list, e.g. `cutoffs$umi_cf_hi <- Inf` means not doing filtering on outliers caused by high total UMI value.

```R
#manually adjust your thresholds
cutoffs$umi_cf_hi <- Inf  # only filter on low total number of UMI

# do filtering
eset.sel <- preMICA.filtering(SparseEset = eset.12k, cutoffs = cutoffs) 
```
 
### Normalization and transformation
In scMINER package, we don't offer functions to perform data normalizaton. You can use your own prefered normalization method. However, **we usually  do CPM and log2 transformation (required) for MICA input**.

```R
norm = 1e6 
exp.norm <- sweep(exprs(eset.sel), 2, norm/unname(Matrix::colSums(exprs(eset.sel))), '*')

# log transformation
# Required for MICA
exp.log2 <- log(exp.norm+1,base=2)

# save as SparseEset
eset.log2 <- CreateSparseEset(data=exp.log2, meta.data = pData(eset.sel), feature.data = fData(eset.sel), add.meta = F)

```


## Step 2: Perform clustering analysis via MICA
{: .d-inline-block :}

MICA was implemented in Python. If you would like to install MICA, please refer to [MICA github page](https://github.com/jyyulab/MICA). There are several handlers for you to choose in MICA for better visualization. A more comprehensive tutorial could be found under [Clustering with MICA](./MICA.md) tab . **Here we suggest saving your working directory prior to running MICA**. 

### Generate MICA input and command
After reviewing all visualizations and finished filtering, you can go ahead and generate clustering (MICA) input with function `generateMICAinput`. This function takes an expression matrix as input, and outputs a cell by gene .txt or .h5 file. Please note that **you should always feed MICA the log or log2 transformed data**.

```R
generateMICAinput(d = exp.log2 ,filename="PBMC12k_MICA_input.txt")

# clean your working environment
rm(exp.log2);rm(exp.norm);
```


We also offer a function called `generate_MICA_cmd ` to help write MICA command in a shell script. In order to run MICA on LSF, you need to set `host=lsf`,  define `queue = [your queue]` (required), and `config_file=[path to your customized config file]`(optional, default as NULL). In `num_cluster`, you can specify a vector of number of K to perform clustering analysis for different number of cluster simultaneously. 

```R
generateMICAcmd(save_sh_at = "./",
                input_file = "./PBMC12k_MICA_input.txt",
                project_name = "PBMC12k",
                num_cluster = c(8,9,10,12,13,14,15),
                output_path= "./",
                host = "lsf",
                queue = "standard",
                output_path = "./")
```



## Step 3: Cell type annotation after clustering
{: .d-inline-block :}

### Read MICA output
{: no_toc }

After clustering via MICA, with function `readMICAoutput`, you can load MICA output (in .txt) as well as your input expression matrix in R to an `expressionSet` object. `ExpressionSet` is the major data structure we used for downstream analysis in R.

> **Note: All functions are designed compatible for both expressionSet and SparseExpressionSet** 

Users shall start with one particular MICA membership and study your optimal number of cluster with cell type signatures. By setting `load_ClusterRes` as TRUE, clustering label will be saved under `eset$ClusterRes`.

```R
eset <- readMICAoutput(eset = eset.log2, load_ClusterRes = TRUE, output_file = "MICA/PBMC12k_k8_tsne_ClusterMem.txt")
```

In order to visualize MICA labels or other metadata on tSNE/UMAP coordinates, your can use function `MICAplot`. Users are required to specify `X` and `Y` coordinates in this function. Other meta data could also be visualized by changing `label` handler. This function outputs a ggplot object. 

```R
MICAplot(input_eset = eset,
		   X = "X", Y="Y", # which meta variable was treated as x or y coordinates
		   color_by = "ClusterRes", pct = 0.5)
```

<center><img src="../plots/3_0_MICA_k8.png" alt="MICA" width="600"></center>

### Marker gene visualization 
{: no_toc }

Picked marker genes could be visualized on t-SNE scatterplot, violin plot or heatmap via function `feature_highlighting`, `feature_vlnplot` and `feature_heatmap`. This will not only help cluster annotation, but identify optimal number of clusters as well.

```R
gn.sel<-c("CD3D","CD27","IL7R","SELL","CCR7","IL32","GZMA",
          "GZMK","DUSP2","CD8A","GZMH","GZMB","CD79A","CD79B","CD86","CD14")

p <- feature_highlighting(input_eset = eset, target = gn.sel, 
		feature="geneSymbol",ylabel = "log2Exp", x="X",y="Y",pct.size = 0.5)
```

<center><img src="../plots/3_1_gene_highlighting.png" alt="Scatterplot" width="800"></center>

```R
p <- feature_vlnplot(input_eset=eset,target=gn.sel,feature = "geneSymbol",
	group_by = "ClusterRes",ylabel = "log2Exp",ncol = 4)
```
<center><img src="../plots/3_2_gene_highlighting_vlnplot.png" alt="violinplot" width="800"></center>


```R
feature_heatmap(input_eset = eset, target = gn.sel, group_by = "ClusterRes",
			 	save_plot = TRUE,width = 6,height = 6,
             name = "log2_expression",plot_name="./GeneHeatmap.png")
```
<center><img src="../plots/3_3_Marker_heatmaps.png" alt="heatmaps" width="800"></center>



### Assign cell type to cluster
{: no_toc }

In order to help assign cell types to each cluster in a more systemmatic way, we introduced `marker_bbplot` function. This function calculated cell type scores for each clusters, and visualize scores using bubble plot, with color scale indicates marker score while circle(bubble) stand for sizes. However, this fucntion requires a pre-defined marker gene lists as input, here we curated a list of well-known marker genes of 9 common immune celltypes as `ref`. **Users are required to follow below header format in order to run this function properly**.

```R
ref <- read.xlsx("Immune_signatures.xlsx")
head(ref)
> head(ref)
  celltype markers weight
1   NaiveT    SELL      1
2   NaiveT    CCR7      1
3     Tmem    IL7R      1
4     Tmem    CD27      1
5     Tmem    IL32      1
6     Tmem    GZMA     -1

draw.marker.bbp(ref=ref,eset=eset,width = 6,height=4, feature = "geneSymbol",group_name = "ClusterRes",
                   save_plot = TRUE, plot_name = "plots/MICA_cluster_score.png")

```
<center><img src="../plots/3_4_MICA_cluster_score.png" alt="Markerbbp" width="800"/></center>



Before diving into network generation section, please assign your celltype as factors in your expression set. Assign cell type as factor will help keep current order of your clusters in following visualizations. **Please do not include "_" in your cell type names since it will cause mis-parsing in later analysis**.

```R
indx<-factor(x=c("NaiveT","Tmem","CD8em","CD8eff","NK","Bcell","DC","Mo"),
             levels=c("NaiveT","Tmem","CD8em","CD8eff","NK","Bcell","DC","Mo"))
eset$celltype <- indx[eset$ClusterRes]

```


## Step 4: Network generation via SJARACNe
{: no_toc }

### Generate SJARACNe input
Before generating cell type/group/cluster specific network, group information should be stored under `pData([your_expressionSet])`. An R function called `generateSJAracneInput` will help to partition input expression matrix and perform essential filtering(filter out 0 expressed genes in cluster) to ensure a reliable network construction. `funcType` is required to specify what kind of network to generate. If set `funcType="TF"`, a reference transcription factor list will be loaded automatically without manual input. However, you do need to define species information for your data using under `ref`.

This function will help create one directory for each group/cell type, containing required input for SJARACNe such as filtered expression matrix in .exp format and filtered TF list in .txt format. 

```R
generateSJARACNeInput(
	input_eset = eset, funcType = "TF", 
	ref = "hg",  #human
	wd.src = "SJARACNE",  #Output directory
	group_tag = "celltype")
```

### Run SJARACNe

SJARACNe works as a separate module which was implemented in python, please refer to [SJARACNe](https://github.com/jyyulab/SJARACNe) github page for installation and basic usage. Please save your working directory before running SJARACNe. 

Here we provide an example to run SJARACNe for all celltypes/clusters. After SJARACNe was sucessfully completed, you will be able to get one network for each cell and functional type.

```
indir = ~/PBMC12K/SJARACNE_PBMC12K/

for i in $(ls -d */ | cut -f1 -d'/');do
sjaracne ${i} $indir/${i}/*.exp $indir/${i}/tf/*.txt $indir/${i}/tf/ --c_threshold 0.01;
echo ${i} has been submitted!;
done
```


## Step 5: Identify cell type specific hidden driver
{: no_toc }

Identify hidden driver from content-based network is the key step in scMINER to help understand your scRNA-seq data, and provide biological insight. 

### Calculate activity
{: no_toc }
Activity calculation is the basis of driver estimation in scMINER. To infer driver activity, expression profile of their targets are intergrated via function `GetActivityFromSJARACNe`. This function takes SJARACNe output path and expression set as input, and return an activity set as well as structured network files if set `save_network_files=TRUE`. **Please note that this function could only work if you used `generateSJARACNeInput` to create SJARACNe input directory and files.**

Since scRNA-seq data are extremly sparse and noisy, please set `activity.method` as `'unweighted'`. 

```R
acs.12k <- GetActivityFromSJARACNe(
    SJARACNe_output_path ="SJARACNE/",
    SJARACNe_input_eset = eset,
    activity.method="unweighted", # we highly recommend using 'unweighted' as activity calculation method
    activity.norm=TRUE, 
    group_tag = "celltype", # which group was used to partition expression profiles
    save_network_file=TRUE, # whether or not save network for each group
    save_path="./networks/") #default as false, but recommended to be TRUE
```


### Driver estimation by differential activity analysis
{: no_toc }

The function `get.DA` was designed to perform differnetial activity analysis from SJARACNe inferred activity matrix. In this function, two-sided student's t-test will be performed to compare mean activity from one cell type V.S. the others. It will return a data frame that includes all TF occurred in original data. Statistics such as t.statistics, p-value, 95%CI, etc. are outputed to help identify hidden drivers. You can save it to file in order to check them manually. 

```R
DAG_result <- get.DA(input_eset = acs.12k,group_tag = "celltype")
```

We also offer a function called `get.Topdrivers` to help picking top drivers for each cell type. You can specify `n` as maximum number of top drivers to pick, and `degree_filter` to restrict number of targets. 

```R
TF_list <-get.Topdrivers(DAG_result = DAG_result,
                             celltype = levels(acs.12k$celltype), # ensure cluster order
                             n = 5, degree_filter = c(50,600))
```


In scMINER, we provide a handful of visualizations to compare driver activity from different cell type/clusters. Here we demo two basic functions: `feature_heatmap` and `feature_vlnplot`. These functions could be used on either expression and activty matrix.

```R
feature_heatmap(input_eset = acs.12k, target = TF_list, group_tag = "celltype",feature = "geneSymbol",
                width = 6,height = 6, save_plot=TRUE, cluster_rows = FALSE,
                name = "Activity",plot_name="./21_TopTFHeatmap.png")
```
<center><img src="../plots/4_1_TopTFHeatmap.png" alt="Driver heatmap" width="550"></center>

```R
#check postive controls
p<-feature_vlnplot(input_eset = acs.12k,feature = "geneSymbol",target=c("LEF1","TCF7","BATF","TBX21","IRF8","SPIB","BATF3","CEBPA"),
                    ylabel = "Activity",group_by = "celltype", ncol=2)
```
<center><img src="../plots/4_2_Known_MR_vlnplot.png" alt="Driver heatmap" width="600"></center>


In order to perform more advanced network analysis utilizing SJARACNe generated cell type specific networks, please infer detailed guidance under [Advanced analysis](../Advanced analysis/PBMC-12k-network) tab.



---

## R session Info
```R
> sessionInfo()
```