---
layout: default
title: Advanced analysis
nav_order: 6

---

# Advanced analysis
{:.no_toc}
Here we demonstrate our advanced downstream analysis pipeline using PBMC (10x genmomics) scRNA-seq data after following driver estimation tutorial under [step by step user guide](./PBMC-12k.md).

## Network visualization
scMINER incorporates a handful of network visualization/exploration function adapted from [NetBID2](https://jyyulab.github.io/NetBID/), a powerful tool for data-driven network-based bayesian Inference of drivers. scMINER also offered several wrappers of basic visualization functions in NetBID2 for better usability. 

### Load networks
You can retrieve your network and store them in a network structure with function `NetBID2::get.SJAracne.network`. 
```R
net1 <- NetBID2::get.SJAracne.network('SJARACNE/NaiveT_8469_8468_4421/tf/SJARACNE_NaiveT_8469_8468_4421/SJARACNE_out.final/consensus_network_ncol_.txt')
net2 <- NetBID2::get.SJAracne.network('SJARACNE/Tmem_8489_8488_2482/tf/SJARACNE_Tmem_8489_8488_2482/SJARACNE_out.final/consensus_network_ncol_.txt')
```

Or, if you followed our analysis pipeline under [step by step user guide](./PBMC-12k.md), you should be able to load 
your network files under `./networks` folder:
```R
load("../plots/NaiveT_TF.network")
```

### Single network visualization
In scMINER, you can visualize your driver and its targets by function `draw.network`. It was adapted from function `draw.targetNet` and `draw.targetNet.TWO` from `NetBID2`. This function can help visualize a driver's targets as well as the relationship(edge) between source and target genes, by taking Mutual information as edge weight, and spearman correlation as direction.

```R
draw.network(net1 = net,src1 = "LEF1", #driver name
	ifSymbol = TRUE, #if your source driver name is a Symbol 
	ifWeighted = TRUE, #if plot edge with weight and direction
	pdf_name = "LEF1.TF_in_NaiveT_network.pdf",
	n_layer=4)

```
<center><img src="../plots/LEF1.TF_in_NaiveT_network.png" alt="drawing" width="600"></center>


### Subnetwork structure visualization between two networks
You can also use `draw.network` function to visualize two networks and their subnetwork structure. This could be used for: 
- Identify common targets from two top driver from the same network
- Identify network rewiring event of same driver in different cell type network.   


Here below is an example for later case:

```R
draw.network(net1 = net1,net2 = net2,
		src1 = "BATF",src2="BATF", 
		source1_z=-3, source2_z=4,
		ifSymbol = TRUE,ifWeighted = TRUE,
		pdf_name = "BATF.TF_in_2_network.pdf",
		n_layer=4)
```
<center><img src="../plots/BATF.TF_in_2_network.png" alt="drawing" width="600"></center>



## Biological function anlaysis for drivers 

### Gene set overlap with targets visualized by bubble plot
When picking candidate hidden drivers, it would be extremly helpful if we could identify the potential biological pathways this driver regulates. With SJARACNe infered network, we can assess as well as  visualizethe overlap between knowledge-based gene sets and driver's targets via function `draw.bubblePlot`. This function returns a bubble plot indiating results from Fisher exact test. 

Before using`draw.bubblePlot` function, you have to load genesets in your working environment by function `gs.preload()`

```R
gs.preload(use_spe='Homo sapiens',update=FALSE)
```

Then you can use function `TopMasterRegulator` to pull out top hidden driver candidates from your Differential activity analysis results. Or write your own fucntions to hand pick candidate to visualize. Here we provide an example of using function `TopMasterRegulator`.

```R 
TF_list <- TopMasterRegulator(DAG_result = DAG_result,
                              celltype="NaiveT",
                              n = 10, degree_filter = c(50,800))
```
Next generate your ID to symbol conversion table, since all gene sets are curated at gene symbol level. Here in our data, we used ensembl_id as our default ID for network construction. In order to match your ID with gene symbols, you can use function: 

```R
tbl <- get_IDtransfer2symbol2type(from_type = "ensembl_gene_id",
		use_genes = rownames(eset.12k),ignore_version = TRUE)
```

Here we provide an example for ploting overlap between target list of top drivers in Naive T cells, and knowledge based gene sets from "Hallmark","KEGG" and "GO".

```R
draw.bubblePlot(driver_list = TF_list,
                show_label = DAG_result[TF_list,"geneSymbol"],
                Z_val = DAG_result[TF_list,"Z_NaiveT"],
                driver_type = NULL,
                target_list = net1$target_list,
                transfer2symbol2type = tbl,
                bg_list = fData(eset.12k)$geneSymbol,
                min_gs_size = 50, max_gs_size = 600, 
                top_geneset_number=8,top_driver_number=10,use_gs = c("H","C5","CP:KEGG"),
                pdf_file = 'NaiveT_bubblePlot.pdf',
                main ='Bubbleplot for top driver targets in NaiveT')

```
<center><img src="../plots/NaiveT_bubblePlot.png" alt="gsbbp"></center>



---
```
> sessionInfo()
```
---