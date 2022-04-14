# scMINER Guided Analysis on activity-based clustering

### Generate bulk network by using MetaCell

To overcome the sparcity of scRNA-seq data and facilitate unbiased activity-based clustering, We start generating pseudo-bulk network. 

```R
# Create a MetaCell object with raw count matrix of scRNA-seq (gene x cells)
library(metacell)
if(!dir.exists("testdb")) dir.create("testdb/")
scdb_init("testdb/", force_reinit=T)
mcell_import_scmat_tsv("test", fn="correctPBMC.csv",dset_nm = "pbmc20K")
mat = scdb_mat("test")

```

### Create membership for each cell based on metacell 

```R
## Examine QC, remove unwanted genes and lowly expressed genes ##
if(!dir.exists("figs1")) dir.create("figs1/")
scfigs_init("figs1/")
mcell_plot_umis_per_cell("test")
mat = scdb_mat("test")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("^IGJ", nms, v=T),
             grep("^IGH",nms,v=T),
             grep("^IGK", nms, v=T),
             grep("^IGL", nms, v=T))
bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MT-MR", nms, v=T), grep("^MT-ND", nms, v=T),"NEAT1","TMSB4X", "TMSB10", ig_genes))
mcell_mat_ignore_genes(new_mat_id="test", mat_id="test", bad_genes, reverse=F)
mcell_mat_ignore_small_cells("test", "test",200)

mcell_add_gene_stat(gstat_id="test", mat_id="test", force=T)
mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=100, T_top3=2)
mcell_plot_gstats(gstat_id="test", gset_id="test_feats")

## Clustering into metacells ##
mcell_coclust_from_graph_resamp(
  coc_id="test_coc500",
  graph_id="test_graph",
  min_mc_size=20,
  p_resamp=0.75, n_resamp=500)

mcell_mc_from_coclust_balanced(
  coc_id="test_coc500",
  mat_id= "test",
  mc_id= "test_mc",
  K=30, min_mc_size=30, alpha=2)

mcell_mc_split_filt(new_mc_id="test_mc_f",
                    mc_id="test_mc",
                    mat_id="test",
                    T_lfc=10, plot_mats=F)
                    

```
The membership is stored in mc_f object.

### Generate eset object for network construction
```R
countmat<-read.csv("correctPBMC.csv")
esetNet <- CreateSparseEset(data=countmat,add.meta = T)
norm = 1e6 
exp.norm <- sweep(exprs(esetNet), 2, norm/unname(Matrix::colSums(exprs(esetNet))), '*')
exp.log2 <- log(exp.norm+1,base=2)

eset.log2 <- CreateSparseEset(data=exp.log2, 
                              meta.data = pData(esetNet), 
                              feature.data = fData(esetNet), 
                              add.meta = F)
                              
eset.log2<-eset.log2[,names(mc_f@mc)]
eset.log2$group<-mc_f@mc

```
### Create Pseudobulk eset
```R
countmat<-read.csv("correctPBMC.csv")
esetNet <- CreateSparseEset(data=countmat,add.meta = T)
norm = 1e6 
exp.norm <- sweep(exprs(esetNet), 2, norm/unname(Matrix::colSums(exprs(esetNet))), '*')
exp.log2 <- log(exp.norm+1,base=2)

eset.log2 <- CreateSparseEset(data=exp.log2, 
                              meta.data = pData(esetNet), 
                              feature.data = fData(esetNet), 
                              add.meta = F)
                              
eset.log2<-eset.log2[,names(mc_f@mc)]
eset.log2$group<-mc_f@mc

CreatePseudobulkEset<-function(eset.log2, group){
  cpmdf<-data.frame(t(exprs(eset.log2)))
  cpmdf$grp<-group
  cpmdf2<-aggregate(. ~ grp, data = cpmdf, FUN = mean)
  rownames(cpmdf2)<-cpmdf2$grp
  cpmdf2$grp<-NULL
  cpmdf2<-t(as.matrix(cpmdf2))
  eset<-CreateSparseEset(data=cpmdf2,add.meta = F) 
  return(eset)
}
esetPseudobulkNet<-CreatePseudobulkEset(eset.log2 = eset.log2,group = eset.log2$group)


```

### Generate Network
```R
esetPseudobulkNet$group<-"metacell"
generateSJARACNeInput1(
  eset =esetPseudobulkNet, funcType = "SIG", 
  ref = "hg",  
  wd.src = "SJAR/SJAR_SIG",  #Output directory
  group_tag = "group")

generateSJARACNeInput1(
  eset = esetPseudobulkNet, funcType = "TF", 
  ref = "hg",  
  wd.src = "SJAR/SJAR_TF",  #Output directory
  group_tag = "group")
  
```

