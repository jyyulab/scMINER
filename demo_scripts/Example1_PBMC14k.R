library(scMINER)

############### PART-I: Load in dataset & Run MICA ###############

### Step0: Create directory
project_main_dir <- '../test'
project_name <- 'PBMC14KDS'
scminer.par <- scMINER::scMINER.dir.create(project_main_dir = project_main_dir,
                                           project_name = project_name)

### Step1: Load in dataset
## Option1: load data from standard 10x genomics files
demo_dir <- system.file('PBMC14KDS_DemoDataSet/DATA/10X/',package = "scMINER")
pbmc.14k.DS.eset <- scMINER::readscRNAseqData(file = demo_dir,
                                              is.10x = T,
                                              CreateSparseEset = T,
                                              add.meta = T)

## Option2: load data from plain matrix text file
# demo_file <- ''
# pbmc.14k.DS <- readscRNAseqData(file=demo_file,
#                                is.10x = F,
#                                CreateSparseEset = F,
#                                add.meta = F,
#                                sep='\t')
# pbmc.14k.DS.eset <- scMINER::CreateSparseEset(data = t(as.matrix(pbmc.14k.DS)),
#                                              add.meta = T)

## Option3: load data from Seurat object file (seruat --> SparseEset)
# meta.data<-SeuratObject@meta.data
# feature.data<-data.frame(rownames(SeuratObject@assays$RNA@data))
# colnames(feature.data)<-"geneSymbol"
# rownames(feature.data)<-feature.data$geneSymbol
# eset<-CreateSparseEset(data=SeuratObject@assays$RNA@data,meta.data = meta.data,feature.data = feature.data,add.meta = F)

### Step2: QC for the SparseEset and filtration
cutoffs <- scMINER::draw.scRNAseq.QC(SparseEset = pbmc.14k.DS.eset,
                            project.name = project_name,
                            plot.dir = scminer.par$out.dir.QC,
                            group = "group",      # this indicate which meta data information will be use in x axis to group violin plots
                            output.cutoff = TRUE) # whether or not to output suggested cutoffs

## Filter for the SparseEset
pbmc.14k.DS.eset.filter <- scMINER::preMICA.filtering(SparseEset = pbmc.14k.DS.eset,
                                             cutoffs = cutoffs)
norm_value <- 1e6
exp.norm <- sweep(exprs(pbmc.14k.DS.eset.filter), 2,
                  norm_value/unname(Matrix::colSums(exprs(pbmc.14k.DS.eset.filter))), '*') # normalization
exp.log2 <- log(exp.norm + 1, base = 2) # log transformation (required by MICA)

## save as SparseEset
pbmc.14k.DS.eset.log2 <- CreateSparseEset(data = exp.log2,
                                       meta.data = pData(pbmc.14k.DS.eset.filter),
                                       feature.data = fData(pbmc.14k.DS.eset.filter),
                                       add.meta = F)
scminer.par$out.dir.DATA_eset <- sprintf('%s/pbmc.14k.DS.eset.log2.RData',
                                         scminer.par$out.dir.DATA)
save(pbmc.14k.DS.eset.log2,file=scminer.par$out.dir.DATA_eset)

### Step3: Prepare MICA input and run MICA
scminer.par$out.dir.MICA_input <- sprintf('%s/PBMC14KDS_MICA_input.h5ad',
                                          scminer.par$out.dir.MICA)
MICA.cmd <- generateMICAinput(eset = pbmc.14k.DS.eset.log2 ,
                              filepath = scminer.par$out.dir.MICA_input,
                              scminer.par = scminer.par)
scminer.par$MICA.cmd <- MICA.cmd

### Step4: Run MICA: mica mds -i ../test/PBMC14KDS/MICA//PBMC14KDS_MICA_input.h5ad -o ../test/PBMC14KDS/MICA/ -pn PBMC14KDS -nc 4 5 6 7 8 9 10 -dk 19
# system(scminer.par$MICA.cmd)

### Save scminer.par
save(scminer.par,file=scminer.par$out.dir.DATA_par)

############### PART-II: Cell type annotation based on MICA output ###############
### Step0: load scminer.par into R session:
# e.g load("../test/PBMC14KDS/DATA/scminer.par.RData")

### Step1: Load MICA output clustering results
## load SparseEset RData
# load(scminer.par$out.dir.DATA_eset) ## if re-start R-session
## or load from R package
# demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#                         package = "scMINER")
# load(demo_file)

## read in MICA output file (choosing a clustering result based on known cell type signatures or silhouette scores), the file name can be modified
scminer.par$out.dir.MICA_output <- sprintf('%s/clustering_umap_euclidean_19.txt',
                                           scminer.par$out.dir.MICA)

## or use file from R package
# scminer.par$out.dir.MICA_output <- system.file('PBMC14KDS_DemoDataSet/MICA/clustering_umap_euclidean_19.txt',
#             package = "scMINER")

pbmc.14k.DS.eset.log2 <- scMINER::readMICAoutput(eset = pbmc.14k.DS.eset.log2,
                                                 load_ClusterRes = TRUE,
                                output_file = scminer.par$out.dir.MICA_output)

## draw MICA results
scMINER::MICAplot(input_eset = pbmc.14k.DS.eset.log2, X = "X", Y = "Y",
                  color_by = "ClusterRes", pct = 0.5)

### Step2: Cell type annotation and visualization
genes_of_interest <-c("CD3D", "CD27", "IL7R",
                      "SELL", "CCR7", "IL32",
                      "GZMA", "GZMK", "DUSP2",
                      "CD8A", "GZMH", "GZMB",
                      "CD79A", "CD79B", "CD86", "CD14")
## UMAP scatter plot
scMINER::feature_highlighting(input_eset = pbmc.14k.DS.eset.log2,
                              target = genes_of_interest,
                              feature = "geneSymbol",
                              ylabel = "log2Exp",
                              x = "X", y = "Y",
                              pct.size = 0.5)

## Violoin plot
scMINER::feature_vlnplot(input_eset = pbmc.14k.DS.eset.log2,
                         target = genes_of_interest,
                         feature = "geneSymbol",
                         group_by = "ClusterRes",
                         ylabel = "log2Exp", ncol = 4)

## Heatmap plot
scMINER::feature_heatmap(input_eset = pbmc.14k.DS.eset.log2,
                         target = genes_of_interest,
                         group_name = "ClusterRes",
                         save_plot = FALSE,
                         width = 6, height = 6,
                         name = "log2Exp")

## Cell type annotation
## Draw a bubble plot
markers_file <- system.file('PBMC14KDS_DemoDataSet/DATA/',
                  'Immune_signatures.xlsx',
                  package = "scMINER")
markers <- openxlsx::read.xlsx(markers_file)
draw.marker.bbp(ref = markers, input_eset = pbmc.14k.DS.eset.log2,
                width = 6, height = 4, feature = "geneSymbol",
                group_name = "ClusterRes", save_plot = FALSE)

indx <- factor(x=c("Monocyte", "CD4TCM", "NK", "Bcell"),
               levels=c("Monocyte", "CD4TCM", "NK", "Bcell"))
pbmc.14k.DS.eset.log2$celltype <- indx[pbmc.14k.DS.eset.log2$ClusterRes]

## save updated SparseEset
save(pbmc.14k.DS.eset.log2,file=scminer.par$out.dir.DATA_eset)

### Save scminer.par
save(scminer.par,file=scminer.par$out.dir.DATA_par)

############### PART-III: Network generation via SJARACNe ###############
### Step0: load scminer.par into R session:
# e.g load("../test/PBMC14KDS/DATA/scminer.par.RData")
## load SparseEset RData
# load(scminer.par$out.dir.DATA_eset) ## if re-start R-session
## or load from R package
# demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#                               package = "scMINER")
# load(demo_file)

## Step1: Generate SJARACNe input
SJAR.cmd.tf <- generateSJARACNeInput(
  input_eset = pbmc.14k.DS.eset.log2,
  funcType = "TF",
  ref = "hg",  # human
  wd.src = scminer.par$out.dir.SJAR,  # output directory
  group_name = "celltype")
SJAR.cmd.sig <- generateSJARACNeInput(
  input_eset = pbmc.14k.DS.eset.log2,
  funcType = "SIG",
  ref = "hg",  # human
  wd.src = scminer.par$out.dir.SJAR,  # output directory
  group_name = "celltype")
scminer.par$SJAR.cmd.tf <- SJAR.cmd.tf
scminer.par$SJAR.cmd.sig <- SJAR.cmd.sig
scminer.par$SJAR.group.name <- "celltype"

## Step2: run SJARACNe
# system(scminer.par$SJAR.cmd.tf$local)
# system(scminer.par$SJAR.cmd.sig$local)

### Save scminer.par
save(scminer.par,file=scminer.par$out.dir.DATA_par)

############### PART-IV: Identify cell-type-specific hidden drivers via MINIE ###############
### Step0: load scminer.par into R session:
# e.g load("../test/PBMC14KDS/DATA/scminer.par.RData")
## load SparseEset RData
# load(scminer.par$out.dir.DATA_eset) ## if re-start R-session
## or load from R package
# demo_file <- system.file('PBMC14KDS_DemoDataSet/DATA/pbmc.14k.DS.eset.log2.RData',
#                         package = "scMINER")
# load(demo_file)
## set the out.dir.SJAR to the demo output file
# scminer.par$out.dir.SJAR <- system.file('PBMC14KDS_DemoDataSet/SJAR/',
#                                              package = "scMINER")

### Step1: Calculate activity from SJARACNe output
acs.14k.tf <- GetActivityFromSJARACNe(
  SJARACNe_output_path = scminer.par$out.dir.SJAR,
  SJARACNe_input_eset = pbmc.14k.DS.eset.log2,
  activity.method="unweighted", # we highly recommend using 'unweighted' as activity calculation method
  activity.norm=TRUE,
  group_name = scminer.par$SJAR.group.name, # which group was used to partition expression profiles
  save_network_file=TRUE, # whether or not save network for each group
  functype="tf",
  save_path=scminer.par$out.dir.SJAR)

acs.14k.sig <- GetActivityFromSJARACNe(
  SJARACNe_output_path = scminer.par$out.dir.SJAR,
  SJARACNe_input_eset = pbmc.14k.DS.eset.log2,
  activity.method="unweighted", # we highly recommend using 'unweighted' as activity calculation method
  activity.norm=TRUE,
  group_name = scminer.par$SJAR.group.name, # which group was used to partition expression profiles
  save_network_file=TRUE, # whether or not save network for each group
  functype="sig",
  save_path=scminer.par$out.dir.SJAR)

# save activity to RData file (optional)
scminer.par$out.dir.AC <- sprintf('%s/%s_Activity.RData',scminer.par$out.dir.DATA,
                                  scminer.par$SJAR.group.name)
AC_eset <- list(AC.TF=acs.14k.tf,AC.SIG=acs.14k.sig)
save(AC_eset,file=scminer.par$out.dir.AC)

### Step2: Driver estimation by differential activity analysis
DAG_result_tf <- get.DA(input_eset = acs.14k.tf,
                        group_name = scminer.par$SJAR.group.name)
DAG_result_sig <- get.DA(input_eset = acs.14k.sig,
                         group_name = scminer.par$SJAR.group.name)

celltype <- levels(pData(pbmc.14k.DS.eset.log2)[,scminer.par$SJAR.group.name])
TF_list <- get.Topdrivers(DAG_result = DAG_result_tf,
                          celltype = celltype,
                          # ensure cluster order
                          n = 5, degree_filter = c(50, 600))

### Step3: Check positive controls
p <- feature_vlnplot(input_eset = acs.14k.sig,
                     feature = "fn",
                     target=c("CD27", "IL7R","CCR7",'GZMA','GZMK','DUSP2','GZMH','GZMB'),
                     ylabel = "Activity",
                     group_by = scminer.par$SJAR.group.name, ncol=2)
p

#### Step4: update parameter object
save(scminer.par,file=scminer.par$out.dir.DATA_par)




