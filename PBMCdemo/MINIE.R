#Wrapped up pipeline with two input:
# MICA input expression matrix
# MICA output txt file
source("MINIE_FUN.r")

eset.demo<-readMICAoutput(input_file="PBMC_Demo_MICA_input_mini.txt", load_clust_label=TRUE,
			   output_file="scMINER_PBMCdemo/scMINER_MICA/scMINER_PBMCdemo_MDS_4/scMINER_MICA_out/PBMCdemo.ggplot.txt")

# highlighitng causal marker genes
gn.sel <- c("GZMK","GZMH","GZMA","CCR7","CD8A","SELL")
p <- gene_highlighting(input_eset=eset.demo, target = gn.sel, title.size = 8)
p
ggsave(p,filename = "TcellMakershighlight.png",width=7,height=3.5,units = "in",dpi=300)

p<-gene_vlnplot(eset.demo,target=gn.sel,group_tag = "label")
p
ggsave(p,filename = "TcellMakersVlnPlot.png",width=7,height=3.5,units = "in",dpi=300)

hmp<-gene_heatmap(eset = eset.demo,target = gn.sel,group_tag = "label",save_plot = TRUE,width = 6,height = 6,save_plot=TRUE,
                  name = "log2_expression",plot_name="./GeneHeatmap.png")
draw(hmp)

# Determine cell type
# Curate a list of good markers in .RObj
# Assign cluster label automatically (Activity based)

ref<-read.xlsx("../../Immunology references/TF_Immune_cells.xlsx",sheet = "Markers_hg19")
hmp<-AssignCellTypes.Hmp(ref=ref,eset=eset.demo,save_plot = TRUE)

celltype<-c("MemoryT","NaiveT","CD8em","CD8eff")
eset.demo$celltype <- celltype[eset.demo$label]

# Generate SJARACNE input from eset 
tf_sigs.hg <-openxlsx::read.xlsx("../Ref/TF_SIG/human/tf_sigs.updated.201806.xlsx")
tf.ref <- dplyr::filter(tf_sigs.hg,grepl('TF',funcType))$geneSymbol

generateSJAracneInput(eset=eset.demo,tf.ref=tf.ref,wd.src="Sjaracne/",group_tag="celltype")

##################
###need editing
acs.demo<-GetActivityFromSJARANCE(SJaracne_output_path="Sjaracne/",SJaracne_input_eset=eset.demo,
						  activity.method="unweighted",activity.norm=TRUE,group_tag = "celltype",
						  save_network_file=FALSE, save_path=NA)

# find a way to Define highly variable Master regulators
res <- FindHVG(eset = acs.demo,group_tag = "celltype",print_screen = FALSE)
#gene_vlnplot(eset=acs.demo,target = resTop10,group_tag = "celltype",ylabel = "Activity",title.size = 5,ncol = 3)
write.xlsx(res,file="AnovaTest.xlsx")

# t-test 
res <- FindDAG(eset=acs.demo,group_tag = "celltype")

