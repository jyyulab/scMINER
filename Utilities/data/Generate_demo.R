# Online github tutorial
# Using extracted T cell data
# Generate truncated demo data

indx.Tcell<- filter(MICA_Result,label%in%c(1,2,3,4))
dim(indx.Tcell) # 9422

ID.500<-sample(indx.Tcell$ID,size=500)

Tcell.500.MICA<-filter(MICA_Result,ID%in%ID.500)

input.500<-input[,ID.500]

save(input.500,file="PBMC_demo_500")

mica.input<-t(input.500)
mica.input<-as.data.frame(mica.input)

mica.input$ID<-colnames(input.500)
mica.input<-mica.input[,c("ID",rownames(input.500))]
write.table(mica.input,file = "PBMC_Demo_MICA_input.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

##########
#only preserve the variable genes 
load("/Volumes/yu3grp/scRNASeq/yu3grp/TracyQian/PBMC/02_MICA/PBMC_12k/Combined_PBMC_12k_hg19_scMINER_input.Rdata")
#Obj.500<-Obj
#Obj.500@data<-Obj@data[,ID.500]
#FindVariableGenes(Obj.500)
#VG <-intersect(Obj.500@var.genes,rownames(input.500))

VG <-intersect(Obj@var.genes,rownames(input.500))

c("GZMA","CCR7","GZMK","GZMH","IL32","SELL")%in%VG #spot checking

input.500.VG <- input.500[VG,]

mica.input<-t(input.500.VG)
mica.input<-as.data.frame(mica.input)

mica.input$ID <- colnames(input.500.VG)
mica.input<-mica.input[,c("ID",rownames(input.500.VG))]
write.table(mica.input,file = "PBMC_Demo_MICA_input_mini.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)



