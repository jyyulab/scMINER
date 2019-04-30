load("~/Box Sync/Misc/PBMC/latest/pbmc12k_eset")
tmp<-exprs(eset.12k)
tmp.smtx<-Matrix(tmp,sparse = TRUE)
fd<-fData(eset.12k)
pd<-pData(eset.12k)


library(Biobase)
setClass( "scMINEReSet",
          contains = "ExpressionSet",
          slots = c(raw.count.data="matrix"),
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c(classVersion("ExpressionSet"), scMINEReSet = "1.0.0" )))
)

cds <- new( "scMINER-Eset",
            assayData = assayDataNew( "environment", exprs=tmp.smtx),
            phenoData= new("AnnotatedDataFrame",data=pd),
            featureData= new("AnnotatedDataFrame",data=fd))


