


#################################################################################################
#                                                                                               #
#   Load annotation files for TCGA/TARGET data and Cell lines from CMP                          #
#                                                                                               #
#                                                                                               #
#################################################################################################

cmp<-read.csv(paste0(inputdata,"model_list_latest.csv"),header=T,stringsAsFactors = F)
cmp2<-cmp
cmp2$model_id<-cmp2$BROAD_ID
MASTER_LIST<-rbind(cmp,cmp2)
AllMut<-readRDS(file=paste0(inputdata,"ConsMut.Rds"))
colnames(AllMut)<-gsub(".","-",colnames(AllMut),fixed=T)
load(paste0(inputdata,"metadata.rdata"))

CIannot<-read.csv(paste0(inputdata,"Celligner_info.csv"),header=T,stringsAsFactors = F)
TCGACMPmap<-read.csv(paste0(inputdata,"TCGA_CMP_TissueMap.csv"),header=T,stringsAsFactors = F)




#################################################################################################
#                                                                                               #
#   Plots of results from Clustering                                                            #
#                                                                                               #
#                                                                                               #
#################################################################################################
load(file=paste0(inputdata,"/umapLocs.rdata"))
rownames(umapLocs)<-make.names(rownames(umapLocs))
split.by=NULL
shape.by=NULL
dims=1:2
pt.size=2
umapLocs<-data.frame(cbind(umapLocs,meta.data[rownames(umapLocs),c("lineage","type")]),stringsAsFactors = FALSE)
expMat<-matrix(0,4,nrow(umapLocs),dimnames=list(letters[1:4],rownames(umapLocs)))

seu_obj <- Seurat::CreateSeuratObject(expMat,
                                      min.cells = 0,
                                      min.features = 0,
                                      meta.data = umapLocs,
                                      assay="RNA")
seuDR<- CreateDimReducObject(embeddings = as.matrix(umapLocs[,1:2]),key="UMAP_")
seu_obj@reductions[["umap"]]<-seuDR
Idents(seu_obj)<-"lineage"
plotdp<-Seurat::DimPlot(seu_obj, reduction = 'umap', group.by = 'lineage',shape.by="type", pt.size = 0.25,cells.highlight=rownames(umapLocs)[umapLocs$type=="CL"],sizes.highlight = 0.5) + ggplot2::theme(legend.position = 'none')+
                 theme(axis.text.x = element_blank(),axis.text.y=element_blank())

load(paste0(inputdata,"TissueColours.Rdata"))

pdf(paste0(outputdata,"SeuratPlotBC_highlightCL.pdf"),useDingbats = FALSE)
print(plotdp)
dev.off()

plotdp<-Seurat::DimPlot(seu_obj, reduction = 'umap', group.by = 'lineage',shape.by="type", pt.size = 0.25) + ggplot2::theme(legend.position = 'none')+
  theme(axis.text.x = element_blank(),axis.text.y=element_blank())



pdf(paste0(outputdata,"SeuratPlotBC_LineageColL.pdf"),useDingbats = FALSE)
print(plotdp)
dev.off()

#################################################################################################
#                                                                                               #
#   Create feature labels for pan cancer clusters                                               #
#                                                                                               #
#                                                                                               #
#################################################################################################


mdata<-meta.data[meta.data$type=="CL",]
mdata<-mdata[!is.na(mdata$type),]
PancanCelligner<-paste0("PCcelligner",unique(meta.data$cluster),"_Expr")
PCmat<-matrix(0,nrow=length(PancanCelligner),ncol=length(unique(mdata$sampleID)))
dimnames(PCmat)<-list(PancanCelligner,unique(mdata$sampleID))
meta.data$clusterlabel<-paste0("PCcelligner",meta.data$cluster,"_Expr")


for(i in 1:nrow(mdata)){
  # print(mdata[i,"sampleID"])
  
    PCmat[mdata[i,"cluster"],mdata[i,"sampleID"]]<-1

  
}

PCmat0<-PCmat
temp<-CLnameMapping(AllMut,PCmat,"col",annotation=MASTER_LIST)
PCmat<-temp$refdata
colnames(PCmat)<-gsub(".","-",colnames(PCmat),fixed=T)

save(PCmat,file=paste0(outputdata,'PCmatCons2.rdata'))

#################################################################################################
#                                                                                               #
#   Check how many clusters can be tested with the available CRISPR data                        #
#                                                                                               #
#                                                                                               #
#################################################################################################


load(paste0(inputdata,"manifestPostCLGROUP.Rdata"))

set1<-unique(manifest$INSTITUTE_ID)
Bid<-cmp[match(set1,cmp$model_id),"BROAD_ID"]
Sid<-cmp[match(set1,cmp$BROAD_ID),"model_id"]
allModels<-setdiff(unique(c(set1,Bid,Sid)),NA)
PCpancan<-PCmat[,colnames(PCmat)%in%allModels]

cat(paste("Number models with a Celligner Group:", ncol(PCpancan)),file=paste0(outputdata,"/CellignerAnalysis.txt"),append=F)
cat("\n",file=paste0(outputdata,"/CellignerAnalysis.txt"),append=T)

cat(paste("Number clusters testable for Celligner Groups PANCAN:", sum(rowSums(PCpancan,na.rm=T)>2&rowSums(PCpancan,na.rm=T)<ncol(PCpancan)-2) ),file=paste0(outputdata,"/CellignerAnalysis.txt"),append=T)
cat("\n",file=paste0(outputdata,"/CellignerAnalysis.txt"),append=T)

#################################################################################################
#                                                                                               #
#   Calculate the number of cell lines matching tumor lineage                                   #
#                                                                                               #
#                                                                                               #
#################################################################################################


load(paste0(inputdata,"CellignerData.Rdata"))

meta.data$sampleID<-rownames(meta.data)
colnames(tCI)<-make.names(colnames(tCI))
tCI<-as.matrix(tCI)

#load in precomputed as can take time to calculate
#tumor_CL_cor<-calc_tumor_CL_cor(tCI,meta.data)
#save(tumor_CL_cor,file="/Volumes/GoogleDrive/My Drive/Priorityv2//tumor_CL_cor.Rdata")
load(paste0(inputdata,"tumor_CL_cor.Rdata"))
unknownT<-meta.data[meta.data$lineage=="unknown","sampleID"]
tumor_CL_cor<-tumor_CL_cor[setdiff(rownames(tumor_CL_cor),unknownT),]
tumor_CL_cor<-tumor_CL_cor[,setdiff(colnames(tumor_CL_cor),unknownT)]


cl_tumour_classes<-get_cell_line_tumor_class(tumor_CL_cor,meta.data,removeUnknown = FALSE)

classfreq<-cell_line_tumor_class_plot(cl_tumour_classes,meta.data,tumor_CL_cor,paste0(outputdata,"/CellignerMatch.pdf"))


