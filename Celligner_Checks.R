


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
#   EMT status of transcriptional profiles                                                      #
#                                                                                               #
#                                                                                               #
#################################################################################################



load(file=paste0(inputdata,"/AllEMT.Rdata"))
names(AllEMT)<-make.names(names(AllEMT))
EMTpdata<-umapLocs
EMTpdata$EMTscore<-AllEMT[rownames(EMTpdata)]
EMTbreaks<-seq(-2,2,1)
EMTpdata$EMTgroup<-1
EMTpdata[EMTpdata$EMTscore<2,"EMTgroup"]<-2
EMTpdata[EMTpdata$EMTscore<1,"EMTgroup"]<-3
EMTpdata[EMTpdata$EMTscore<0,"EMTgroup"]<-4
EMTpdata[EMTpdata$EMTscore<(-1),"EMTgroup"]<-5
EMTpdata[EMTpdata$EMTscore<(-2),"EMTgroup"]<-6
EMTpdata$cluster<-meta.data[rownames(EMTpdata),"cluster"]
EMTpdata<-EMTpdata[!is.na(EMTpdata$cluster),]

pdf(paste0(outputdata,"Umap_type.pdf"),useDingbats = FALSE)
print(ggplot(EMTpdata,aes(x=UMAP_1,y=UMAP_2,color=type))+geom_point(aes(color=type),size=0.2,show.legend = TRUE)+scale_colour_manual(values=c("CL"="red","tumour"="grey"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")))
dev.off()

pdf(paste0(outputdata,"Umap_lineage.pdf"),useDingbats = FALSE,width=10,height=5)
print(ggplot(EMTpdata,aes(x=UMAP_1,y=UMAP_2,color=lineage))+geom_point(aes(color=lineage),size=0.2,show.legend = TRUE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")))
dev.off()

pdf(paste0(outputdata,"Umap_EMT.pdf"),useDingbats = FALSE)
print(ggplot(EMTpdata,aes(x=UMAP_1,y=UMAP_2,color=EMTscore))+geom_point(aes(color=EMTscore),size=0.2)+scale_colour_viridis_c()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")))
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



#load in precomputed as can take time to calculate
#tumor_CL_cor<-calc_tumor_CL_cor(tCI,meta.data)
#save(tumor_CL_cor,file="/Volumes/GoogleDrive/My Drive/Priorityv2//tumor_CL_cor.Rdata")
load(paste0(inputdata,"tumor_CL_cor.Rdata"))
meta.data$lineage<-tolower(meta.data$lineage)
unknownT<-meta.data[meta.data$lineage=="unknown","sampleID"]
tumor_CL_cor<-tumor_CL_cor[setdiff(rownames(tumor_CL_cor),unknownT),]
tumor_CL_cor<-tumor_CL_cor[,setdiff(colnames(tumor_CL_cor),unknownT)]
rownames(meta.data)<-meta.data$sampleID
rownames(meta.data)<-make.names(rownames(meta.data))
meta.data$sampleID<-rownames(meta.data)
meta.emt<-meta.data
meta.emt$lineage<-NA

meta.emt[rownames(EMTpdata),"lineage"]<-EMTpdata$EMTgroup

cl_tumour_classesE<-get_cell_line_tumor_class(tumor_CL_cor,meta.emt)

classfreqE<-cell_line_tumor_class_plot(cl_tumour_classesE,meta.emt,tumor_CL_cor,paste0(outputdata,"/CellignerMatch_EMT.pdf"))
CEm<-as.matrix(classfreqE)
diag(CEm)
cl_tumour_classes<-get_cell_line_tumor_class(tumor_CL_cor,meta.data,removeUnknown = FALSE)

classfreq<-cell_line_tumor_class_plot(cl_tumour_classes,meta.data,tumor_CL_cor,paste0(outputdata,"/CellignerMatch.pdf"))
cfmat<-as.matrix(classfreq)
cfval<-diag(cfmat)

all.equal(names(cl_tumour_classes),names(cl_tumour_classesE))
Assignments<-cbind(cl_tumour_classes,cl_tumour_classesE)
Assignments<-cbind(Assignments,meta.data[rownames(Assignments),"lineage"])
Assignments<-cbind(Assignments,EMTpdata[rownames(Assignments),'EMTgroup'])
matchLineage<-Assignments[,1]==Assignments[,3]
matchEMT<-Assignments[,2]==Assignments[,4]
Assignments<-cbind(Assignments,matchLineage,matchEMT)
matchEither<-as.logical(Assignments[,5])|as.logical(Assignments[,6])
Assignments<-cbind(Assignments,matchEither)


sum(matchEither)/nrow(Assignments)
