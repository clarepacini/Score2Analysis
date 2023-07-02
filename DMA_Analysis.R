
#################################################################################################
#                                                                                               #
#   Load annotation data for dependency/feature association (DFA) analysis                      #
#                                                                                               #
#                                                                                               #
#################################################################################################

load(paste0(inputdata,'07_EssMatrix_bDepletionsB2.Rdata'))
cmp<-read.csv(paste0(inputdata,"model_list_latest.csv"),header=T,stringsAsFactors = FALSE)

load(file=paste0(inputdata,"paralogList.Rdata"))
load(file=paste0(inputdata,"cancerDrivers.rdata"))

load(paste0(inputdata,"pcgene.rdata"))

LoFgenes<-cancerDrivers[[1]][which(cancerDrivers[[2]]=="LoF")]

#################################################################################################
#                                                                                               #
#   Load results from regression analyses                                                        #
#                                                                                               #
#                                                                                               #
#################################################################################################



load(file=paste0(inputdata,"/DFRP.Rdata"))
cat(paste("Number tests performed: ",nrow(DFRP)),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=F)


#################################################################################################
#                                                                                               #
#   Generate three groups of DFAs with different levels of stringency                           #
#                                                                                               #
#                                                                                               #
#################################################################################################


DFRa<-DFRP[DFRP$FDR<1&DFRP$f2>0.1&abs(DFRP$logFC)>1,]


colnames(DFRa)[colnames(DFRa)=="Depleted Gene"]<-"TARGET"
colnames(DFRa)[colnames(DFRa)=="depAssoc"]<-"ASSOCIATION_EFFECT"
save(DFRa,file=paste0(outputdata,"/DFRa.Rdata"))
write.table(DFRa,file=paste0(outputdata,"/SupplementaryTable7.tsv"),row.names=F,quote=F,sep="\t")

cat(paste("Number genes with strongest marker: ",length(unique(DFRa$TARGET))),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)

DFRa$M<-unlist(sapply(DFRa$FEATURE,function(x) strsplit(x,"_",fixed=T)[[1]][1]))
DFRa$PPI_min<-unlist(sapply(1:nrow(DFRa),function(x) ifelse(DFRa[x,"M"]==DFRa[x,"TARGET"],0,6)))


DFRe<-DFRP[DFRP$FDR<1&abs(DFRP$logFC)>0.5,]

colnames(DFRe)[colnames(DFRe)=="Depleted Gene"]<-"TARGET"
colnames(DFRe)[colnames(DFRe)=="depAssoc"]<-"ASSOCIATION_EFFECT"
save(DFRe,file=paste0(outputdata,"/DFRe_Extended.Rdata"))
write.table(DFRe,file=paste0(outputdata,"/SupplementaryTable9.tsv"),row.names=F,quote=F,sep="\t")

cat(paste("Number of dependencies with asscociation found permissive: ",length(unique(DFRe$TARGET))),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)

cat(paste("Percent of dependencies with asscociation found permissive: ",length(unique(DFRe$TARGET))/500),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)

NumberMarkers<-table(DFRe$TARGET)
summary(as.vector(NumberMarkers))
DFRe$M<-unlist(sapply(DFRe$FEATURE,function(x) strsplit(x,"_",fixed=T)[[1]][1]))
DFRe$PPI_min<-unlist(sapply(1:nrow(DFRe),function(x) ifelse(DFRe[x,"M"]==DFRe[x,"TARGET"],0,6)))
cat(paste("Median number of feature, permissive set all: ",median(table(DFRe$TARGET))),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)
cat(paste("Minimum number of feature, permissive set all: ",min(table(DFRe$TARGET))),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)
cat(paste("Maximum number of feature, permissive set all: ",max(table(DFRe$TARGET))),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)


#################################################################################################
#                                                                                               #
#   Plot of number of markers per dependency                                                    #
#                                                                                               #
#                                                                                               #
#################################################################################################
markercount<-as.vector(table(DFRe$TARGET))
names(markercount)<-names(table(DFRe$TARGET))
textvals<-markercount[markercount>=200]
plotmc<-data.frame(mc=markercount,label=names(markercount))
plotmc[plotmc$mc<200,"label"]<-""
pointvals<-data.frame(x=textvals,y=0.0001,label=names(textvals))
pdf(paste0(outputdata,"/Permissive_NmarkersPerDep.pdf"),useDingbats = FALSE,width=7,height=4)

ggplot(data=plotmc,aes(x=mc))+geom_density()+ggrepel::geom_label_repel(data=pointvals,aes(x=x,y=y,label=label),max.overlaps = 20,size=3)+theme_bw()+xlab("Number markers")
dev.off()

#################################################################################################
#                                                                                               #
#   Assign DFA groups to significant associations in different sets                             #
#                                                                                               #
#                                                                                               #
#################################################################################################

DFRdataAll<-get_DMA(DFRa,plist,LoFgenes,"FEATURE",cancerDrivers,splitMarkers = FALSE)
DFRdata<-DFRdataAll$DFAresults
useDFAgroups<-c("SelfAddiction_Mutation","Transcriptional_Subtype","SelfAddiction_Expression",
  "Paralog","SL","AddictionND","OncogenicAddiction_Act","SelfAddiction_Variant",
  "SelfAddiction_CN","LoF_mutation_SL","LoF_other_SL","SelfAddiction_Protein","MSI")
cat(paste("Percent of asscociations with DFA: ",sum(DFRdata$SLgroup%in%useDFAgroups)/nrow(DFRdata)),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)
cat(paste("Number of addiction DFA: ",sum(DFRdata$SLgroup%in%c("Addiction","SelfAddiction_Mutation","SelfAddiction_Expression",
                                                               "Addiction","OncogenicAddiction_Act","SelfAddiction_Variant",
                                                               "SelfAddiction_CN","SelfAddiction_Protein"))),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)
cat(paste("Number of DFA: ",sum(DFRdata$SLgroup%in%useDFAgroups)),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)

save(DFRdata,file=paste0(outputdata,"/DFRdata.Rdata"))


DFRsplit<-DFRdataAll$TargetTypeMatrix
DFRallInfo<-cbind(DFRdata,DFRsplit)
write.table(DFRallInfo,file=paste0(outputdata,"/SupplementaryTable8.tsv"),row.names=F,quote=F,sep="\t")

save(DFRsplit,file=paste0(outputdata,"/DFRsplit.Rdata"))
DFRadd<-which(rowSums(DFRsplit[,c("Target_Activating","Target_cancerDriverAct","SelfAddiction_CN","SelfAddiction_Expression","SelfAddiction_Mutation","SelfAddiction_Variant","SelfAddiction_Protein")])>0)

DFRaddD<-which(rowSums(DFRsplit[,c("Target_cancerDriverAct","SelfAddiction_CN","SelfAddiction_Expression","SelfAddiction_Mutation","SelfAddiction_Variant","SelfAddiction_Protein")])>0)

DFRaddND<-which(DFRsplit[,c("Target_Activating")]==1)

DFRaddDA<-which(DFRsplit[,c("Target_cancerDriverAct")]==1)

DFRSelfadd<-which(rowSums(DFRsplit[,c("SelfAddiction_CN","SelfAddiction_Expression","SelfAddiction_Mutation","SelfAddiction_Variant","SelfAddiction_Protein")])>0)

DFRComp<-which(rowSums(DFRsplit[,c("MSI","Tsubtype")])>0)

DFRSL<-which(rowSums(DFRsplit[,c("LoFSLmut","LoFSLother","Paralog","SL")])>0)


vennset<-list()
vennset[["Addiction"]]<-unique(rownames(DFRsplit)[DFRadd])
vennset[["SyntheticLethal"]]<-unique(rownames(DFRsplit)[DFRSL])
vennset[["Composite"]]<-unique(rownames(DFRsplit)[DFRComp])
venn.diagramCP(vennset,paste0(outputdata,"Venn3Groups.pdf"),imagetype="pdf",height=5,width=5,units="in")
vennAdd<-list()
vennAdd[["OtherAddiction"]]<-unique(rownames(DFRsplit)[DFRaddND])
vennAdd[["DriverGeneAddiction"]]<-unique(rownames(DFRsplit)[DFRaddDA])
vennAdd[["SelfAddiction"]]<-unique(rownames(DFRsplit)[DFRSelfadd])
venn.diagramCP(vennAdd,paste0(outputdata,"VennAddiction.pdf"),imagetype="pdf",height=5,width=5,units="in")

vennSL<-list()
vennSL[["OtherSL"]]<-unique(rownames(DFRsplit)[DFRsplit[,"SL"]==1])
vennSL[["LoFmutationSL"]]<-unique(rownames(DFRsplit)[DFRsplit[,"LoFSLmut"]==1])
vennSL[["DriverGeneSL"]]<-unique(rownames(DFRsplit)[DFRsplit[,"LoFSLother"]==1])
vennSL[["Paralogs"]]<-unique(rownames(DFRsplit)[DFRsplit[,"Paralog"]==1])
venn.diagramCP(vennSL,paste0(outputdata,"VennSL.pdf"),imagetype="pdf",height=5,width=5,units="in")



DFRdataAll3<-get_DMA(DFRe,plist,LoFgenes,"FEATURE",cancerDrivers,splitMarkers = FALSE)
DFRdata3<-DFRdataAll3$DFAresults
save(DFRdata3,file=paste0(outputdata,"/DFRdata3.Rdata"))

DFRsplit3<-DFRdataAll3$TargetTypeMatrix
save(DFRsplit3,file=paste0(outputdata,"/DFRsplit3.Rdata"))
vennset3<-list()
vennset3[["Addiction"]]<-unique(rownames(DFRsplit3)[rowSums(DFRsplit3[,c("Target_Activating","Target_cancerDriverAct","SelfAddiction_CN","SelfAddiction_Expression","SelfAddiction_Mutation","SelfAddiction_Variant","SelfAddiction_Protein")])>0])
vennset3[["SyntheticLethal"]]<-unique(rownames(DFRsplit3)[rowSums(DFRsplit3[,c("LoFSLmut","LoFSLother","Paralog","SL")])>0])
vennset3[["Composite"]]<-unique(rownames(DFRsplit3)[rowSums(DFRsplit3[,c("MSI","Tsubtype")])>0])
venn.diagramCP(vennset3,paste0(outputdata,"Venn3Groups3.pdf"),imagetype="pdf",height=5,width=5,units="in")
vennAdd3<-list()
vennAdd3[["OtherAddiction"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"Target_Activating"]==1])
vennAdd3[["DriverGeneAddiction"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"Target_cancerDriverAct"]==1])
vennAdd3[["SelfAddiction"]]<-unique(rownames(DFRsplit3)[rowSums(DFRsplit3[,c("SelfAddiction_CN","SelfAddiction_Expression","SelfAddiction_Mutation","SelfAddiction_Variant","SelfAddiction_Protein")])>0])
venn.diagramCP(vennAdd3,paste0(outputdata,"VennAddiction3.pdf"),imagetype="pdf",height=5,width=5,units="in")

vennSL3<-list()
vennSL3[["OtherSL"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"SL"]==1])
vennSL3[["LoFmutationSL"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"LoFSLmut"]==1])
vennSL3[["DriverGeneSL"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"LoFSLother"]==1])
vennSL3[["Paralogs"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"Paralog"]==1])
venn.diagramCP(vennSL3,paste0(outputdata,"VennSL3.pdf"),imagetype="pdf",height=5,width=5,units="in")

