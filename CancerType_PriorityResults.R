

#################################################################################################
#                                                                                               #
#   Load Priority Results for Cancer Type analysis                                              #
#   Basic overview of number of targets in different tractability groups                        #                                                                                           
#                                                                                               #
#################################################################################################
load(file=paste0(inputdata,"TissueTypeColors.Rdata"))


GLOBAL<-read.table(file=paste0(inputdata,'Table_S9_allPriority_WithPPI_all.txt'),sep="\t",header=T,stringsAsFactors = F)

ScoreID<-paste(GLOBAL$ctype,GLOBAL$TARGET,sep="-")

GLOBAL$ID<-ScoreID


GLOBAL$numberBMs<-unlist(sapply(GLOBAL$MARKER,function(x) length(unique(strsplit(x,"//",fixed=T)[[1]]))))

cat(paste("Percent of targets, with only 1 marker:",sum(GLOBAL$numberBMs==1)/nrow(GLOBAL)),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=F)
cat(paste("Percent of targets, with 2 marker:",sum(GLOBAL$numberBMs==2)/nrow(GLOBAL)),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
#number of targets with a marker:
cat(paste("Percent targets with a marker:",sum(GLOBAL$MARKER!="N/A")/nrow(GLOBAL)),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Number of targets, with non-zero RWR score:",length(unique(GLOBAL[GLOBAL$RWRscore!=0,"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

Wm<-sum(GLOBAL$MARKER!="N/A")/nrow(GLOBAL)*100
Na<-sum(GLOBAL$RWRscore=="100")/nrow(GLOBAL)*100
Nb<-sum(GLOBAL$RWRscore=="75")/nrow(GLOBAL)*100
Nc<-sum(GLOBAL$RWRscore=="50")/nrow(GLOBAL)*100
Nd<-sum(GLOBAL$RWRscore=="25")/nrow(GLOBAL)*100

bardata<-data.frame(Set=c("With Marker","Network A","Network B","Network C","Network D"),Proportion=c(Wm,Na,Nb,Nc,Nd),stringsAsFactors = FALSE)
pdf(paste0(outputdata,"Figure3c.pdf"))
ggplot(bardata,aes(x=Set,y=Proportion,fill=Set))+geom_bar(stat="identity")+coord_flip()+theme_bw()
dev.off()
GLOBAL[which(GLOBAL$MARKER==""),"MARKER"]<-NA


cat(paste("Number Cancer types with priority target:",length(unique(GLOBAL$ctype))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Total number of unique targets, cancer type basis:",length(unique(GLOBAL$TARGET))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Number Cancer types with priority target, tractability group 2:",length(unique(GLOBAL[GLOBAL$GROUP==2,"ctype"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Total number of unique targets, cancer type basis, tract group 2:",length(unique(GLOBAL[GLOBAL$GROUP==2,"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Total number of unique targets, other disease indication:",length(unique(GLOBAL[GLOBAL$INDICATION=="other Disease","TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Unique targets, cancer type basis, other disease indication:",paste(unique(GLOBAL[GLOBAL$INDICATION=="other Disease","TARGET"]),collapse=",")),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

cat(paste("Number of targets group 1:",length(unique(GLOBAL[GLOBAL$GROUP==1,"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Number of targets group 2:",length(unique(GLOBAL[GLOBAL$GROUP==2,"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Number of targets group 3:",length(unique(GLOBAL[GLOBAL$GROUP==3,"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)


#################################################################################################
#                                                                                               #
#   Annotate results with DFA groups                                                            #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################


load(file=paste0(inputdata,"cancerDrivers.rdata"))
load(file=paste0(inputdata,"paralogList.Rdata"))
LoFgenes<-cancerDrivers[[1]][which(cancerDrivers[[2]]=="LoF")]


DFAdata<-get_DMA(GLOBAL,plist,LoFgenes,markerColumn="MARKER",cancerDrivers)
save(DFAdata,file=paste0(outputdata,"/DFAdata.Rdata"))
GdfSplit<-DFAdata$DFAresults
save(GdfSplit,file=paste0(outputdata,"/GdfSplit.Rdata"))

TTM<-DFAdata$TargetTypeMatrix

(colSums(TTM)/nrow(TTM))*100
TTall<-TTM
TTM<-TTM[,c(1:11,13:14)]
category<-colSums(TTall)
cat(paste("Proportion of targets, with an activating DFR all:",
          sum(category[c("SelfAddiction_Mutation","SelfAddiction_Expression",
          "SelfAddiction_CN","SelfAddiction_Variant","SelfAddiction_Protein","Target_cancerDriverAct","Target_Activating")])/nrow(TTall)),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Number of DMA, with an activating DFR all:",
          length(grep("SelfAddiction_Mutation|SelfAddiction_Expression|
                         SelfAddiction_CN|SelfAddiction_Variant|SelfAddiction_Protein|OncogenicAddiction_Act|AddictionND",GdfSplit$SLgroup))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

uT<-GdfSplit[grep("SelfAddiction_Mutation|SelfAddiction_Expression|
                         SelfAddiction_CN|SelfAddiction_Variant|SelfAddiction_Protein|OncogenicAddiction_Act|AddictionND",GdfSplit$SLgroup),"TARGET"]
cat(paste("Number of targets, with an activating DFR all:",
          length(unique(uT))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Number of DMAs, self addiction, unique targets in drivers:",length(unique(GdfSplit[grep("SelfAddiction_Mutation|SelfAddiction_Variant|SelfAddiction_CN|SelfAddiction_Protein|SelfAddiction_Expression",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

cat(paste("Number unique targets, addiction non-driver:",length(unique(GdfSplit[grep("AddictionND",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)


cat(paste("Proportion of targets, with SL all:",sum(category[c("Paralog","LoFSLmut","LoFSLother","SyntheticLethal")])/nrow(TTall)),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

cat(paste("Number of targets, SL :",length(unique(GdfSplit[grep("LoF_other_SL|LoF_mutation_SL|SyntheticLethal|Paralog",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

cat(paste("Number of targets, SL lof  drivers:",length(unique(GdfSplit[grep("LoF_other_SL|LoF_mutation_SL",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

cat(paste("Number of targets, Paralogs:",length(unique(GdfSplit[grep("Paralog",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)


cat(paste("Proportion of targets, with Composite all:",sum(category[c("MSI","Tsubtype")])/nrow(TTall)),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Number of targets, Composite:",length(unique(GdfSplit[grep("MSI|Transcriptional_Subtype",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)


uTDriver<-intersect(GdfSplit[grep("OncogenicAddiction_Act",GdfSplit$SLgroup),"TARGET"],cancerDrivers[[1]])


vennset<-list()
vennset[["Addiction"]]<-unique(uT)
vennset[["SyntheticLethal"]]<-unique(GdfSplit[grep("LoF_other_SL|LoF_mutation_SL|SyntheticLethal|Paralog",GdfSplit$SLgroup),"TARGET"])
vennset[["Composite"]]<-unique(GdfSplit[grep("MSI|Transcriptional_Subtype",GdfSplit$SLgroup),"TARGET"])
venn.diagramCP(vennset,paste0(outputdata,"Venn3GroupsPT.pdf"),imagetype="pdf",height=5,width=5,units="in")
vennAdd<-list()
vennAdd[["OtherAddiction"]]<-unique(GdfSplit[grep("AddictionND",GdfSplit$SLgroup),"TARGET"])
vennAdd[["DriverGeneAddiction"]]<-unique(uTDriver)
vennAdd[["SelfAddiction"]]<-unique(GdfSplit[grep("SelfAddiction_Mutation|SelfAddiction_Variant|SelfAddiction_CN|SelfAddiction_Protein|SelfAddiction_Expression",GdfSplit$SLgroup),"TARGET"])
venn.diagramCP(vennAdd,paste0(outputdata,"VennAddictionPT.pdf"),imagetype="pdf",height=5,width=5,units="in")

vennSL<-list()
vennSL[["OtherSL"]]<-unique(GdfSplit[grep("SyntheticLethal",GdfSplit$SLgroup),"TARGET"])
vennSL[["LoFmutationSL"]]<-unique(GdfSplit[grep("LoF_mutation_SL",GdfSplit$SLgroup),"TARGET"])
vennSL[["DriverGeneSL"]]<-unique(GdfSplit[grep("LoF_other_SL",GdfSplit$SLgroup),"TARGET"])
vennSL[["Paralogs"]]<-unique(unique(GdfSplit[grep("Paralog",GdfSplit$SLgroup),"TARGET"]))
venn.diagramCP(vennSL,paste0(outputdata,"VennSLPT.pdf"),imagetype="pdf",height=5,width=5,units="in")

#################################################################################################
#                                                                                               #
#   Priority plots for each cancer type individually and together                               #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################

load(paste0(inputdata,'/39_allMarkers_all.Rdata'))
load(paste0(inputdata,'/priority_threshold_all.Rdata'))
th<-priority_threshold
MARKERclass<-allMarkers
allctypes<-unique(allMarkers$ANALYSIS)

#Super priority Plot;
#GLOBAL<-read.table(sep='\t',file=paste0(inputdata,"Table_TotalRes.txt"),header=T,stringsAsFactors = F)
CTC<-TissueTypeColors
names(CTC)<-make.names(names(TissueTypeColors))
TOTRES<-GLOBAL[,c("ctype","TARGET","PRIORITY","TRACTABILITY","MARKERCLASS","RWRscore")]
colnames(TOTRES)[4]<-"BUCKET"

superPriorityPlot(TOTRES=TOTRES,allMarkers=allMarkers,plotname=outputdata,TissueColors=CTC,shape="Biomarker",indi="RWR",Priority="PRIORITY")

#################################################################################################
#                                                                                               #
#   Example score results with biomarkers                                                       #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################

load(paste0(inputdata,"07_EssMatrix_bDepletionsB2.Rdata"))
load(paste0(inputdata,"07_EssMatrix_qnorm_corrected_logFCs.RData"))


load(paste0(inputdata,"EXPpcBC.Rdata"))
ExpBEM<-EXPpcBC

ProtBEM<-readRDS(file=paste0(inputdata,"ProteomicMarkers.Rds"))
discBEM<-readRDS(file=paste0(inputdata,'ConsBEM.Rds'))

cmp<-read.csv(paste0(inputdata,"model_list_latest.csv"),header=T,stringsAsFactors = F)
cmp2<-cmp
cmp2$model_id<-cmp2$BROAD_ID
cmp<-rbind(cmp,cmp2)

load(paste0(inputdata,"pcgene.rdata"))

#SETDB1 with PIWIL4 in AML
plot_Biomarkers("SETDB1","PIWIL4","Exp",ExpBEM,ctype="Acute Myeloid Leukemia",logFC=qnorm_corrected_logFCs,annot=cmp,pointcol = TissueTypeColors["Acute.Myeloid.Leukemia"],outputdir = outputdata)

#DLL1 with APOBEC3B in Rhabdomyosarcoma
plot_Biomarkers("DLL1","APOBEC3B","Exp",ExpBEM,ctype="Rhabdomyosarcoma",logFC=qnorm_corrected_logFCs,annot=cmp,pointcol = TissueTypeColors["Rhabdomyosarcoma"],outputdir = outputdata)

#RAB18 with UDGH
plot_Biomarkers("RAB18","UGDH","Prot",ProtBEM,ctype="Kidney Carcinoma",logFC=qnorm_corrected_logFCs,annot=cmp,pointcol = TissueTypeColors["Kidney.Carcinoma"],outputdir =outputdata)

plot_Biomarkers("MITF","MITF","Exp",ExpBEM,ctype="Melanoma",logFC=qnorm_corrected_logFCs,annot=cmp,pointcol = TissueTypeColors["Melanoma"],outputdir = outputdata)
plot_Biomarkers("CAB39","REG4","Exp",ExpBEM,ctype="Gastric Carcinoma",logFC=qnorm_corrected_logFCs,annot=cmp,pointcol = TissueTypeColors["Gastric.Carcinoma"],outputdir = outputdata)
plot_Biomarkers("STX4","RAB25","Exp",ExpBEM,ctype="Ovarian Carcinoma",logFC=qnorm_corrected_logFCs,annot=cmp,pointcol = TissueTypeColors["Ovarian.Carcinoma"],outputdir = outputdata)
plot_Biomarkers("STXBP3","RAB25","Exp",ExpBEM,ctype="Ovarian Carcinoma",logFC=qnorm_corrected_logFCs,annot=cmp,pointcol = TissueTypeColors["Ovarian.Carcinoma"],outputdir = outputdata)
MetBEM<-readRDS("~/Library/CloudStorage/GoogleDrive-cp16@sanger.ac.uk/My Drive/PaperPriorityv2/GeneratePipelineInputs/MetaboliteData.Rds")

plot_Biomarkers("MUC20","kynurenine","Met",MetBEM,ctype="Non Small Cell Lung Carcinoma",logFC=qnorm_corrected_logFCs,annot=cmp,pointcol = TissueTypeColors["Non.Small.Cell.Lung.Carcinoma"],outputdir = outputdata)


#################################################################################################
#                                                                                               #
#   Compare to the result from Project Score                                                    #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################

ScoreMarkers<-read.csv(paste0(inputdata,"Score_Markers.csv"),header=T,stringsAsFactors = FALSE)


GlobalS2<-GLOBAL
load(paste0(inputdata,"15_allPriorityTargets.RData"))
allgenes<-unique(GLOBAL$TARGET)

load(paste0(inputdata,"allSymbol.Rdata"))

allMap<-makeNameMap(allgenes,allSymbol)
genenames<-updateNames(GLOBAL$TARGET,allMap,FALSE)
GLOBAL$TARGET<-genenames
TissueScore<-unique(ScoreMarkers$Cancer.Tissue.Type)
names(TissueScore)<-c("Bone","Breast","Central.Nervous.System","Large.Intestine","Esophagus",
                      "Haematopoietic.and.Lymphoid","Head.and.Neck","Peripheral.Nervous.System",
                      "Ovary","PanCancer","Stomach","Lung","Head.and.Neck","Pancreas","Lung")
TissueScore<-TissueScore[!TissueScore=="PanCancer"]
GlobalCL<-GLOBAL[GLOBAL$ctype!="PANCAN",]
CoreGenes<-list()
PrevPriority<-list()
KeepScore<-list()
load(paste0(inputdata,"10_PANCANCER_coreFitness_genes.Rdata"))
#Showing exact FDR for the same ones tested
for(i in 1:length(TissueScore)){
  load(paste0(inputdata,"09_ADM_",names(TissueScore)[i],"_coreFitnessGenes.Rdata"))
  CoreGenes[[TissueScore[i]]]<-coreFitnessGenes
  PrevPriority[[TissueScore[i]]]<-intersect(coreFitnessGenes,GlobalCL[GlobalCL$ctype==TissueScore[i],"TARGET"])
  KeepScore[[TissueScore[i]]]<-GlobalCL[GlobalCL$ctype==TissueScore[i]&!GlobalCL$TARGET%in%union(coreFitnessGenes,PanCancerCoreFitnessGenes),]
}

ScoreRemain<-do.call(rbind,KeepScore)
ScoreRemain$ctype<-make.names(ScoreRemain$ctype)
BothCtype<-intersect(GlobalS2$ctype,ScoreRemain$ctype)
ScoreRemain$usectype<-ScoreRemain$ctype
GlobalS2$usectype<-GlobalS2$ctype
GlobalS2[GlobalS2$ctype=="Ewing.s.Sarcoma","usectype"]<-"Bone"
GlobalS2[GlobalS2$ctype=="Glioma","usectype"]<-"Central.Nervous.System"
GlobalS2[GlobalS2$ctype=="Glioblastoma","usectype"]<-"Central.Nervous.System"
GlobalS2[GlobalS2$ctype=="Esophageal.Squamous.Cell.Carcinoma","usectype"]<-"Esophagus"
GlobalS2[GlobalS2$ctype=="Acute.Myeloid.Leukemia","usectype"]<-"Haematopoietic.and.Lymphoid"
GlobalS2[GlobalS2$ctype=="B.Cell.Non.Hodgkin.s.Lymphoma","usectype"]<-"Haematopoietic.and.Lymphoid"
GlobalS2[GlobalS2$ctype=="Plasma.Cell.Myeloma","usectype"]<-"Haematopoietic.and.Lymphoid"
GlobalS2[GlobalS2$ctype=="Non.Small.Cell.Lung.Carcinoma","usectype"]<-"Lung.Adenocarcinoma"
GlobalS2[GlobalS2$ctype=="Small.Cell.Lung.Carcinoma","usectype"]<-"Lung.Adenocarcinoma"

GlobalS2$origctype<-"Yes"
GlobalS2[GlobalS2$usectype%in%c("Bone","Central.Nervous.System","Esophagus",
                                "Haematopoietic.and.Lymphoid","Lung.Adenocarcinoma"),"origctype"]<-"No"

GlobalS2$ID<-paste(GlobalS2$TARGET,GlobalS2$usectype,sep="-")
ScoreRemain$ID<-paste(ScoreRemain$TARGET,ScoreRemain$ctype,sep="-")

AllTargets<-unique(c(GlobalS2$ID,ScoreRemain$ID))

load(paste0(inputdata,"PRIORITY_vectors_Score1.RData"))
PvecS1<-PRIORITY_vectors
names(PvecS1)<-make.names(names(PvecS1))
allS1<-NULL
for(i in 1:length(PvecS1)){
  Data<-PvecS1[[i]]
  ctype<-names(PvecS1)[i]
  allMap<-makeNameMap(names(Data),allSymbol)
  genenames<-updateNames(names(Data),allMap,FALSE)
  names(Data)<-genenames
  temp<-data.frame(ctype=ctype,Target=genenames,PRIORITY=Data,ID=paste(genenames,ctype,sep="-"),origctype="Yes",stringsAsFactors = FALSE)
  allS1<-rbind(allS1,temp)
}

load(paste0(inputdata,"PvectorsL3.Rdata"))
PvecS2<-PvectorsL3

names(PvecS2)[names(PvecS2)=="Ewing.s.Sarcoma"]<-"Bone"
names(PvecS2)[names(PvecS2)=="Glioma"]<-"Central.Nervous.System"
names(PvecS2)[names(PvecS2)=="Glioblastoma"]<-"Central.Nervous.System"
names(PvecS2)[names(PvecS2)=="Esophageal.Squamous.Cell.Carcinoma"]<-"Esophagus"
names(PvecS2)[names(PvecS2)=="Acute.Myeloid.Leukemia"]<-"Haematopoietic.and.Lymphoid"
names(PvecS2)[names(PvecS2)=="B.Cell.Non.Hodgkin.s.Lymphoma"]<-"Haematopoietic.and.Lymphoid"
names(PvecS2)[names(PvecS2)=="Plasma.Cell.Myeloma"]<-"Haematopoietic.and.Lymphoid"
names(PvecS2)[names(PvecS2)=="Non.Small.Cell.Lung.Carcinoma"]<-"Lung.Adenocarcinoma"
names(PvecS2)[names(PvecS2)=="Small.Cell.Lung.Carcinoma"]<-"Lung.Adenocarcinoma"

allS2<-NULL
for(i in unique(names(PvecS2))){
  sell<-which(names(PvecS2)==i)
  if(length(sell)==1){
  Data<-PvecS2[[i]]
  ctype<-i
  octype<-"Yes"
  if(ctype%in%c("Bone","Central.Nervous.System","Esophagus","Haematopoietic.and.Lymphoid","Lung.Adenocarcinoma")){
    octype<-"No"
  }
  temp<-data.frame(ctype=ctype,Target=names(Data),PRIORITY=Data,ID=paste(names(Data),ctype,sep="-"),origctype=octype,stringsAsFactors = FALSE)
  allS2<-rbind(allS2,temp)
  }else{
    Data<-NULL
    allg<-unique(unlist(lapply(PvecS2[sell],function(x) names(x))))
    for(k in sell){
      if(is.null(Data)){
        Data<-rep(NA,length(allg))
        names(Data)<-allg
        Data[names(PvecS2[[k]])]<-PvecS2[[k]]
        
      }else{
        temp<-rep(NA,length(allg))
        names(temp)<-allg
        temp[names(PvecS2[[k]])]<-PvecS2[[k]]

        t2<-cbind(Data,temp)
        Data<-apply(t2,1,function(x) max(x,na.rm=T))
        names(Data)<-allg
      }
    }
    ctype<-i
    octype<-"Yes"
    if(ctype%in%c("Bone","Central.Nervous.System","Esophagus","Haematopoietic.and.Lymphoid","Lung.Adenocarcinoma")){
      octype<-"No"
    }
    temp<-data.frame(ctype=ctype,Target=names(Data),PRIORITY=Data,ID=paste(names(Data),ctype,sep="-"),origctype=octype,stringsAsFactors = FALSE)
    allS2<-rbind(allS2,temp)
  }
}


load(paste0(inputdata,"/L3all.Rdata"))
uctype<-L3all[,"ct"]

uctype[uctype=="Ewing.s.Sarcoma"]<-"Bone"
uctype[uctype=="Glioma"]<-"Central.Nervous.System"
uctype[uctype=="Glioblastoma"]<-"Central.Nervous.System"
uctype[uctype=="Esophageal.Squamous.Cell.Carcinoma"]<-"Esophagus"
uctype[uctype=="Acute.Myeloid.Leukemia"]<-"Haematopoietic.and.Lymphoid"
uctype[uctype=="B.Cell.Non.Hodgkin.s.Lymphoma"]<-"Haematopoietic.and.Lymphoid"
uctype[uctype=="Plasma.Cell.Myeloma"]<-"Haematopoietic.and.Lymphoid"
uctype[uctype=="Non.Small.Cell.Lung.Carcinoma"]<-"Lung.Adenocarcinoma"
uctype[uctype=="Small.Cell.Lung.Carcinoma"]<-"Lung.Adenocarcinoma"



CompareRes<-data.frame(Target=AllTargets,Score=allS1[match(AllTargets,allS1$ID),"PRIORITY"],
                       Score2=allS2[match(AllTargets,allS2$ID),"PRIORITY"],
                       origctype=allS2[match(AllTargets,allS2$ID),"origctype"],L3class=L3all[match(AllTargets,L3all[,"ID"]),"BMclass"])
CompareRes[is.na(CompareRes$Score),"Score"]<-(-1)
CompareRes[is.na(CompareRes$Score2),"Score2"]<-(-1)
CompareRes[is.na(CompareRes$origctype),"origctype"]<-"Yes"
CompareRes$ctype<-unlist(sapply(CompareRes$Target,function(x) strsplit(x,"-",fixed=T)[[1]][2]))
CompareRes[CompareRes$ctype==1,"ctype"]<-"Lung.Adenocarcinoma"

pdf(paste0(outputdata,"/Compare_Priority_ScorevScore2.pdf"),useDingbats = FALSE,width=8,height=4)
print(ggplot(data=CompareRes,aes(x=Score2,y=Score,colour=ctype))+geom_point()+theme_bw())
dev.off()

Score2Only<-CompareRes[CompareRes$Score==(-1),]
ggplot(data=CompareRes,aes(x=ctype))+geom_bar()+coord_flip()
Set<-rep(1,nrow(CompareRes))
Set[CompareRes$Score>=41&CompareRes$Score2<36]<-2
Set[CompareRes$Score<41&CompareRes$Score2>=36]<-3

CompareRes<-cbind(CompareRes,Set)
#ggplot(data=MarkerData,aes(x=Source,fill=MARKER_TYPE))+geom_bar(position="fill")+ylab("Proportion of all markers")+scale_fill_manual(values=MarkerColours)
#Set 1 - priority in both
#Set 2 - priority Score only
#Set 3 - priority Score 2 only
CompareRes$Set<-factor(CompareRes$Set,levels=c(2,1,3))

coul <- brewer.pal(6, "Set3") 
names(coul)<-unique(CompareRes$L3class)
CompareRes$L3class<-factor(CompareRes$L3class,levels=c("0","A","B","C","D","NA"))
pdf(paste0(outputdata,"/Barchart_NetworkScores.pdf"),useDingbats = FALSE)
ggplot(data=CompareRes,aes(x=Set,fill=L3class))+geom_bar(position="fill")+theme_bw()+scale_fill_manual(values=coul)
dev.off()
#show best p-value for each marker on each one:
Score1anova<-read.delim(paste0(inputdata,"/all_CS_hits_Score1.txt"),header=T,stringsAsFactors = FALSE,row.names=NULL)
Score2ct<-unique(GlobalS2$ctype)

load(file=paste0(inputdata,"allHits_S2.Rdata"))


allHits_S1<-read.delim(paste0(inputdata,"/all_CS_hits_Score1.txt"),header=T,stringsAsFactors = FALSE,row.names=NULL)
logpval<-(-1*log(allHits_S1[,"FEATURE_ANOVA_pval"],10))
allHits_S1<-cbind(allHits_S1,logpval)
#show where cut-off is for Score and Score2
#pull markers from GlobalS2 and max class in Score_Markers

load(file=paste0(inputdata,"/GlobalS2.Rdata"))
load(file=paste0(inputdata,"/GlobalS1.Rdata"))

CompareRes$ScoreLpval<-as.numeric(GlobalS1[match(CompareRes$Target,GlobalS1$ID),"logpval"])
CompareRes$Score2Lpval<-as.numeric(GlobalS2[match(CompareRes$Target,GlobalS2$ID),"logpval"])

#show  number of marker types for Score2  and type i.e. large number of expression markers.
GlobalS2[GlobalS2$MARKERCLASS%in%c("D","N/A"),"MARKERCLASS"]<-""
MarkerData2<-GlobalS2[,c("MARKERCLASS","MARKER","MARKER_TYPE")]
MarkerData1<-GlobalS1[,c("Marker.Class","MARKER")]
MarkerData1$MARKER_TYPE<-""
m1<-grep("mut",MarkerData1$MARKER)
m2<-grep("_cna",MarkerData1$MARKER)
m3<-grep("MSI_status",MarkerData1$MARKER)
MarkerData1[m1,"MARKER_TYPE"]<-"mut"
MarkerData1[m2,"MARKER_TYPE"]<-"cna"
MarkerData1[m3,"MARKER_TYPE"]<-"Composite"
MarkerData1[intersect(m1,m2),"MARKER_TYPE"]<-c("mut,cna")
MarkerData1[intersect(m1,m3),"MARKER_TYPE"]<-c("mut,composite")
MarkerData1[intersect(m2,m3),"MARKER_TYPE"]<-c("cna,composite")
colnames(MarkerData1)[1]<-"MARKERCLASS"
MarkerData1$Source<-"Score"
MarkerData2$Source<-"Score2"

MarkerData<-rbind(MarkerData1,MarkerData2)
MarkerData[MarkerData$MARKER_TYPE=="N/A","MARKER_TYPE"]<-""
MarkerData[MarkerData$MARKER_TYPE=="","MARKER_TYPE"]<-"None"
MarkerColours<-rep("blue",length(unique(MarkerData$MARKER_TYPE)))


# Classic palette BuPu, with 4 colors
coul <- brewer.pal(10, "RdYlBu") 

# Add more colors to this palette :
coul <- colorRampPalette(coul)(25)
MarkerColours<-coul[1:21]
names(MarkerColours)<-unique(MarkerData$MARKER_TYPE)
MarkerColours["None"]<-"grey"
MarkerData$MARKER_TYPE<-factor(MarkerData$MARKER_TYPE,levels=c("None",sort(setdiff(names(MarkerColours),"None"))))
pdf(paste0(outputdata,"/BarChart_markerTypes_Compare.pdf"),useDingbats = FALSE)
print(ggplot(data=MarkerData,aes(x=Source,fill=MARKER_TYPE))+geom_bar(position="fill")+ylab("Proportion of all markers")+scale_fill_manual(values=MarkerColours))
dev.off()

CompareRes$NewCT<-1
CompareRes[CompareRes$ctype%in%c("Esophagus","Lung.Adenocarcinoma","Haematopoietic.and.Lymphoid","Central.Nervous.System"),"NewCT"]<-2
CompareRes[CompareRes$ctype%in%c("Rhabdomyosarcoma","Neuroblastoma","Melanoma",
                                 "Kidney.Carcinoma","Hepatocellular.Carcinoma",
                                 "Endometrial.Carcinoma","Cervical.Carcinoma",
                                 "Bladder.Carcinoma","Biliary.Tract.Carcinoma"),"NewCT"]<-3




UseT_S1<-GlobalS1[GlobalS1$ID%in%CompareRes$Target,]



load(file=paste0(inputdata,"/allplotFC.Rdata"))
load(file=paste0(inputdata,"/allplotFC1.Rdata"))

pdf(paste0(outputdata,"/FC_values_Score2.pdf"),useDingbats = FALSE,width=4,height=2.5)
print(ggplot(data=allplotFC,aes(x=IDs,y=value,colour=value.1))+geom_point(size=0.4)+scale_color_brewer(palette="YlOrRd")+theme_bw())
dev.off()
allplotFC1$value.1<-as.factor(allplotFC1$value.1)
pdf(paste0(outputdata,"/FC_values_Score1.pdf"),useDingbats = FALSE,width=4,height=2.5)
print(ggplot(data=allplotFC1,aes(x=IDs,y=value,colour=value.1))+geom_point(size=0.4)+scale_color_brewer(palette="YlOrRd")+theme_bw())
dev.off()

allplotFC$FS<-13.2
allplotFC[allplotFC$value.1==2,"FS"]<-22.8
allplotFC[allplotFC$value.1==3,"FS"]<-32.4
allplotFC[allplotFC$value.1==4,"FS"]<-42


allplotFC1$FS<-8.75
allplotFC1[allplotFC1$value.1==2,"FS"]<-17.5
allplotFC1[allplotFC1$value.1==3,"FS"]<-26.25
allplotFC1[allplotFC1$value.1==4,"FS"]<-35
allplotFC1[allplotFC1$value.1==5,"FS"]<-43.75

S2_FS<-allplotFC%>%group_by(.,IDs)%>%summarise(meanScore=mean(FS),numberCL=length(FS))
S2_FS$ctype<-unlist(sapply(S2_FS$IDs,function(x) strsplit(x,'-',fixed=T)[[1]][2]))
S1_FS<-allplotFC1%>%group_by(.,IDs)%>%summarise(meanScore=mean(FS),numberCL=length(FS))
S1_FS$ctype<-unlist(sapply(S1_FS$IDs,function(x) strsplit(x,'-',fixed=T)[[1]][2]))
S1_FS[S1_FS$ctype==1,"ctype"]<-"Lung.Adenocarcinoma"

testDiffFS<-wilcox.test(S2_FS$meanScore,S1_FS$meanScore)
print(testDiffFS$p.value)
S2_FS$source<-"S2"
S1_FS$source<-"S1"
FS_compare<-rbind(S2_FS,S1_FS)
pdf(paste0(outputdata,"/FS_values_Compare.pdf"),useDingbats = FALSE,width=3,height=4)
print(ggplot(data=FS_compare,aes(x=source,y=meanScore,fill=source))+geom_boxplot()+ylim(0,45)+ylab("Fitness score")+theme_bw())
dev.off()
pdf(paste0(outputdata,"/FS_values_Score2.pdf"),useDingbats = FALSE,width=8,height=4)
print(ggplot(data=S2_FS,aes(x=IDs,y=meanScore,fill=ctype))+geom_bar(stat="identity")+ylim(0,45)+ylab("Fitness score"))
dev.off()
pdf(paste0(outputdata,"/FS_nCL_Score2.pdf"),useDingbats = FALSE,width=8,height=4)
print(ggplot(data=S2_FS,aes(x=IDs,y=numberCL,fill=ctype))+geom_bar(stat="identity")+ylim(0,50)+ylab("Number scoring cell lines"))
dev.off()

pdf(paste0(outputdata,"/FS_values_Score1.pdf"),useDingbats = FALSE,width=8,height=4)
print(ggplot(data=S1_FS,aes(x=IDs,y=meanScore,fill=ctype))+geom_bar(stat="identity")+ylim(0,45)+ylab("Fitness score"))
dev.off()
pdf(paste0(outputdata,"/FS_nCL_Score1.pdf"),useDingbats = FALSE,width=8,height=4)
print(ggplot(data=S1_FS,aes(x=IDs,y=numberCL,fill=ctype))+geom_bar(stat="identity")+ylim(0,50)+ylab("Number scoring cell lines"))
dev.off()


S2_FS$NB<-0
S2_FS$NS<-0
load(paste0(inputdata,"manifestCL.Rdata"))
load(paste0(inputdata,"TissueColours.Rdata"))


manifest$Source<-"Broad"
manifest[grep("SIDM",manifest$INSTITUTE_ID),"Source"]<-"Sanger"

mt<-unique(manifest[,c("CANCER_TYPE","Source","INSTITUTE_ID")])
mt<-data.frame(Ctype=unlist(mt$CANCER_TYPE),Source=unlist(mt$Source),ID=unlist(mt$INSTITUTE_ID),stringsAsFactors = FALSE)
Sourcecounts<-mt%>%dplyr::count(Ctype,Source)

Bcounts<-Sourcecounts[Sourcecounts[,"Source"]=="Broad",]
Scounts<-Sourcecounts[Sourcecounts[,"Source"]=="Sanger",]

S2_FS$NB<-Bcounts[match(S2_FS$ctype,Bcounts$Ctype),"n"]
S2_FS$NS<-Scounts[match(S2_FS$ctype,Scounts$Ctype),"n"]

S2_new<-melt(S2_FS[,c("IDs","NB","NS")])
S2_new[is.na(S2_new$value),"value"]<-0
pdf(paste0(outputdata,"/FS_nCL_Source.pdf"),useDingbats = FALSE,width=8,height=4)

print(ggplot(data=S2_new,aes(x=IDs,y=value,fill=variable))+geom_bar(stat="identity",position="stack")+ylim(0,80)+ylab("Number scoring cell lines")+scale_fill_manual(values=c("pink","purple")))
dev.off()

