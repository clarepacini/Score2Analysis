

#################################################################################################
#                                                                                               #
#   Load Priority Results for Cancer Type analysis                                              #
#   Basic overview of number of targets in different tractability groups                        #                                                                                           
#                                                                                               #
#################################################################################################
load(file=paste0(inputdata,"TissueTypeColors.Rdata"))
#supplementary Table 12:
GLOBAL<-read.table(file=paste0(inputdata,'Table_S9_allPriority_WithPPI_all.txt'),sep="\t",header=T,stringsAsFactors = F)



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

LoFgenes<-cancerDrivers[[1]][which(cancerDrivers[[2]]=="LoF")]

load(file=paste0(inputdata,"paralogList.Rdata"))


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
cat(paste("Number of DMAs, self addiction, unique targets:",length(unique(GdfSplit[grep("SelfAddiction_Mutation|SelfAddiction_Variant|SelfAddiction_CN|SelfAddiction_Protein|SelfAddiction_Expression",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

cat(paste("Number unique targets, addiction non-driver:",length(unique(GdfSplit[grep("AddictionND",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)


cat(paste("Proportion of targets, with SL all:",sum(category[c("Paralog","LoFSLmut","LoFSLother","SL")])/nrow(TTall)),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

cat(paste("Number of targets, SL :",length(unique(GdfSplit[grep("LoF_other_SL|LoF_mutation_SL|SL|Paralog",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

cat(paste("Number of targets, SL lof  drivers:",length(unique(GdfSplit[grep("LoF_other_SL|LoF_mutation_SL",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)

cat(paste("Number of targets, Paralogs:",length(unique(GdfSplit[grep("Paralog",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)


cat(paste("Proportion of targets, with Composite all:",sum(category[c("MSI","Tsubtype")])/nrow(TTall)),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)
cat(paste("Number of targets, Composite:",length(unique(GdfSplit[grep("MSI|Transcriptional_Subtype",GdfSplit$SLgroup),"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)


ParalogData<-DFAdata[[1]][DFAdata[[1]]$SLgroup=="Paralog",]
cat(paste("Unique Number of targets with paralog marker:",length(unique(ParalogData[,"TARGET"]))),file=paste0(outputdata,"/PriorityResultsCT.txt"),sep="\n",append=T)



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
TOTRES<-read.table(sep='\t',file=paste0(inputdata,"Table_TotalRes.txt"),header=T,stringsAsFactors = F)
CTC<-TissueTypeColors
names(CTC)<-make.names(names(TissueTypeColors))
TOTRES<-GLOBAL[,c("ctype","TARGET","PRIORITY","TRACTABILITY","MARKERCLASS","RWRscore")]
colnames(TOTRES)[4]<-"BUCKET"
superPriorityPlot(TOTRES=TOTRES,allMarkers=allMarkers,plotname=outputdata,TissueColors=CTC,shape="Biomarker",indi="RWR")


