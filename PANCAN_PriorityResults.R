

load(file=paste0(inputdata,"TissueTypeColors.Rdata"))

GLOBAL<-read.table(file=paste0(inputdata,'Table_S9_allPriority_WithPPI_allPANCAN.txt'),sep="\t",header=T,stringsAsFactors = F)



load(paste0(inputdata,'priority_threshold_all_PANCAN.RData'))
load(paste0(inputdata,'39_allMarkers_all_PANCAN.Rdata'))
load(paste0(inputdata,'39_allMarkers_all_PANCAN.Rdata'))
load(paste0(inputdata,'/Pvectors_PANCAN.RData'))
load(paste0(inputdata,'/PAI_PANCAN.RData'))
load(paste0(inputdata,'/PAB_PANCAN.RData'))
#################################################################################################
#                                                                                               #
#   Annotate results with DFA groups                                                            #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################
load(paste0(inputdata,"cancerDrivers.Rdata"))
LoFgenes<-cancerDrivers[[1]][which(cancerDrivers[[2]]=="LoF")]

load(file=paste0(inputdata,"paralogList.Rdata"))


DFAdata<-get_DFAgroup(GLOBAL,plist,LoFgenes,markerColumn="MARKER",cancerDrivers)
save(DFAdata,file=paste0(outputdata,"/DFAdataPANCAN.Rdata"))

ParalogData<-DFAdata[[1]][DFAdata[[1]]$SLgroup=="Paralog",]
cat(paste("Unique Number of targets with paralog marker:",length(unique(ParalogData[,"TARGET"]))),file=paste0(outputdata,"/PriorityResultsPANCAN.txt"),sep="\n",append=T)

#################################################################################################
#                                                                                               #
#   PANCAN priority plots                                                                       #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################

cvec<-c("PANCAN"="#85933B")
th<-priority_threshold
MARKERclass<-allMarkers


TOTRES<-GLOBAL[,c("ctype","TARGET","PRIORITY","TRACTABILITY","MARKERCLASS","RWRscore")]
colnames(TOTRES)[4]<-"BUCKET"

superPriorityPlot(TOTRES=TOTRES,allMarkers=allMarkers,plotname=outputdata,TissueColors=cvec,shape="Biomarker",indi="RWR",plotsuffix="SuperPriorityPlotPC.pdf",Priority="PRIORITY")


