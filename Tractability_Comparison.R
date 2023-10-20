#################################################################################################
#                                                                                               #
#   Comparison of tractability buckets between Score and Score 2                                #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################

GLOBAL<-read.table(file=paste0(inputdata,'Table_S9_allPriority_WithPPI_all.txt'),sep="\t",header=T,stringsAsFactors = F)

PANCAN<-read.table(file=paste0(inputdata,'Table_S9_allPriority_WithPPI_allPANCAN.txt'),sep="\t",header=T,stringsAsFactors = F)
#tractability_both_manualAmenden_Nov18_EG.RData
load("~/Library/CloudStorage/GoogleDrive-cp16@sanger.ac.uk/My Drive/Score2Submission/InputData/tractability_both_manualAmenden_Nov18_EG.RData")
tract_old<-tractability_both
load("~/Library/CloudStorage/GoogleDrive-cp16@sanger.ac.uk/My Drive/Score2Submission/InputData/tractability_both.Rdata")

allTargets<-unique(c(GLOBAL$TARGET,PANCAN$TARGET))

tractData<-data.frame(target=allTargets,
                      tractOld=tract_old[match(allTargets,rownames(tract_old)),"min_bucket"],
                      tractNew=tractability_both[match(allTargets,rownames(tractability_both)),"min_bucket"],stringsAsFactors = FALSE)
tractData$tractDiff<-tractData$tractOld-tractData$tractNew
tdata<-table(tractData$tractDiff)



tdata2<-data.frame(diff=names(tdata),count=as.numeric(tdata),stringsAsFactors=FALSE)
tdata2$sign<-"positive"
tdata2[tdata2$diff<0,"sign"]<-"negative"
tdata2[tdata2$diff<0,"count"]<-(-1*as.numeric(tdata2[tdata2$diff<0,"count"]))
#get example labels:
tdata2$Examples<-""
tdata2[tdata2$diff==-6,"Examples"]<-tractData[which(tractData$tractDiff==-6),"target"]
tdata2[tdata2$diff==-2,"Examples"]<-paste(tractData[which(tractData$tractDiff==-2),"target"],collapse=", ")
tdata2[tdata2$diff==-1,"Examples"]<-paste(tractData[which(tractData$tractDiff==-1),"target"][c(5,8,9)],collapse=", ")
tdata2[tdata2$diff==1,"Examples"]<-paste(tractData[which(tractData$tractDiff==1),"target"][c(1,2,3)],collapse=", ")
tdata2[tdata2$diff==2,"Examples"]<-paste(tractData[which(tractData$tractDiff==2),"target"][c(3)],collapse=", ")
tdata2[tdata2$diff==3,"Examples"]<-paste(tractData[which(tractData$tractDiff==3),"target"][c(1,6,7)],collapse=", ")
tdata2[tdata2$diff==4,"Examples"]<-paste(tractData[which(tractData$tractDiff==4),"target"][c(1,2,3)],collapse=", ")
tdata2[tdata2$diff==5,"Examples"]<-paste(tractData[which(tractData$tractDiff==5),"target"][c(1,2)],collapse=", ")
tdata2[tdata2$diff==6,"Examples"]<-paste(tractData[which(tractData$tractDiff==6),"target"][c(1)],collapse=", ")

pdf(paste0(outputdata,"tractability_comparison.pdf"))
ggplot(tdata2, aes(x=diff, y=count, fill=sign)) + 
  geom_bar(stat="identity", position="identity")+coord_flip()+ylab("Number of genes")+
  xlab("Increase tractability in Score2")+scale_x_discrete(limits=as.character(seq(-6,6,1)))+
geom_text(aes(label=Examples), position=position_dodge(width=0.5), vjust=0.5,hjust=-0.05,size=3)+
  theme_bw()
dev.off()
#+geom_vline(xintercept=0,colour="black")
