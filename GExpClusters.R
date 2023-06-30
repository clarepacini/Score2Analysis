
load(file=paste0(inputdata,"CellOut_metadata.rdata"))
#Total number of patient samples per group and number of groups found by Seurat
cmp<-read.csv(paste0(inputdata,"model_list_latest.csv"),header=T,stringsAsFactors = F)
cmp2<-cmp
cmp2$model_id<-cmp$BROAD_ID
cmp<-rbind(cmp,cmp2)

#read in the parameter selection and info:
SeuratParams<-read.csv(paste0(inputdata,"SeuratParams.csv"),header=T,stringsAsFactors = FALSE)


load(file=paste0(inputdata,"/Nclust.Rdata"))

ClusterCompare<-data.frame(Group=rep(rownames(Nclust),3),SNN=Nclust[,1],DBscan=Nclust[,2],Spectrum=Nclust[,3])
ClusterCompare<-reshape2::melt(ClusterCompare,id="Group")
MC<-apply(Nclust,1,median)
names(MC)<-rownames(Nclust)
orderC<-names(MC)[order(MC,decreasing=F)]
ClusterCompare[,"Group"]<-factor(ClusterCompare[,"Group"],levels=orderC)
pdf(paste0(outputdata,"ClusteringNumberClusters.pdf"),useDingbats = FALSE,width=6,height=4)
print(ggplot(data=ClusterCompare,aes(x=Group,y=value,fill=variable))+geom_bar(stat="identity",position="dodge")+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=3))+coord_flip())
dev.off()

SPdata<-SeuratParams[SeuratParams$Group%in%rownames(Nclust),]
SP<-SPdata[,"L1"]
names(SP)<-SPdata[,"Group"]
orderCT<-names(SP)[order(SP,decreasing=T)]
SPdata[,"Group"]<-factor(SPdata[,"Group"],levels=orderCT)
SPdata$NChbd<-Nclust[SPdata[,"Group"],2]
SPdata$NCsnn<-Nclust[SPdata[,"Group"],1]
SPdata$NCspectrum<-Nclust[SPdata[,"Group"],3]


pdf(paste0(outputdata,"SeuratClusteringResults.pdf"),useDingbats = FALSE,width=4,height=4)
print(ggplot(data=SPdata,aes(x=Group,y=L1))+geom_point(colour='red')+geom_point(data=SPdata,aes(x=Group,y=M1),colour='blue')+
  geom_point(data=SPdata,aes(x=Group,y=H1),colour="purple")+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)))
dev.off()

pdf(paste0(outputdata,"SeuratClusteringNumberClusters.pdf"),useDingbats = FALSE,width=4,height=3)
print(ggplot(data=SPdata,aes(x=Group,y=NCsnn))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)))
dev.off()


pdf(paste0(outputdata,"HBscanClusteringResults.pdf"),useDingbats = FALSE,width=4,height=4)
print(ggplot(data=SPdata,aes(x=Group,y=L3))+geom_point(colour='red')+geom_point(data=SPdata,aes(x=Group,y=M3),colour='blue')+
  geom_point(data=SPdata,aes(x=Group,y=H3),colour="purple")+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)))
dev.off()

pdf(paste0(outputdata,"HBscanClusteringNumberClusters.pdf"),useDingbats = FALSE,width=4,height=3)
print(ggplot(data=SPdata,aes(x=Group,y=NChbd))+geom_bar(stat="identity")+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6))+scale_y_continuous(breaks=c(0,2,4,6,8),labels=c("0","2","4","6","8")))
dev.off()



load(file=paste0(inputdata,"accuracy.Rdata"))
accuracyDF<-data.frame(accuracy=unlist(accuracy[,"Accuracy"]),Group=as.character(accuracy[,6]),stringsAsFactors = FALSE)
accuracyDF[,"Group"]<-factor(accuracyDF[,"Group"],levels=orderCT)

pdf(paste0(outputdata,"SVMaccuracyResults.pdf"),useDingbats = FALSE,width=4,height=6)
print(ggplot(data=accuracyDF,aes(x=Group,y=accuracy))+geom_point(colour='orange')+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)))
dev.off()



load(file=paste0(inputdata,"perfRFDF.Rdata"))
pdf(paste0(outputdata,"RFaccuracyResults.pdf"),useDingbats = FALSE,width=4,height=4)
print(ggplot(data=perfRFDF,aes(x=Group,y=accuracy))+geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=6)))
dev.off()



