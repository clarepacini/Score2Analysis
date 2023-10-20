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
#   Generate two groups of DFAs with different levels of stringency                           #
#                                                                                               #
#                                                                                               #
#################################################################################################


DFRa<-DFRP[DFRP$FDR<1&DFRP$f2>0.1&abs(DFRP$logFC)>1,]

ScoreM<-unique(DFRP$`Depleted Gene`)
colnames(DFRa)[colnames(DFRa)=="Depleted Gene"]<-"TARGET"
colnames(DFRa)[colnames(DFRa)=="depAssoc"]<-"ASSOCIATION_EFFECT"
save(DFRa,file=paste0(outputdata,"/DFRa.Rdata"))
write.table(DFRa,file=paste0(outputdata,"/SupplementaryTable7.tsv"),row.names=F,quote=F,sep="\t")

cat(paste("Number genes with strongest marker: ",length(unique(DFRa$TARGET))),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)

DFRa$M<-unlist(sapply(DFRa$FEATURE,function(x) strsplit(x,"_",fixed=T)[[1]][1]))
DFRa$Mtype<-unlist(sapply(DFRa$FEATURE,function(x) strsplit(x,"_",fixed=T)[[1]][length(strsplit(x,"_",fixed=T)[[1]])]))
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

deps<-unique(DFRe$TARGET)
temp<-NULL
for(i in deps){
  sel<-DFRe[which(DFRe$TARGET==i),]
  if(nrow(sel)==1){
    temp<-rbind(temp,sel)
  }else{
    selfeat<-which.min(sel$FDR)
    temp<-rbind(temp,sel[selfeat,])
  }
}
DFReb<-temp
#################################################################################################
#                                                                                               #
#   Power analysis of DMAs                                                                      #
#                                                                                               #
#                                                                                               #
#################################################################################################


DFRa$MarkerGroup<-"D"
DFRa[abs(DFRa$logFC)>1.5,"MarkerGroup"]<-"C"
DFRa[abs(DFRa$logFC)>2,"MarkerGroup"]<-"B"
DFRa[abs(DFRa$logFC)>2.5,"MarkerGroup"]<-"A"
Rsquarevals<-DFRa%>%group_by(MarkerGroup,Mtype)%>%summarise(mean(Rsquare))
Pvals<-DFRa%>%group_by(MarkerGroup,Mtype)%>%summarise(mean(P.Value))
#this says if the variables fit to a particular r2 gives the sample size required to achieve a level of power
#use the R2 from the regression output
powerlevels<-seq(0.01,0.999,by=0.001)
classaCN<-data.frame(pwr=powerlevels,omic="CN",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[1,3]), k = 1, power = x, alpha = as.numeric(Pvals[1,3]))$n))
classaComp<-data.frame(pwr=powerlevels,omic="Comp",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[2,3]), k = 1, power = x, alpha = as.numeric(Pvals[1,3]))$n))
classaExpr<-data.frame(pwr=powerlevels,omic="Expr",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[3,3]), k = 1, power = x, alpha = as.numeric(Pvals[1,3]))$n))
classaProt<-data.frame(pwr=powerlevels,omic="Prot",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[4,3]), k = 1, power = x, alpha = as.numeric(Pvals[1,3]))$n))
classbMet<-data.frame(pwr=powerlevels,omic="Met",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[8,3]), k = 1, power = x, alpha = as.numeric(Pvals[1,3]))$n))
classbMut<-data.frame(pwr=powerlevels,omic="Mut",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[9,3]), k = 1, power = x, alpha = as.numeric(Pvals[1,3]))$n))
classbvar<-data.frame(pwr=powerlevels,omic="Var",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[11,3]), k = 1, power = x, alpha = as.numeric(Pvals[1,3]))$n))

classaCN2<-data.frame(pwr=powerlevels,omic="CN2",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[18,3]), k = 1, power = x, alpha = as.numeric(Pvals[18,3]))$n))
classaComp2<-data.frame(pwr=powerlevels,omic="Comp2",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[19,3]), k = 1, power = x, alpha = as.numeric(Pvals[19,3]))$n))
classaExpr2<-data.frame(pwr=powerlevels,omic="Expr2",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[20,3]), k = 1, power = x, alpha = as.numeric(Pvals[20,3]))$n))
classaProt2<-data.frame(pwr=powerlevels,omic="Prot2",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[23,3]), k = 1, power = x, alpha = as.numeric(Pvals[23,3]))$n))
classbMet2<-data.frame(pwr=powerlevels,omic="Met2",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[21,3]), k = 1, power = x, alpha = as.numeric(Pvals[21,3]))$n))
classbMut2<-data.frame(pwr=powerlevels,omic="Mut2",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[22,3]), k = 1, power = x, alpha = as.numeric(Pvals[22,3]))$n))
classbvar2<-data.frame(pwr=powerlevels,omic="Var2",n=sapply(powerlevels,function(x) pwrss.f.reg(r2 = as.numeric(Rsquarevals[24,3]), k = 1, power = x, alpha = as.numeric(Pvals[24,3]))$n))

#add in the lowest class ones as well because increase the sample size....

Pwrdata<-rbind(classaCN,classaComp,classaExpr,classaProt,classbMet,classbMut,classbvar,
               classaCN2,classaComp2,classaExpr2,classaProt2,classbMet2,classbMut2,classbvar2)
par(mfrow=c(1,1))
pdf(paste0(outputdata,"/Power_DMAs.pdf"),useDingbats = FALSE,height=4,width=6)
ggplot(data=Pwrdata,aes(x=n,y=pwr,col=omic))+geom_line()+theme_bw()+ylim(0,1)+xlim(0,930)+
  geom_vline(xintercept = c(538,579,887,909,917,930),alpha=0.5)
dev.off()

#################################################################################################
#                                                                                               #
#   Compare DMAs to Pan-cancer NMF signatures and EMT scores                                    #
#                                                                                               #
#                                                                                               #
#################################################################################################
load(paste0(inputdata,"/07_EssMatrix_qnorm_corrected_logFCs.RData"))

load(paste0(inputdata,"H.Rdata"))


Dep137<-qnorm_corrected_logFCs[unique(DFRa$TARGET),]
geneAnnotation<-read.delim(paste0(inputdata,'/protein-coding_gene.txt'),sep = '\t',header=TRUE,stringsAsFactors = FALSE)
rownames(geneAnnotation)<-geneAnnotation[,"symbol"]
TOTRES<-AllRegression(paste0(outputdata,'LM_results'),Dep137,H,inparallel=FALSE,
                      geneAnnotation,biomarkertype = "CExpr",Tissue=NULL,MSI=NULL,markerIsFactor=FALSE,scalingFeat=FALSE)

sum(TOTRES$FDR<1)
table(TOTRES[TOTRES$FDR<1,"FEATURE"])
length(unique(TOTRES[TOTRES$FDR<1,"Depleted Gene"]))

pdf(paste0(outputdata,"/SuppFERMT2_EMT.pdf"))
print(SummaryPlot("FERMT2","EMT","Exp",H,ctype="PANCAN",logFC=qnorm_corrected_logFCs,annot=cmp,cellcol = NULL,SigMat=H,NULL))
dev.off()
pdf(paste0(outputdata,"/SuppATP1B3_Metabolic.pdf"))
print(SummaryPlot("ATP1B3","EMT","Exp",H,ctype="PANCAN",logFC=qnorm_corrected_logFCs,annot=cmp,cellcol = NULL,SigMat=H,NULL))
dev.off()
pdf(paste0(outputdata,"/SuppNAMPT_EMT.pdf"))
print(SummaryPlot("NAMPT","EMT","Exp",H,ctype="PANCAN",logFC=qnorm_corrected_logFCs,annot=cmp,cellcol = NULL,SigMat=H,NULL))
dev.off()

pdf(paste0(outputdata,"/SuppWWTR1_Metabolic.pdf"))
print(SummaryPlot("WWTR1","EMT","Exp",H,ctype="PANCAN",logFC=qnorm_corrected_logFCs,annot=cmp,cellcol = NULL,SigMat=H,NULL))
dev.off()
pdf(paste0(outputdata,"/SuppSTXBP3_Apoptosis.pdf"))
print(SummaryPlot("STXBP3","EMT","Exp",H,ctype="PANCAN",logFC=qnorm_corrected_logFCs,annot=cmp,cellcol = NULL,SigMat=H,NULL))
dev.off()
pdf(paste0(outputdata,"/SuppJUN_Apoptosis.pdf"))
print(SummaryPlot("JUN","EMT","Exp",H,ctype="PANCAN",logFC=qnorm_corrected_logFCs,annot=cmp,cellcol = NULL,SigMat=H,NULL))
dev.off()


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
#   Plot of dependencies and markers                                                            #
#                                                                                               #
#                                                                                               #
#################################################################################################
DFRP<-cbind("PANCAN",DFRP)
DFRP[,1]<-unlist(sapply(DFRP$FEATURE,function(x) strsplit(x,"_",fixed=TRUE)[[1]][length(strsplit(x,"_",fixed=TRUE)[[1]])]))
PlotColours<-c("lightblue","orange","darkgreen","yellow","purple","pink","darkblue")
names(PlotColours)<-unique(DFRP[,1])
getpoints<-globalVolcano2(DFRP,TissueColors=PlotColours,hitids=which(DFRP$FDR<1&DFRP$f2>0.1&abs(DFRP$logFC)>1),ctype="All",pvalcol = "P.Value",signCol =NA,effectSize = "logFC",legend=FALSE,FDRthresh=1,DEPcol="Depleted Gene")
pdf(paste0(outputdata,"Volcano_PanCanDMA.pdf"),height=6,width=5,useDingbats = FALSE)
print(globalVolcano2(DFRP,TissueColors=PlotColours,hitids=which(DFRP$FDR<1&DFRP$f2>0.1&abs(DFRP$logFC)>1),ctype="All",idpoints=getpoints$points,pvalcol = "P.Value",signCol = NA,effectSize = "logFC",legend=FALSE,FDRthresh=1,DEPcol="Depleted Gene"))
dev.off()

png(paste0(outputdata,"Volcano_PanCanDMA.png"),height=6,width=5,units="in", res=300, pointsize=12)
print(globalVolcano2(DFRP,TissueColors=PlotColours,hitids=which(DFRP$FDR<1&DFRP$f2>0.1&abs(DFRP$logFC)>1),ctype="All",idpoints=getpoints$points,pvalcol = "P.Value",signCol = NA,effectSize = "logFC",legend=FALSE,FDRthresh=1,DEPcol="Depleted Gene"))

dev.off()
#################################################################################################
#                                                                                               #
#   Compare to previous results                                                                 #
#                                                                                               #
#                                                                                               #
#################################################################################################

DataSetCol<-c("darkgrey","orange","lightblue","darkblue")
names(DataSetCol)<-c("RNAi","Score","Score2","Score2_2")

RNAi<-read.csv("~/Library/CloudStorage/GoogleDrive-cp16@sanger.ac.uk/My Drive/Score2Submission/InputData/RNAi_SixSigma.csv",header=T,stringsAsFactors = FALSE)

ScorePancanRes<-read.csv("~/Library/CloudStorage/GoogleDrive-cp16@sanger.ac.uk/My Drive/Score2Submission/InputData/PanCancerPriority_Score.csv",header=T,stringsAsFactors = FALSE)
ScorePancanRes<-ScorePancanRes[1:92,1:30]

load(file=paste0(outputdata,"/normLRTresnq.Rdata"))
LRT<-normLRTres$LRT

#Scatter plot showing LRT values for genes in different sets.
SixSigma<-unique(RNAi[RNAi$is.six.sigma,"Gene.dependency"])
ScoreG<-unique(ScorePancanRes$Gene)
Score2<-unique(DFRP$`Depleted Gene`)

#need to update SixSigma and ScoreG to update to new gene symbol 
allgenes<-unlist(c(SixSigma,ScoreG))
allgenes<-unique(allgenes)
load(paste0(inputdata,"allSymbol.Rdata"))

allMap<-makeNameMap(allgenes,allSymbol)
genenames<-updateNames(SixSigma,allMap,FALSE)
SixSigma<-genenames
genenames<-updateNames(ScoreG,allMap,FALSE)
ScoreG<-genenames
pdf(paste0(outputdata,"/LRT_allPapers.pdf"),useDingbats = FALSE,width=6,height=3)

plot(density(LRT[intersect(SixSigma,names(LRT))]),col=DataSetCol["RNAi"],main="normLRT of targets in different datasets",lty=1,lwd=2)
#lines(density(LRT[intersect(ScoreG,names(LRT))]),col=DataSetCol["Score"])
lines(density(LRT[Score2]),col=DataSetCol["Score2"],lty=1,lwd=2)
dev.off()

load(paste0(inputdata,'10_PANCANCER_coreFitness_genes.RData'))
load(paste0(inputdata,'10_PANCANCER_CommonEssentialsAUC.RData'))

#Plot showing numbers tested overlapping

vennsetR<-list()
vennsetR[["CoreEssentials"]]<-PanCancerCoreFitnessGenes
vennsetR[["Score2_500"]]<-Score2
vennsetR[["RNAi"]]<-SixSigma
venn.diagramCP(vennsetR,paste0(outputdata,"Venn3Groups_RNAi.pdf"),imagetype="pdf",height=5,width=5,units="in")
vennsetS<-list()
vennsetS[["CoreEssentials"]]<-PanCancerCoreFitnessGenes
vennsetS[["Score2_500"]]<-Score2
vennsetS[["Score"]]<-ScoreG
venn.diagramCP(vennsetS,paste0(outputdata,"Venn3Groups_Score.pdf"),imagetype="pdf",height=5,width=5,units="in")

FCdataS1<-read.delim(paste0(inputdata,"/01a_qnorm_corrected_logFCs.tsv"),sep="\t",header=T,row.names=1)

load(file=paste0(inputdata,"/normLRTresnqS1.Rdata"))

LRTS1<-normLRTresS1$LRT
Top500s2<-read.table(file=paste0(outputdata,"/Top500normLRT.tsv"),sep='\t',header=T)
s2Thresh<-min(LRT[Top500s2[,1]])
sum(LRTS1>s2Thresh)
pdf(paste0(outputdata,"/LRT_Score_Score2.pdf"),useDingbats = FALSE,width=6,height=3)

plot(density(LRTS1),col=DataSetCol["Score"],main="normLRT of genes in different datasets",lwd=2,xlim=c(0,1800))
#lines(density(LRT[intersect(ScoreG,names(LRT))]),col=DataSetCol["Score"])
lines(density(LRT),col=DataSetCol["Score2"],lwd=2)
abline(v=s2Thresh)
dev.off()

#################################################################################################
#                                                                                               #
#   Comparison of DMA oncogene associations in different studies                                #
#                                                                                               #
#                                                                                               #
#################################################################################################

RNAimarkers<-read.csv(paste0(inputdata,'RNAi_Oncogene.csv'),header=T,stringsAsFactors = FALSE)
setdiff(RNAimarkers$Gene.dependency.marker,cancerDrivers[[1]])

FoundRNAi<-unique(RNAimarkers[RNAimarkers$Rank<100,"Gene.dependency.marker"])
FoundRNAi<-FoundRNAi[!is.na(FoundRNAi)]
ResRNAi<-RNAimarkers[RNAimarkers$Rank<100,]

RNAigroups<-read.csv(paste0(inputdata,'RNAi_MDP.csv'),header=T,stringsAsFactors = FALSE)
RNAi6s<-RNAigroups[RNAigroups$is.six.sigma=="TRUE",]
RNAirel<-RNAi6s[RNAi6s$MDP.feature.group%in%c("Mutation driven","Expression driven","CYCLOPS"),]
RNAirel<-RNAirel[RNAirel$FDR<0.05,]
RNAionco<-intersect(RNAirel$Gene.dependency,cancerDrivers[[1]])
DFRP$M<-unlist(sapply(DFRP$FEATURE,function(x) strsplit(x,"_",fixed=T)[[1]][1]))
DFRP$PPI_min<-unlist(sapply(1:nrow(DFRP),function(x) ifelse(DFRP[x,"M"]==DFRP[x,"Depleted Gene"],0,6)))
DFRP<-DFRP[order(DFRP$FDR),]
splitByDep<-split.data.frame(DFRP,f=DFRP$`Depleted Gene`)

rankOnco<-lapply(splitByDep,function(x) x[which(x[,"PPI_min"]==0),"FDR"])
rank100<-unlist(lapply(rankOnco,function(x) ifelse(length(x)>0,x<1,FALSE)))
intersect(names(rank100[rank100]),RNAimarkers$Gene.dependency.marker)
Score2onco<-intersect(names(rank100[rank100]),cancerDrivers[[1]])

ScoreMarkers<-read.csv(paste0(inputdata,"Score_MarkerClass.csv"),header=T,stringsAsFactors = FALSE)
ScoreMarkers<-ScoreMarkers[ScoreMarkers$Target%in%cancerDrivers[[1]],]
CheckOnco<-apply(ScoreMarkers,1,function(x) length(grep(x["Target"],x["Cancer.Functional.Event"]))>0)
ScoreOnco<-unique(ScoreMarkers[CheckOnco,"Target"])

Allonco<-unique(c(ScoreOnco,Score2onco,RNAionco))

vennsetO<-list()
vennsetO[["RNAi"]]<-RNAionco
vennsetO[["Score2"]]<-Score2onco
vennsetO[["Score"]]<-ScoreOnco
venn.diagramCP(vennsetO,paste0(outputdata,"Venn3Groups_Onco.pdf"),imagetype="pdf",height=5,width=5,units="in")
RO<-setdiff(RNAionco,c(Score2onco,ScoreOnco))
setdiff(RNAionco,Score2onco)
intersect(RNAionco,ScoreOnco)
intersect(RNAionco,intersect(Score2onco,ScoreOnco))
Sonly<-setdiff(ScoreOnco,c(Score2onco,RNAionco))
setdiff(intersect(Score2onco,ScoreOnco),RNAionco)
setdiff(intersect(Score2onco,RNAionco),ScoreOnco)
S2O<-setdiff(Score2onco,union(ScoreOnco,RNAionco))

ROnly<-RNAirel[RNAirel$Gene.dependency%in%RO,]
RNAirel[RNAirel$Gene.dependency=="MET",]
#get marker types for each set- In RNAi and Score2 only
RS2<-setdiff(intersect(Score2onco,RNAionco),ScoreOnco)
ROnco<-RNAirel[RNAirel$Gene.dependency%in%RS2,]
allOnco<-lapply(splitByDep,function(x) x[which(x[,"PPI_min"]==0),])
oncoS<-allOnco[RS2]

oncoS<-lapply(oncoS,function(x) x[which(x[,"FDR"]<1),])

oncoS<-do.call('rbind',oncoS)
plotOnco<-oncoS[,c("Depleted Gene","FEATURE","FDR")]
colnames(plotOnco)<-c("Gene","Feature","FDR")
plotOnco$Feature<-unlist(sapply(plotOnco$Feature,function(x) strsplit(x,"_",fixed=T)[[1]][2]))
plotOnco$rank<-rank(plotOnco$FDR)
plotOnco2<-ROnco[,c("Gene.dependency","MDP.feature.group","FDR")]
colnames(plotOnco2)<-c("Gene","Feature","FDR")
plotOnco2[plotOnco2$Feature=="Expression driven","Feature"]<-"Expr"
plotOnco2[plotOnco2$Feature=="Mutation driven","Feature"]<-"mut"
plotOnco2[plotOnco2$Feature=="CYCLOPS","Feature"]<-"CN"
plotOnco2$rank<-rank(plotOnco2$FDR)
plotOnco<-cbind(plotOnco,source="Score2")
plotOnco2<-cbind(plotOnco2,source="RNAi")
pOnco<-rbind(plotOnco,plotOnco2)

#try using corrplot package:
pOnco$row<-paste(pOnco$source,pOnco$Feature,sep="-")
plotMat<-matrix(NA,nrow=6,ncol=15)
rownames(plotMat)<-c("Score2-Expr","RNAi-Expr","Score2-CN","RNAi-CN","Score2-mut","RNAi-mut")
colnames(plotMat)<-RS2
for(i in 1:nrow(pOnco)){
  plotMat[pOnco[i,"row"],pOnco[i,"Gene"]]<-(-1*pOnco[i,"rank"])
}
#plotMat<-(-log(plotMat,10))
plotMat[is.na(plotMat)]<-0
plotMat[plotMat!=0]<-1
pdf(paste0(outputdata,"/CompareMarkers_OncoScore2RNAi.pdf"),useDingbats = FALSE)
corrplot(plotMat,is.cor = FALSE,p.mat=plotMat,sig.level=2,insig="blank",tl.col="black")
dev.off()

#score2 only:
oncoSO<-allOnco[S2O]
oncoSO<-lapply(oncoSO,function(x) x[which(x[,"FDR"]<1),])
oncoSO<-do.call('rbind',oncoSO)

plotMat2<-matrix(NA,3,47)
rownames(plotMat2)<-c("Expr","CN","mut")
colnames(plotMat2)<-S2O

plotOncoS2<-oncoSO[,c("Depleted Gene","FEATURE","FDR")]
colnames(plotOncoS2)<-c("Gene","Feature","FDR")
plotOncoS2$Feature<-unlist(sapply(plotOncoS2$Feature,function(x) strsplit(x,"_",fixed=T)[[1]][2]))

#marker types all three
All<-intersect(intersect(Score2onco,RNAionco),ScoreOnco)
ROncoA<-RNAirel[RNAirel$Gene.dependency%in%All,]

oncoSA<-allOnco[All]
oncoSA<-lapply(oncoSA,function(x) x[which(x[,"FDR"]<1),])

oncoSA<-do.call('rbind',oncoSA)
plotOncoA<-oncoSA[,c("Depleted Gene","FEATURE","FDR")]
colnames(plotOncoA)<-c("Gene","Feature","FDR")
plotOncoA$Feature<-unlist(sapply(plotOncoA$Feature,function(x) strsplit(x,"_",fixed=T)[[1]][length(strsplit(x,"_",fixed=T)[[1]])]))

plotOnco2A<-ROncoA[,c("Gene.dependency","MDP.feature.group","FDR")]
colnames(plotOnco2A)<-c("Gene","Feature","FDR")
plotOnco2A[plotOnco2A$Feature=="Expression driven","Feature"]<-"Expr"
plotOnco2A[plotOnco2A$Feature=="Mutation driven","Feature"]<-"mut"
plotOnco2A[plotOnco2A$Feature=="CYCLOPS","Feature"]<-"CN"

SM<-ScoreMarkers[CheckOnco,]
plotOnco3A<-SM[SM$Target%in%All,]
plotOnco3A$Feature<-c("mut","mut","mut","mut","CN","CN")
plotOnco3A<-plotOnco3A[,c("Target","Feature","CFE.ANOVA.FDR")]
colnames(plotOnco3A)<-c("Gene","Feature","FDR")
plotOncoA<-cbind(plotOncoA,source="Score2")
plotOnco2A<-cbind(plotOnco2A,source="RNAi")
plotOnco3A<-cbind(plotOnco3A,source="Score")
pOncoA<-rbind(plotOncoA,plotOnco2A,plotOnco3A)

SM[SM$Target=="MET",]
SM[SM$Target%in%Sonly,]
#try using corrplot package:
pOncoA$row<-paste(pOncoA$source,pOncoA$Feature,sep="-")
plotMatA<-matrix(NA,nrow=10,ncol=6)
rownames(plotMatA)<-c("Score2-Expr","RNAi-Expr","Score-Expr","Score2-CN","RNAi-CN","Score-CN","Score2-mut","RNAi-mut","Score-mut"
                      ,"Score2-var")
colnames(plotMatA)<-All
for(i in 1:nrow(pOncoA)){
  plotMatA[pOncoA[i,"row"],pOncoA[i,"Gene"]]<-(pOncoA[i,"FDR"])
}
#plotMat<-(-log(plotMat,10))
plotMatA[is.na(plotMatA)]<-0
plotMatA[plotMatA!=0]<-1
pdf(paste0(outputdata,"/CompareMarkers_OncoAll.pdf"),useDingbats = FALSE)
corrplot(plotMatA,is.cor = FALSE,p.mat=plotMatA,sig.level=2,insig="blank",tl.col="black")
dev.off()

#marker types Score, Score2
SS2<-setdiff(intersect(Score2onco,ScoreOnco),RNAionco)

oncoSS2<-allOnco[SS2]
oncoSS2<-lapply(oncoSS2,function(x) x[which(x[,"FDR"]<1),])

oncoSS2<-do.call('rbind',oncoSS2)
plotOncoSS2<-oncoSS2[,c("Depleted Gene","FEATURE","FDR")]
colnames(plotOncoSS2)<-c("Gene","Feature","FDR")
plotOncoSS2$Feature<-unlist(sapply(plotOncoSS2$Feature,function(x) strsplit(x,"_",fixed=T)[[1]][length(strsplit(x,"_",fixed=T)[[1]])]))



plotOnco32<-SM[SM$Target%in%SS2,]
plotOnco32$Feature<-c("CN","CN","CN")
plotOnco32<-plotOnco32[,c("Target","Feature","CFE.ANOVA.FDR")]
colnames(plotOnco32)<-c("Gene","Feature","FDR")
plotOncoSS2<-cbind(plotOncoSS2,source="Score2")

plotOnco32<-cbind(plotOnco32,source="Score")
pOncoSS2<-rbind(plotOncoSS2,plotOnco32)

#try using corrplot package:
pOncoSS2$row<-paste(pOncoSS2$source,pOncoSS2$Feature,sep="-")
plotMatSS2<-matrix(NA,nrow=4,ncol=3)
rownames(plotMatSS2)<-c("Score2-Expr","Score-Expr","Score2-CN","Score-CN"
                      )
colnames(plotMatSS2)<-SS2
for(i in 1:nrow(pOncoSS2)){
  plotMatSS2[pOncoSS2[i,"row"],pOncoSS2[i,"Gene"]]<-(pOncoSS2[i,"FDR"])
}
#plotMat<-(-log(plotMat,10))
plotMatSS2[is.na(plotMatSS2)]<-0
plotMatSS2[plotMatSS2!=0]<-1
pdf(paste0(outputdata,"/CompareMarkers_OncoScoreScore2.pdf"),useDingbats = FALSE)
corrplot(plotMatSS2,is.cor = FALSE,p.mat=plotMatSS2,sig.level=2,insig="blank",tl.col="black")
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
  "Paralog","SyntheticLethal","AddictionND","OncogenicAddiction_Act","SelfAddiction_Variant",
  "SelfAddiction_CN","LoF_mutation_SL","LoF_other_SL","SelfAddiction_Protein","MSI")
cat(paste("Percent of asscociations with DFA: ",sum(DFRdata$SLgroup%in%useDFAgroups)/nrow(DFRdata)),file=paste0(outputdata,"/DMAlandscape.txt"),sep="\n",append=T)
cat(paste("Number of addiction DFA: ",sum(DFRdata$SLgroup%in%c("Addiction","SelfAddiction_Mutation","SelfAddiction_Expression",
                                                               "AddictionND","OncogenicAddiction_Act","SelfAddiction_Variant",
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

DFRSL<-which(rowSums(DFRsplit[,c("LoFSLmut","LoFSLother","Paralog","SyntheticLethal")])>0)


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
vennSL[["OtherSL"]]<-unique(rownames(DFRsplit)[DFRsplit[,"SyntheticLethal"]==1])
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
vennset3[["SyntheticLethal"]]<-unique(rownames(DFRsplit3)[rowSums(DFRsplit3[,c("LoFSLmut","LoFSLother","Paralog","SyntheticLethal")])>0])
vennset3[["Composite"]]<-unique(rownames(DFRsplit3)[rowSums(DFRsplit3[,c("MSI","Tsubtype")])>0])
venn.diagramCP(vennset3,paste0(outputdata,"Venn3Groups3.pdf"),imagetype="pdf",height=5,width=5,units="in")
vennAdd3<-list()
vennAdd3[["OtherAddiction"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"Target_Activating"]==1])
vennAdd3[["DriverGeneAddiction"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"Target_cancerDriverAct"]==1])
vennAdd3[["SelfAddiction"]]<-unique(rownames(DFRsplit3)[rowSums(DFRsplit3[,c("SelfAddiction_CN","SelfAddiction_Expression","SelfAddiction_Mutation","SelfAddiction_Variant","SelfAddiction_Protein")])>0])
venn.diagramCP(vennAdd3,paste0(outputdata,"VennAddiction3.pdf"),imagetype="pdf",height=5,width=5,units="in")

vennSL3<-list()
vennSL3[["OtherSL"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"SyntheticLethal"]==1])
vennSL3[["LoFmutationSL"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"LoFSLmut"]==1])
vennSL3[["DriverGeneSL"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"LoFSLother"]==1])
vennSL3[["Paralogs"]]<-unique(rownames(DFRsplit3)[DFRsplit3[,"Paralog"]==1])
venn.diagramCP(vennSL3,paste0(outputdata,"VennSL3.pdf"),imagetype="pdf",height=5,width=5,units="in")

