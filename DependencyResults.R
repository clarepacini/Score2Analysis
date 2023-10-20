
data('EssGenes.DNA_REPLICATION_cons')
data('EssGenes.KEGG_rna_polymerase')
data('EssGenes.PROTEASOME_cons')
data('EssGenes.ribosomalProteins')
data('EssGenes.SPLICEOSOME_cons')
data('BAGEL_essential')
data('BAGEL_nonEssential')
#################################################################################################
#                                                                                               #
#   Load manifests, cell line annotation data and get overview of screened models               #
#                                                                                               #
#                                                                                               #
#################################################################################################

load(paste0(inputdata,"manifestCL.Rdata"))
load(paste0(inputdata,"TissueColours.Rdata"))


manifest$Source<-"Broad"
manifest[grep("SIDM",manifest$INSTITUTE_ID),"Source"]<-"Sanger"
manifest$LIBRARY<-"Avana"
manifest[grep("SIDM",manifest$INSTITUTE_ID),"LIBRARY"]<-"KY"


cmp<-read.csv(paste0(inputdata,'model_list_latest.csv'),header=T,stringsAsFactors = FALSE)

MASTER_LIST<-cmp
cmp2<-cmp
cmp2$model_id<-cmp2$BROAD_ID
cmp<-rbind(cmp,cmp2)


allIDs<-unique(manifest$INSTITUTE_ID)
allIDs[allIDs%in%cmp$BROAD_ID]<-cmp[match(allIDs[allIDs%in%cmp$BROAD_ID],cmp$BROAD_ID),"model_id"]
manifest$CMP_ID<-manifest$INSTITUTE_ID
manifest[manifest$INSTITUTE_ID%in%cmp$BROAD_ID,"CMP_ID"]<-cmp[match(manifest[manifest$INSTITUTE_ID%in%cmp$BROAD_ID,"INSTITUTE_ID"],cmp$BROAD_ID),"model_id"]
manifest$CellLineName<-cmp[match(manifest$CMP_ID,cmp$model_id),"model_name"]

AllCrispr<-manifest[,c("INSTITUTE_ID","TISSUE","CANCER_TYPE")]

write.csv(manifest[,c("SAMPLE_ID","CellLineName","INSTITUTE_ID","LIBRARY","TISSUE","CANCER_TYPE","PROJECTS","FILENAMES","Source","CMP_ID")],file=paste0(outputdata,"SupplementaryTable1.csv"),quote=F)


cat(paste("Number of Cell lines:",length(allIDs)),file=paste0(outputdata,"DataOverview.txt"),append=F,sep="\n")
cat(paste("Number of cancer types:",length(unique(AllCrispr$CANCER_TYPE))),file=paste0(outputdata,"DataOverview.txt"),append=T,sep="\n")


AllCrispr<-unique(AllCrispr)
Parentterms<-sapply(unique(AllCrispr$CANCER_TYPE),function(x) AllCrispr[match(x,AllCrispr$CANCER_TYPE),"TISSUE"])
NumberCT<-sapply(unique(AllCrispr$CANCER_TYPE),function(x) sum(AllCrispr$CANCER_TYPE==x))
Parentterms<-gsub("."," ",Parentterms,fixed=T)
names(Parentterms)<-gsub("."," ",names(Parentterms),fixed=T)

df<-data.frame(name=Parentterms,type=names(Parentterms),value=NumberCT,stringsAsFactors = FALSE)%>%
  mutate(name = as.factor(name) %>% fct_reorder(value, sum)) %>%
  arrange(name, value) %>%
  mutate(type = as.factor(type) %>% fct_reorder2(name, value))
df$type<-factor(df$type,levels=c(levels(df$type),"Other"))
df[df$value<10,"type"]<-"Other"
rownames(df)<-NULL
df<-df%>%group_by(name,type)%>%summarise(value=sum(value))%>%as.data.frame(.)



CancerTypeColours<-unique(AllCrispr$CANCER_TYPE)
names(CancerTypeColours)<-CancerTypeColours
CancerTypeColours<-TissueColours[AllCrispr[match(names(CancerTypeColours),AllCrispr$CANCER_TYPE),"TISSUE"]]
names(CancerTypeColours)<-unique(AllCrispr$CANCER_TYPE)
save(CancerTypeColours,file=paste0(outputdata,"CancerTypeColours.Rdata"))

#################################################################################################
#                                                                                               #
#   Tissue-specific and pancancer core fitness genes overview                                   #
#                                                                                               #
#                                                                                               #
#################################################################################################


load(paste0(inputdata,'/09_ADM_analysisTypes.RData'))


ctypes<-sort(summary(as.factor(ADManalysisTypes)),decreasing=TRUE)
ctypes<-names(which(ctypes>=8))
load(paste0(inputdata,'/07_EssMatrix_qnorm_corrected_logFCs.RData'))
load(paste0(inputdata,'/07_EssMatrix_bDepletionsB2.Rdata'))
bDepletions<-bDepletionsB2
load(paste0(inputdata,'/10_PANCANCER_coreFitness_genes.RData'))
load(paste0(inputdata,'/10_PANCANCER_CommonEssentialsAUC.RData'))



load(paste0(inputdata,"/curated_BAGEL_essential.rdata"))
load(paste0(inputdata,"/curated_BAGEL_nonEssential.rdata"))
load(paste0(inputdata,"/BAGEL_v2_ESSENTIAL_GENES.rdata"))
load(paste0(inputdata,"/histones.RData"))
allgenes<-unique(c(EssGenes.DNA_REPLICATION_cons,EssGenes.KEGG_rna_polymerase,EssGenes.PROTEASOME_cons,
                   EssGenes.ribosomalProteins,EssGenes.SPLICEOSOME_cons,BAGEL_essential,BAGEL_nonEssential,
                   curated_BAGEL_essential,curated_BAGEL_nonEssential,histones))

load(paste0(inputdata,"/allSymbol.Rdata"))

allMap<-makeNameMap(allgenes,allSymbol)
EssGenes.DNA_REPLICATION_cons<-updateNames(EssGenes.DNA_REPLICATION_cons,allMap,uset=FALSE)
EssGenes.KEGG_rna_polymerase<-updateNames(EssGenes.KEGG_rna_polymerase,allMap,uset=FALSE)
EssGenes.PROTEASOME_cons<-updateNames(EssGenes.PROTEASOME_cons,allMap,uset=FALSE)
EssGenes.ribosomalProteins<-updateNames(EssGenes.ribosomalProteins,allMap,uset=FALSE)
EssGenes.SPLICEOSOME_cons<-updateNames(EssGenes.SPLICEOSOME_cons,allMap,uset=FALSE)
BAGEL_essential<-updateNames(BAGEL_essential,allMap,uset=FALSE)
BAGEL_nonEssential<-updateNames(BAGEL_nonEssential,allMap,uset=FALSE)
curated_BAGEL_nonEssential<-updateNames(curated_BAGEL_nonEssential,allMap,uset=FALSE)
curated_BAGEL_essential<-updateNames(curated_BAGEL_essential,allMap,uset=FALSE)
histones<-updateNames(histones,allMap,uset=FALSE)

BEM<-readRDS(file=paste0(inputdata,'ConsBEM.Rds'))
BEMall<-BEM
BEMallM<-BEMall[grep("_mut",rownames(BEMall)),]
TestFeatAll<-rowSums(BEMall,na.rm=T)>4+(rowSums(BEMall,na.rm=T)<(ncol(BEMall)-4))
TestFeatAllM<-rowSums(BEMallM,na.rm=T)>4+(rowSums(BEMallM,na.rm=T)<(ncol(BEMallM)-4))

NC<-vector()

contextSpecific<-list()
CoreContextSpecific<-list()
for (i in 1:length(ctypes)){
  
  ctype<-ctypes[i]
  print(ctype)
  
  load(paste(inputdata,'/09_ADM_',ctype,'_coreFitnessGenes.Rdata',sep=''))
  cellLines<-names(which(ADManalysisTypes==ctypes[i]))
  NC[i]<-length(cellLines)
  
  depleted<-names(which(rowSums(bDepletions[,intersect(cellLines,colnames(bDepletions))])>2))
  

  depleted<-setdiff(depleted,c(PanCancerCoreFitnessGenes,coreFitnessGenes,EssGenes.DNA_REPLICATION_cons,EssGenes.PROTEASOME_cons,EssGenes.SPLICEOSOME_cons,
                               EssGenes.KEGG_rna_polymerase,EssGenes.ribosomalProteins,curated_BAGEL_essential,curated_BAGEL_nonEssential,histones))
  ccs<-setdiff(coreFitnessGenes,c(PanCancerCoreFitnessGenes,EssGenes.DNA_REPLICATION_cons,EssGenes.PROTEASOME_cons,EssGenes.SPLICEOSOME_cons,
                                  EssGenes.KEGG_rna_polymerase,EssGenes.ribosomalProteins,curated_BAGEL_essential,curated_BAGEL_nonEssential,histones))
  contextSpecific[[i]]<-depleted
  CoreContextSpecific[[i]]<-ccs
  ContexSpecESS<-depleted
  save(ContexSpecESS,file=paste(outputdata,'12_CS_ess_',ctype,'.RData',sep=''))
  
}

names(contextSpecific)<-ctypes
names(CoreContextSpecific)<-ctypes
allCS<-sort(unique(unlist(contextSpecific)))
allCSCF<-sort(unique(unlist(CoreContextSpecific)))

CSmatrix<-NULL
CSCFmatrix<-NULL
for (i in 1:length(contextSpecific)){
  CSmatrix<-cbind(CSmatrix,is.element(allCS,contextSpecific[[i]]))    
  CSCFmatrix<-cbind(CSCFmatrix,is.element(allCSCF,CoreContextSpecific[[i]])) 
}

rownames(CSmatrix)<-allCS
colnames(CSmatrix)<-ctypes
rownames(CSCFmatrix)<-allCSCF
colnames(CSCFmatrix)<-ctypes

CSmatrix<-CSmatrix+0
CSCFmatrix<-CSCFmatrix+0
CSmatrix<-CSmatrix[order(rowSums(CSmatrix),decreasing=TRUE),order(colSums(CSmatrix),decreasing=TRUE)]

save(CSmatrix,file=paste0(outputdata,'13_CSmatrix.RData'))
cat(paste("Median number core fitness dependencies (>2 CLs):",median(colSums(CSmatrix))),file=paste0(outputdata,"/DependencyResults.txt"),sep="\n",append=T)
BEM<-CLnameMapping(bDepletions,BEM,annot=cmp)$refdata
#Overview of number of markers that can be tested:
NumberVar<-length(grep("_var",rownames(BEM)))
NumberMut<-length(grep("_mut",rownames(BEM)))
NumberTS<-length(grep("_Celligner",rownames(BEM)))+length(grep("Tsubtype",rownames(BEM)))+length(grep("PAM50",rownames(BEM)))+length(grep("ScMod",rownames(BEM)))+length(grep("CMS",rownames(BEM)))
BEMmut<-BEM[grep("_mut",rownames(BEM)),]
TestFeat<-rowSums(BEM,na.rm=T)>4+(rowSums(BEM,na.rm=T)<(ncol(BEM)-4))
TestFeatMut<-rowSums(BEMmut,na.rm=T)>4+(rowSums(BEMmut,na.rm=T)<(ncol(BEMmut)-4))

#################################################################################################
#                                                                                               #
#   Calculate pan-cancer selectivity (normLRT) for all dependencies                             #
#                                                                                               #
#                                                                                               #
#################################################################################################


#this takes a long time to run and so is loaded from pre-calculation
#normLRTres<-normLRT2(qnorm_corrected_logFCs)
#save(normLRTres,file=paste0(outputdata,"/normLRTresnq.Rdata"))
load(file=paste0(inputdata,"/normLRTresnq.Rdata"))
save(normLRTres,file=paste0(outputdata,"/normLRTresnq.Rdata"))
LRT<-normLRTres$LRT
write.table(LRT,file=paste0(outputdata,"/LRTvalues.tsv"),sep='\t',quote=F,col.names=NA)

summary(LRT)
top<-2000
topLRTgenes<-names(sort(LRT,decreasing=T))[1:top]
BinaryDep<-bDepletions[setdiff(rownames(CSmatrix),CommonEssentialsAUC$cfgenes),]

test<-qnorm_corrected_logFCs
test<-(-1*test)
test<-apply(test,1,minmax)
inputdataNMF<-t(test)
save(inputdataNMF,file=paste0(outputdata,"inputdataNMFnq.rdata"))


nDeps<-rowSums(bDepletions,na.rm=T)
LRT<-LRT[setdiff(names(LRT),union(PanCancerCoreFitnessGenes,CommonEssentialsAUC$cfgenes))]
topLRTgenes<-names(sort(LRT,decreasing=T))
write.table(topLRTgenes[1:500],file=paste0(outputdata,"/Top500normLRT.tsv"),sep='\t',quote=F,row.names=FALSE)
lrtplotdata<-data.frame(xval=1:500,yval=sort(LRT,decreasing=T)[1:500],gene=rep("",500),nDep=nDeps[topLRTgenes[1:500]],stringsAsFactors = FALSE)
lrtplotdata$nGroup<-signif(lrtplotdata$nDep,1)
lrtplotdata[lrtplotdata$nDep<10,'nGroup']<-10
lrtplotdata[lrtplotdata$nDep>200,'nGroup']<-200



pdf(paste0(outputdata,"/normLRTranks_colour.pdf"),width=10,height=5)
#,col=nDep
print(ggplot(aes(x=xval,y=yval,label=gene,color=nGroup),data=lrtplotdata)+geom_point(aes(color=nGroup,size=.7))+geom_text_repel(colour="blue",size=4,max.overlaps = 50)+
        ylab("Selectivity: normLRT")+xlab("Rank")+theme(panel.background = element_rect(fill="white",color="grey",linetype="solid",size=0.5),
                                                        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                                                        colour = "grey"), 
                                                        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                                                        colour = "grey"))+
        ylim(0,2000)+
        scale_color_viridis_c())

#+custom_colors)

dev.off()




#zoom to just top 50 genes:
lrtplotdata<-data.frame(xval=1:50,yval=sort(LRT,decreasing=T)[1:50],gene=topLRTgenes[1:50],nDep=nDeps[topLRTgenes[1:50]],stringsAsFactors = FALSE)
lrtplotdata$nGroup<-signif(lrtplotdata$nDep,1)
lrtplotdata[lrtplotdata$nDep<10,'nGroup']<-10
lrtplotdata[lrtplotdata$nDep>200,'nGroup']<-200

pdf(paste0(outputdata,"/normLRTranks_Zoom_Colour.pdf"),width=10,height=5)
#,col=nDep
print(ggplot(aes(x=xval,y=yval,label=gene,color=nGroup),data=lrtplotdata)+geom_point(aes(color=nGroup,size=1.5))+geom_text_repel(colour="blue",size=4,max.overlaps = 50)+
        ylab("Selectivity: normLRT")+xlab("Rank")+theme(panel.background = element_rect(fill="white",color="grey",linetype="solid",size=0.5),
                                                        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                                                                        colour = "white"), 
                                                        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                                                        colour = "white"))+
        scale_color_viridis_c()
)


dev.off()



#################################################################################################
#                                                                                               #
#   Tissue specific core fitness genes analysis                                                 #
#                                                                                               #
#                                                                                               #
#################################################################################################

colnames(CSCFmatrix)<-gsub("."," ",colnames(CSCFmatrix),fixed=T)

idlabs<-colnames(CSCFmatrix)
ncols<-length(unique(idlabs))

labcols<-TissueColours[idlabs]

TC<-TissueColours
names(TC)<-gsub("."," ",names(TC),fixed=T)
labcols<-TC[colnames(CSCFmatrix)]
write.csv(CSCFmatrix,file=paste0(outputdata,"/SupplementaryTable2.csv"),quote=FALSE,row.names=TRUE)

plotdata<-data.frame(NumberCore=colSums(CSCFmatrix),Tissue=colnames(CSCFmatrix),stringsAsFactors = FALSE)
CRISPR<-make.names(colnames(bDepletions))


names(NC)<-ctypes
pdf(paste0(outputdata,"NumberCSE_NumberCL.pdf"),useDingbats = FALSE)
  plot(NC,plotdata$NumberCore,col=TC[plotdata$Tissue],pch=19,cex=2,xlab="Number of Cell lines",ylab="Number context-specific core fitness genes")
  text(35,1000,paste("Spearman Corr:",format(cor(NC,plotdata$NumberCore,method="spearman"),scientific=FALSE)))
dev.off()

plotdata$Tissue<-factor(plotdata$Tissue,levels=plotdata[order(plotdata$NumberCore),"Tissue"])

pdf(paste0(outputdata,"NumberCoreSpecificEssentials.pdf"))
  print(ggplot(plotdata,aes(fill=Tissue,y=NumberCore,x=Tissue))+geom_bar(stat="identity")+coord_flip()+scale_fill_manual(values=labcols)+theme(legend.position="none"))
dev.off()


cat(paste("Median number context specific essentials:",median(plotdata$NumberCore)),file=paste0(outputdata,"/DependencyResults.txt"),sep="\n",append=T)
cat(paste("Range number context specific essentials:",range(plotdata$NumberCore)),file=paste0(outputdata,"/DependencyResults.txt"),sep="\n",append=T)

cat(paste("Number pancancer core essentials:",length(PanCancerCoreFitnessGenes)),file=paste0(outputdata,"/DependencyResults.txt"),sep="\n",append=T)

#################################################################################################
#                                                                                               #
#   Bootstrap Tissue specific core fitness genes analysis                                       #
#                                                                                               #
#                                                                                               #
#################################################################################################



load(file=paste0(inputdata,"NcoreF300.Rdata"))
load(file=paste0(inputdata,"NcoreF500.Rdata"))
load(file=paste0(inputdata,"Ntissue500.Rdata"))
load(file=paste0(inputdata,"Ntissue300.Rdata"))
#plot how many tissue could run analysis for per sampled dataset:
Ntissueplot<-rbind(data.frame(NT=Ntissue300,Nsamples=rep(300,length(Ntissue300))),
                   data.frame(NT=Ntissue500,Nsamples=rep(500,length(Ntissue500))))
Ntissueplot$Nsamples<-as.factor(Ntissueplot$Nsamples)
pdf(paste0(outputdata,"/Boxplot_NumberTissue_Downsample.pdf"),useDingbats = FALSE,width=3,height=4)
print(ggplot(data=Ntissueplot,aes(x=Nsamples,y=NT))+geom_boxplot(outlier.size=0.8)+
  geom_jitter(color="black", size=0.8, alpha=1,width=0.25,height=0)+theme_bw()+xlab("Number cell lines sampled")+ylab("Number tissue identified core-fitness genes")+ylim(0,21)+geom_hline(yintercept=20, linetype="dashed"))
dev.off()

plot300<-NULL
for(i in 1:length(NcoreF300)){
  if(!is.null(NcoreF300[[i]])){
    temp<-data.frame(NcF<-NcoreF300[[i]],Tissue=rep(names(NcoreF300)[i],length(NcoreF300[[i]])))
    plot300<-rbind(plot300,temp)
  }
}
plot300$Tissue<-gsub("."," ",fixed=T,plot300$Tissue)
plot300$Tissue<-factor(plot300$Tissue,levels=plotdata[order(plotdata$NumberCore),"Tissue"])

pdf(paste0(outputdata,"NumberCoreSpecificEssentials_DS300.pdf"),useDingbats = FALSE,width=4,height=6)
print(ggplot(plot300,aes(fill=Tissue,y=NcF....NcoreF300..i..,x=Tissue))+geom_boxplot(outlier.size = 0.4)+
        coord_flip()+scale_fill_manual(values=labcols)+theme_bw()+theme(legend.position="none")+
        geom_jitter(color="black", size=0.4, alpha=0.7,width=0.25,height=0)+ylab("Number tissue-specific core Fitness Genes"))
dev.off()

plot500<-NULL
for(i in 1:length(NcoreF500)){
  if(!is.null(NcoreF500[[i]])){
    temp<-data.frame(NcF<-NcoreF500[[i]],Tissue=rep(names(NcoreF500)[i],length(NcoreF500[[i]])))
    plot500<-rbind(plot500,temp)
  }
}
plot500$Tissue<-gsub("."," ",fixed=T,plot500$Tissue)
plot500$Tissue<-factor(plot500$Tissue,levels=plotdata[order(plotdata$NumberCore),"Tissue"])

pdf(paste0(outputdata,"NumberCoreSpecificEssentials_DS500.pdf"),useDingbats = FALSE,width=4,height=6)
print(ggplot(plot500,aes(fill=Tissue,y=NcF....NcoreF500..i..,x=Tissue))+geom_boxplot(outlier.size = 0.4)+
        coord_flip()+scale_fill_manual(values=labcols)+theme_bw()+theme(legend.position="none")+
        geom_jitter(color="black", size=0.4, alpha=0.7,width=0.25,height=0)+ylab("Number tissue-specific core Fitness Genes"))
dev.off()

#################################################################################################
#                                                                                               #
#   Tissue selective dependency analysis                                                        #
#                                                                                               #
#                                                                                               #
#################################################################################################



tissuetypes<-names(table(AllCrispr$TISSUE))[table(AllCrispr$TISSUE)>9]

n_minus1<-c()
n_minus2<-c()
n_minus3<-c()
CEfc_tissue<-c()
TSGfc_tissue<-list()
RefCF<-list()
ncell_tissue<-c()
allcore<-c()
specificT<-list()
inputFCdata<-qnorm_corrected_logFCs[setdiff(rownames(qnorm_corrected_logFCs),PanCancerCoreFitnessGenes),]
for(i in tissuetypes){
  tissueID<-i
  modelIDs<-AllCrispr[AllCrispr$TISSUE==i,"INSTITUTE_ID"]
  CEfc_tissue[i]<-median(as.vector(qnorm_corrected_logFCs[intersect(rownames(qnorm_corrected_logFCs),PanCancerCoreFitnessGenes),modelIDs]))
  load(paste(inputdata,'/09_ADM_',tissueID,'_coreFitnessGenes.Rdata',sep=''))
  allcore<-c(allcore,coreFitnessGenes)
  sel<-inputFCdata[setdiff(rownames(inputFCdata),coreFitnessGenes),modelIDs]
  tempfc<-data.frame(FC=as.vector(qnorm_corrected_logFCs[intersect(rownames(qnorm_corrected_logFCs),coreFitnessGenes),modelIDs]),tissue=i,type="TSCF",stringsAsFactors = FALSE)
  RefFC<-data.frame(FC=as.vector(qnorm_corrected_logFCs[intersect(rownames(qnorm_corrected_logFCs),BAGEL_essential),modelIDs]),tissue=i,type="REF",stringsAsFactors = FALSE)
  
  TSGfc_tissue[[i]]<-tempfc
  RefCF[[i]]<-RefFC
  ncell_tissue[i]<-length(modelIDs)
  t1<-apply(sel,1,function(x) sum(x<(-1.5),na.rm=T)>1)
  t1<-names(t1[t1>0])
  specificT[[i]]<-t1
  t2<-apply(sel,1,function(x) sum(x<(-2),na.rm=T)>1)
  t2<-names(t2[t2>0])
  t3<-apply(sel,1,function(x) sum(x<(-3),na.rm=T)>1)
  t3<-names(t3[t3>0])
  n_minus1<-c(n_minus1,length(setdiff(t1,t2)))
  n_minus2<-c(n_minus2,length(setdiff(t2,t3)))
  n_minus3<-c(n_minus3,length(t3))
}
names(n_minus1)<-tissuetypes
names(n_minus2)<-tissuetypes
names(n_minus3)<-tissuetypes
fcdata<-data.frame(tissueType=tissuetypes,n_1=n_minus1,n_2=n_minus2,n_3=n_minus3,stringsAsFactors = FALSE)
fcdata<-melt(fcdata)
fcdatan1<-fcdata[fcdata$variable=="n_1",]
fcdata$tissueType<-factor(make.names(fcdata$tissueType),levels=names(ncell_tissue)[order(ncell_tissue,decreasing=F)])
fcdata$variable<-factor(fcdata$variable,levels=c("n_3","n_2","n_1"))

gout4<-ggplot(data=as.data.frame(fcdata),aes(x=tissueType,y=value,fill=variable))+geom_bar(position="stack",stat="identity",colour="grey")+theme_bw()+coord_flip()+scale_fill_manual(values=c("orange","blue","lightblue"),labels=c("Minus3","Minus2","Minus1"))
pdf(paste0(outputdata,"/NumberByFC_tissue.pdf"),useDingbats = FALSE)
print(gout4)
dev.off()

#################################################################################################
#                                                                                               #
#   Cumulative Tissue selective dependency analysis                                             #
#                                                                                               #
#                                                                                               #
#################################################################################################
CumRes<-list()
CumulativePoints<-seq(0,1,0.1)
for(i in tissuetypes){
  tissueID<-i
  modelIDs<-AllCrispr[AllCrispr$TISSUE==i,"INSTITUTE_ID"]
  load(paste(inputdata,'/09_ADM_',tissueID,'_coreFitnessGenes.Rdata',sep=''))
  #allcore<-c(allcore,coreFitnessGenes)
  sel<-inputFCdata[setdiff(rownames(inputFCdata),coreFitnessGenes),modelIDs]
  
  IsDepleted<-(sel<(-1.5))+0
  PropDep<-rowSums(IsDepleted,na.rm=T)/ncol(IsDepleted)
  AtOne<-PropDep[PropDep!=0]
  
  CumRes[[i]]<-data.frame(tissue=i,xvalues=CumulativePoints,NumberGenes=sapply(CumulativePoints,function(x) sum(AtOne>x)),stringsAsFactors = FALSE)
  
  
}

pgenes_10<-lapply(CumRes,function(x) x[1,"NumberGenes"]/sum(x[,"NumberGenes"]))
summary(unlist(pgenes_10))
CdataAll<-do.call('rbind',CumRes)

TColours<-TissueColours
names(TColours)<-make.names(names(TissueColours))
gout5<-ggplot(data=as.data.frame(CdataAll),aes(x=xvalues,y=NumberGenes,colour=tissue))+geom_line()+
  scale_colour_manual(values=TColours[tissuetypes])+
  theme_bw()+xlab("Percent of cell lines depleted")+ylab("Number of Genes")
pdf(paste0(outputdata,"/CumulativeDepletions_tissue.pdf"),useDingbats = FALSE,width=6,height=4)
print(gout5)
dev.off()

#################################################################################################
#                                                                                               #
#   Comparison tissue core versus selective dependencies                                        #
#                                                                                               #
#                                                                                               #
#################################################################################################



corefc<-do.call('rbind',TSGfc_tissue)
Reffc<-do.call('rbind',RefCF)
fcdistn<-rbind(corefc,Reffc)
fcdistn$tissue<-factor(make.names(as.character(fcdistn$tissue)),levels=names(ncell_tissue)[order(ncell_tissue,decreasing=F)])
save(fcdistn,file=paste0(outputdata,"/fcdistn.Rdata"))
labcolmn<-labcols
names(labcolmn)<-make.names(names(labcolmn))
gout5<-ggplot(data=fcdistn,aes(x=tissue,y=FC,fill=type))+geom_violin()+theme_bw()+coord_flip()+scale_fill_manual(values=c("REF"='red',"TSCF"='purple'))
pdf(paste0(outputdata,"/Distn_core_tissue.pdf"),useDingbats = FALSE)
print(gout5)
dev.off()

allsets<-union(allcore,PanCancerCoreFitnessGenes)
number_specific<-lapply(specificT,function(x) setdiff(x,allsets))
onlyone<-sapply(1:length(number_specific),function(x) setdiff(number_specific[[x]],unlist(number_specific[setdiff(1:length(number_specific),x)])))
names(onlyone)<-names(number_specific)
