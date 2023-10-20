
load(paste0(inputdata,"cancerDrivers.rdata"))

load(paste0(inputdata,"/07_EssMatrix_bDepletionsB2.Rdata"))
cmp<-read.csv(paste0(inputdata,"model_list_latest.csv"),header=T,stringsAsFactors = FALSE)
cmp$cancer_type<-make.names(cmp$cancer_type)

load(paste0(outputdata,"DFRdata.Rdata"))
load(paste0(outputdata,"DFRsplit.Rdata"))
load(paste0(inputdata,"PPInet.Rdata"))

load(paste0(inputdata,"PPIigraph.Rdata"))
degreelist<-igraph::degree(PPIigraph)


#################################################################################################
#                                                                                               #
#   Read in RWR biomarker analysis for 500 selective pancancer dependencies.                    #
#   Includes all possible biomarkers in biomarker classes A-D                                   #                                                                                           
#                                                                                               #
#################################################################################################
allrwr<-read.table(paste0(inputdata,"/AllBM_CombinedPCrwr_aovout.tsv"),header=T,stringsAsFactors = FALSE,sep="\t")

allrwr$id<-paste(allrwr$Depleted.Gene,allrwr$FEATURE,sep="-")
allrwr<-allrwr[allrwr$AnovaClass!="None",]
allrwr<-allrwr[!is.na(allrwr$AnovaClass),]

#################################################################################################
#                                                                                               #
#   Annotate most stringent features from DFA analysis for 500 most selective dependencies.     #
#                                                                                               #
#                                                                                               #
#################################################################################################
Usegroups<-c("SyntheticLethal","LoF_other_SL","LoF_mutation_SL","Paralog","SelfAddiction_Mutation","SelfAddiction_Expression","SelfAddiction_Variant","OncogenicAddiction_Act",
             "SelfAddiction_CN","AddictionND")
DFRdataT<-cbind(DFRdata,DFRsplit[,c(1:11,13:14)])
DFRuse<-DFRdataT[rowSums(DFRdataT[,35:47])>0,]
DFRuse$AddictionGroup<-rowSums(DFRuse[,c("SelfAddiction_Mutation","SelfAddiction_Expression","SelfAddiction_Variant","Target_cancerDriverAct",
                                         "SelfAddiction_CN","SelfAddiction_Protein","Target_Activating")])
DFRuse$SLset<-rowSums(DFRuse[,c("SyntheticLethal","Paralog","LoFSLmut","LoFSLother")])
DFRuse$SelfAddictionGroup<-rowSums(DFRuse[,c("SelfAddiction_Mutation","SelfAddiction_Expression","SelfAddiction_Variant",
                                         "SelfAddiction_CN","SelfAddiction_Protein")])
DFRuse$ExcSelfAddictionGroup<-rowSums(DFRuse[,c("Target_cancerDriverAct",
                                         "Target_Activating")])

DFRuse$RWRscore<-0
for(i in 1:nrow(DFRuse)){
  sel<-paste(DFRuse[i,"TARGET"],DFRuse[i,"FEATURE"],sep="-")
  if(sel%in%allrwr$id){
    DFRuse[i,"RWRscore"]<-allrwr[match(sel,allrwr$id),"RWRscore"]
  }
}
rwrfreq<-table(DFRuse[,"RWRscore"])

write.table(DFRuse,file=paste0(outputdata,"/SupplementaryTable9.tsv"),quote=F,sep="\t",row.names=F)

DFRnonzero<-DFRuse[DFRuse$RWRscore!=0,]
cat(paste("Percent DFA with signif RWR score:",sum(rwrfreq[2:5])/sum(rwrfreq)),file=paste0(outputdata,"RWRlandscape.txt"),append=F,sep="\n")
cat(paste("Number dependencies with signif RWR score:",length(unique(DFRnonzero$TARGET))),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")

#################################################################################################
#                                                                                               #
#   Test for differences between Addiction and Synthetic Lethal relationships.                  #
#                                                                                               #
#                                                                                               #
#################################################################################################

cmat<-matrix(0,nrow=2,ncol=5)
rownames(cmat)<-c("Addiction","SyntheticLethal")
colnames(cmat)<-c("100","75","50","25","0")
cmatSA<-cmat
rownames(cmatSA)<-c("SelfAddiction","SyntheticLethal")
cmatNo<-cmat
rownames(cmatNo)<-c("ExcSelfAddiction","SyntheticLethal")
for(i in 1:nrow(DFRuse)){
  slgroup<-DFRuse[i,"SLset"]
  #if(slgroup%in%c("SL","LoF_other_SL","LoF_mutation_SL","Paralog")){
  if(slgroup>0){
    rgroup<-as.character(DFRuse[i,"RWRscore"])
    cmat["SyntheticLethal",rgroup]<-cmat["SyntheticLethal",rgroup]+1
    cmatNo["SyntheticLethal",rgroup]<-cmatNo["SyntheticLethal",rgroup]+1
  }

  slgroup<-DFRuse[i,"AddictionGroup"]
  if(slgroup>0){
    rgroup<-as.character(DFRuse[i,"RWRscore"])
    cmat["Addiction",rgroup]<-cmat["Addiction",rgroup]+1
  }
  slgroup<-DFRuse[i,"SelfAddictionGroup"] 

  if(slgroup>0){
      rgroup<-as.character(DFRuse[i,"RWRscore"])
      cmatSA["SelfAddiction",rgroup]<-cmatSA["SelfAddiction",rgroup]+1
      
      
  }
  slgroup<-DFRuse[i,"ExcSelfAddictionGroup"]
   
  if(slgroup>0){
      rgroup<-as.character(DFRuse[i,"RWRscore"])
      cmatNo["ExcSelfAddiction",rgroup]<-cmatNo["ExcSelfAddiction",rgroup]+1
      
      
    }
  
}


DFAaddSL<-contingencyPlot(cmat,filename="CorrPlot_DFR_AddictionSL",outputdata)
DFAnoSASL<-contingencyPlot(cmatNo,filename="CorrPlot_DFR_NoSelfAddictionSL",outputdata)
cat(paste("Class 1 P-value Addiction SL:",DFAaddSL$pvalTable[1,1]),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")
cat(paste("Class 2 P-value Addiction SL:",DFAaddSL$pvalTable[1,2]),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")
cat(paste("Class 3 P-value Addiction SL:",DFAaddSL$pvalTable[1,3]),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")
cat(paste("Class 1 P-value Addiction SL exc Self:",DFAnoSASL$pvalTable[1,1]),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")

cat(paste("Global P-value Addiction SL:",DFAaddSL$globalPval),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")

cmatpd<-reshape2::melt(cmat)
cmatpd<-group_by(cmatpd, Var1) %>% mutate(percent = value/sum(value))
cmatpd$Var2<-factor(cmatpd$Var2,levels=c("100","75","50","25","0"),labels=c("A","B","C","D","None"))

pdf(paste0(outputdata,"/AddictionSL.pdf"),useDingbats = FALSE)
print(ggplot(data=cmatpd,aes(fill=Var1,x=Var2,y=percent))+geom_bar(stat="identity",position="dodge")+scale_fill_brewer()+theme_bw())
dev.off()
png(paste0(outputdata,"/AddictionSL.png"),width=4,height=4,units="in",res=300)
print(ggplot(data=cmatpd,aes(fill=Var1,x=Var2,y=percent))+geom_bar(stat="identity",position="dodge")+scale_fill_brewer()+theme_bw())
dev.off()

#################################################################################################
#                                                                                               #
#   Look for differences between Driver and non-driver features.                                #
#                                                                                               #
#                                                                                               #
#################################################################################################



DFRuse$IsDriver<-FALSE
DFRuse[DFRuse$M%in%cancerDrivers[[1]],"IsDriver"]<-TRUE

cdmat<-matrix(0,nrow=2,ncol=5)
rownames(cdmat)<-c("Driver","NonDriver")
colnames(cdmat)<-c("100","75","50","25","0")
for(i in 1:nrow(DFRuse)){
  
  if(DFRuse[i,"IsDriver"]){
    rgroup<-as.character(DFRuse[i,"RWRscore"])
    cdmat["Driver",rgroup]<-cdmat["Driver",rgroup]+1
    
  }else{
    rgroup<-as.character(DFRuse[i,"RWRscore"])
    cdmat["NonDriver",rgroup]<-cdmat["NonDriver",rgroup]+1
    
  }
  
}

DFAdriverND<-contingencyPlot(cdmat,filename="CorrPlot_DFR_DriverND",outputdata)
cat(paste("Global P-value Driver Non-driver:",DFAdriverND$globalPval),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")
cat(paste("P-value Driver Non-driver, zero score:",DFAdriverND$pvalTable["Driver","0"]),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")

cdmatpd<-reshape2::melt(cdmat)
cdmatpd<-group_by(cdmatpd, Var1) %>% mutate(percent = value/sum(value))
cdmatpd$Var2<-factor(cdmatpd$Var2,levels=c("100","75","50","25","0"),labels=c("A","B","C","D","None"))

pdf(paste0(outputdata,"/DriverNonDriver.pdf"),useDingbats = FALSE)
print(ggplot(data=cdmatpd,aes(fill=Var1,x=Var2,y=percent))+geom_bar(stat="identity",position="dodge")+scale_fill_brewer(palette=4)+theme_bw())
dev.off()
png(paste0(outputdata,"/DriverNonDriver.png"),width=4,height=4,units="in",res=300)
print(ggplot(data=cdmatpd,aes(fill=Var1,x=Var2,y=percent))+geom_bar(stat="identity",position="dodge")+scale_fill_brewer(palette=4)+theme_bw())
dev.off()

MaxSet<-DFRuse[DFRuse$RWRscore==100,]
s1<-unique(MaxSet[MaxSet$IsDriver,"TARGET"])
s2<-unique(MaxSet[MaxSet$IsDriver==FALSE,"TARGET"])

cat(paste("Unique number targets with driver feature, 100 RWR",length(s1)),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")
cat(paste("Number DFAs with driver feature, 100 RWR",sum(MaxSet$IsDriver)),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")

cat(paste("Number DFAs with non-driver feature, 100 RWR",sum(!MaxSet$IsDriver)),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")
#################################################################################################
#                                                                                               #
#   Test for enrichment of ontology terms in feature set with no RWR score and not a driver     #
#                                                                                               #
#                                                                                               #
#################################################################################################

#the selfaddiction expression shouldnt be in the group. 
DFRzero<-DFRuse[DFRuse$RWRscore==0,]
DFRzero<-DFRzero[!DFRzero$SLgroup%in%c("SelfAddiction_Expression","Paralog"),]
DFRzero<-DFRzero[!DFRzero$IsDriver,]
OtherMarkers<-unique(DFRzero$M)
cat(paste("Number DFAs with non-driver feature, no network score",length(OtherMarkers)),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")

enrichRes<-gost(OtherMarkers,user_threshold = 1e-03)
KEGGenrichNMF<-enrichRes$result[enrichRes$result$source=="KEGG",]
ReactomeEnrichNMF<-enrichRes$result[enrichRes$result$source=="REAC",]
GOBPEnrichNMF<-enrichRes$result[enrichRes$result$source=="GO:BP",]

GOBPEnrichNMF$logpval<-(-1*log(GOBPEnrichNMF$p_value,10))

GOBPEnrichNMF$term_name<-factor(GOBPEnrichNMF$term_name,levels=GOBPEnrichNMF$term_name[order(GOBPEnrichNMF$logpval,decreasing = F)])
pdf(paste0(outputdata,"/GOBPzeroRWR.pdf"),useDingbats = FALSE)
print(ggplot(data=GOBPEnrichNMF,aes(x=term_name,y=logpval,fill="red"))+scale_fill_manual(values=c("darkblue"))+geom_bar(stat="identity")+theme_bw()+coord_flip())
dev.off()
png(paste0(outputdata,"/GOBPzeroRWR.png"),width=4,height=4,units="in",res=300)
print(ggplot(data=GOBPEnrichNMF,aes(x=term_name,y=logpval,fill="purple"))+scale_fill_manual(values=c("darkblue"))+geom_bar(stat="identity")+theme_bw()+coord_flip())
dev.off()

cat(paste("Enrichment p-value cell adhesion",GOBPEnrichNMF[GOBPEnrichNMF$term_name=="cell adhesion","p_value"]),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")
cat(paste("Enrichment p-value cell differentiation",GOBPEnrichNMF[GOBPEnrichNMF$term_name=="cell differentiation","p_value"]),file=paste0(outputdata,"RWRlandscape.txt"),append=T,sep="\n")

#################################################################################################
#                                                                                               #
#   Analyse extended set of features with biomarker class A-D for 500 selective dependencies    #
#                                                                                               #
#                                                                                               #
#################################################################################################


allPC<-read.table(paste0(inputdata,"/AllBM_CombinedPC_aovout.tsv"),header=T,stringsAsFactors = FALSE,sep="\t")

allPC$id<-paste(allPC$Depleted.Gene,allPC$FEATURE,sep="-")
allPC<-allPC[allPC$AnovaClass!="None",]
allPC<-allPC[!is.na(allPC$AnovaClass),]
load(file=paste0(inputdata,"paralogList.Rdata"))
load(file=paste0(inputdata,"cancerDrivers.rdata"))
LoFgenes<-cancerDrivers[[1]][which(cancerDrivers[[2]]=="LoF")]

allPC$TARGET<-allPC$Depleted.Gene
allPC$BMtype<-allPC$type

#################################################################################################
#                                                                                               #
#   Association between different biomarker classes with RWR class A and different features     #
#                                                                                               #
#                                                                                               #
#################################################################################################
allPC$DepDegree<-unlist(degreelist[allPC$TstringID])
allPC$BMdegree<-unlist(degreelist[allPC$StringID])
Ndeps<-rowSums(bDepletionsB2,na.rm=T)
allPC$nDepLines<-Ndeps[allPC$Depleted.Gene]
allPC$M<-unlist(sapply(allPC$FEATURE,function(x) strsplit(x,"_",fixed=T)[[1]][1]))
allPC$selfAssoc<-allPC$Depleted.Gene==allPC$M

PCcont<-matrix(0,nrow=4,ncol=4)
rownames(PCcont)<-c("A","B","C","D")
colnames(PCcont)<-c("100","75","50","25")

for(i in 1:nrow(allPC)){
  Aclass<-allPC[i,"AnovaClass"]
  Rclass<-as.character(allPC[i,"RWRscore"])
  if(!is.na(Aclass)&!is.na(Rclass)&Aclass!="None"){
    PCcont[Aclass,Rclass]<-PCcont[Aclass,Rclass]+1
  }
}


chisquareRes<-contingencyPlot(PCcont,filename="CorrPlot_Map_AovRWRpc",outputdata)
write.table(chisquareRes$pvalTable,file=paste0(outputdata,"/SupplementaryTable10.tsv"),quote=F,sep="\t")
