
#################################################################################################
#                                                                                               #
#   Load data and patient prevalence information                                                #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################


load(file=paste0(inputdata,"GeneP1Upv2.rdata"))
load(file=paste0(inputdata,"GeneP2Upv2.rdata"))
load(file=paste0(inputdata,"GeneP1Downv2.rdata"))
load(file=paste0(inputdata,"GeneP2Downv2.rdata"))



load(file=paste0(inputdata,"TCGAbems.Rdata"))
#add in composite prevalence estimated for the transcriptional subtypes and the progeny subtypes?
load(file=paste0(inputdata,"TsubtypeP.rdata"))
load(file=paste0(inputdata,"TissueTypeColors.Rdata"))

ctypeMapSubtype<-read.csv(paste0(inputdata,"cTypeMapSubtype.csv"),header=F,stringsAsFactors = F)
ctypeMapSubtype$ctype<-make.names(ctypeMapSubtype[,1])

#need to update, need to get the new lists 
load(file=paste0(outputdata,"/GdfSplit.Rdata"))



GdfSplit$M<-unlist(sapply(GdfSplit$MARKER,function(x) strsplit(x,"_",fixed=TRUE)[[1]][1]))
GdfSplit$ID<-paste0(GdfSplit$ID,GdfSplit$M)
GdfSplit$ID<-paste0(GdfSplit$ID,GdfSplit$MARKER_TYPE)
GdfSplit$ASSOCIATION_EFFECT<-gsub(" ","",GdfSplit$ASSOCIATION_EFFECT)

load(file=paste0(inputdata,"CN_gain_PP.Rdata"))
load(file=paste0(inputdata,"CN_loss_PP.Rdata"))

#################################################################################################
#                                                                                               #
#   Convert all patient names to patient IDs                                                    #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################


TCGAmut<-TCGAbems
for(i in 1:length(TCGAbems)){
  temp<-unlist(sapply(colnames(TCGAbems[[i]]),function(x) paste(strsplit(x,"-",fixed=T)[[1]][1:3],collapse="-" )))
  colnames(TCGAmut[[i]])<-temp
}

CNgainMat<-CN_gain_PP
CNlossMat<-CN_loss_PP

for(i in 1:length(CN_gain_PP)){
  temp<-unlist(sapply(colnames(CN_gain_PP[[i]]),function(x) paste(strsplit(x,"-",fixed=T)[[1]][1:3],collapse="-" )))
  colnames(CNgainMat[[i]])<-temp
}
for(i in 1:length(CN_loss_PP)){
  temp<-unlist(sapply(colnames(CN_loss_PP[[i]]),function(x) paste(strsplit(x,"-",fixed=T)[[1]][1:3],collapse="-" )))
  colnames(CNlossMat[[i]])<-temp
}


#################################################################################################
#                                                                                               #
#   Get marker prevalence for each of the priority targets (Cancer Type specific)               #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################

ctypes<-unique(GdfSplit$ctype)

#Get all the patient data available for cancer types included in our data set - this is total pan cancer number of patients.


SampleNos<-list()
for(i in ctypes){
  output<-Get_PatientMarkers(i)
  if(!is.null(output)){
    SampleNos[[i]]<-output$Nsamples
  }
}
TotalSamples<-sum(unlist(SampleNos))
Get_CTtargetPrevalence<-function(GdfInput,ExcludeMarkers=NULL,IncGroup,ExcNoNetwork=FALSE,ScoreThresh=0,incGexp=TRUE,GexpOnly=FALSE){
  PrevMat<-list()
  GdfSplit<-GdfInput
  for(i in ctypes){
    idx<-which(GdfSplit$ctype==i)
    T1<-GdfSplit[idx,]
    
    Tall<-T1
    Tall<-Tall[Tall$GROUP%in%IncGroup,]
    Tall<-Tall[!Tall$MARKER%in%ExcludeMarkers,]
    if(ExcNoNetwork){
      Tall<-Tall[Tall$RWRscore>0,]
    }
    Tall<-Tall[Tall$PRIORITYL3>ScoreThresh,]
 
    tcgaid<-ctypeMapSubtype[match(i,make.names(ctypeMapSubtype[,1])),8]
    #Check whether it is tissue or subtype and then get the correct matrices.
    if(tcgaid%in%names(GeneP1Upv2)){
      ExprInc<-GeneP1Upv2[[tcgaid]]
      ExprDec<-GeneP1Downv2[[tcgaid]]
    }else{
  
      ExprInc<-GeneP2Upv2[[tcgaid]]
      ExprDec<-GeneP2Downv2[[tcgaid]]
    }
    colnames(ExprInc)<-unlist(sapply(colnames(ExprInc),function(x) paste0(strsplit(x,".",fixed=T)[[1]][1:3],collapse="-")))
    colnames(ExprDec)<-unlist(sapply(colnames(ExprDec),function(x) paste0(strsplit(x,".",fixed=T)[[1]][1:3],collapse="-")))
    
    tcgaid<-ctypeMapSubtype[match(i,make.names(ctypeMapSubtype[,1])),3]
    CNgain<-CNgainMat[[tcgaid]]
    CNloss<-CNlossMat[[tcgaid]]
    #mutect, varscan:
    tcgaid<-ctypeMapSubtype[match(i,make.names(ctypeMapSubtype[,1])),3]
    Mutmat1<-TCGAmut[[paste0(tcgaid,".mutect")]]
    Mutmat2<-TCGAmut[[paste0(tcgaid,".varscan")]]
    msamples<-intersect(colnames(Mutmat1),colnames(Mutmat2))
    mgenes<-intersect(rownames(Mutmat1),rownames(Mutmat2))
    Mmat<-Mutmat1[mgenes,msamples]+Mutmat2[mgenes,msamples]
    
    Mmat<-(Mmat==2)+0
    if(length(msamples)>1&!is.null(ExprInc)){

      Pmat<-NULL
      if(GexpOnly){
        try(Pmat<-Get_PrevalenceMatrix(Tall[Tall[,"MARKER_TYPE"]=="expr",],Mmat,CNgain,CNloss,ExprUp=ExprInc,ExprDown=ExprDec))
      }else{
        if(incGexp){
          try(Pmat<-Get_PrevalenceMatrix(Tall,Mmat,CNgain,CNloss,ExprUp=ExprInc,ExprDown=ExprDec))
        }else{
          try(Pmat<-Get_PrevalenceMatrix(Tall,Mmat,CNgain,CNloss))
        }
      }
      PrevMat[[i]]<-Pmat
     
    }else{
      PrevMat[[i]]<-NULL
    
    }
    
  }
  return(PrevMat)

}



PrevMat<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers=NULL,IncGroup=c(1,2),ScoreThresh=0)
PrevMatNoTP53<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(1,2),ScoreThresh = 0)
PrevMatG2<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers=NULL,IncGroup=c(2),ScoreThresh=0)
PrevMatNoTP53G2<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(2),ScoreThresh = 0)
PrevMatnoGE<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers=NULL,IncGroup=c(1,2),ScoreThresh=0,incGexp = FALSE)
PrevMatNoTP53noGE<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(1,2),ScoreThresh = 0,incGexp = FALSE)
PrevMatnoGE2<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers=NULL,IncGroup=c(2),ScoreThresh=0,incGexp = FALSE)
PrevMatNoTP53noGE2<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(2),ScoreThresh = 0,incGexp = FALSE)
PrevMatnoGE1<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers=NULL,IncGroup=c(1),ScoreThresh=0,incGexp = FALSE)
PrevMatNoTP53noGE1<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(1),ScoreThresh = 0,incGexp = FALSE)
PrevMatNoTP53G1<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(1),ScoreThresh = 0,incGexp = TRUE)

Pall<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers=NULL,IncGroup=c(1,2,3),ScoreThresh=0)
Pallno53<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(1,2,3),ScoreThresh=0)
Pallno53noGE<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(1,2,3),ScoreThresh=0,incGexp = FALSE)

GexpG1<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(1),ScoreThresh = 0,GexpOnly=TRUE)
GexpG12<-Get_CTtargetPrevalence(GdfSplit,ExcludeMarkers="TP53_mut",IncGroup=c(1,2),ScoreThresh = 0,GexpOnly=TRUE)
#################################################################################################
#                                                                                               #
#   Get overview statistics for patient prevalence                                              #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################



#All groups:
Get_PatientPrev<-function(Pmat,removeTop=0,Nsamples){
  if(removeTop!=0){
    for(i in 1:length(Pmat)){
      temp<-Pmat[[i]]
      check<-sum(rowSums(temp,na.rm=T)/ncol(temp)<removeTop)>0
      if(sum(check)>0){
        Pmat[[i]]<-Pmat[[i]][check,]
      }else{
        Pmat[[i]]<-matrix(0,nrow=1,ncol=ncol(temp))
      }
    }
    
  }

    NoPatientsOneTarget<-lapply(Pmat,function(x) try(colSums(x,na.rm=T)>0))
    PropPatientperCT<-lapply(NoPatientsOneTarget,function(x) try(sum(x)/length(x)))
    PropPatientpc<-sum(unlist(NoPatientsOneTarget)/Nsamples)
  
  TargetsPP<-lapply(Pmat,function(x) try(colSums(x,na.rm=T)))
  TargetsPPperCT<-lapply(TargetsPP,mean)
  TargetsPPpc<-mean(unlist(TargetsPPperCT))
  return(list(PropCT=PropPatientperCT,PropPC=PropPatientpc,AvgTppCT=TargetsPPperCT,AvgTppPC=TargetsPPpc))
}


AllG<-Get_PatientPrev(Pall,Nsamples=TotalSamples)
AllG_no53<-Get_PatientPrev(Pallno53,Nsamples=TotalSamples)
AllG_no53noGE<-Get_PatientPrev(Pallno53noGE,Nsamples=TotalSamples)


cat(paste("Proportion patients pancan all: ",AllG$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=F)
cat(paste("Proportion patients pancan all, no P53: ",AllG_no53$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Proportion patients pancan all, no P53 no GExp: ",AllG_no53noGE$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)

cat(paste("Avg Targets PP pancan all: ",AllG$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Avg Targets PP pancan all, no P53: ",AllG_no53$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Avg Targets PP pancan all, no P53 no GExp: ",AllG_no53noGE$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)

AllG12<-Get_PatientPrev(PrevMat,Nsamples=TotalSamples)
AllG12_no53<-Get_PatientPrev(PrevMatNoTP53,Nsamples=TotalSamples)
cat(paste("Proportion patients pancan G1,2: ",AllG12$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Proportion patients pancan G1,2, no P53: ",AllG12_no53$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Avg Targets PP pancan G1,2: ",AllG12$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Avg Targets PP pancan G1,2, no P53: ",AllG12_no53$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)

G12_allnoGE<-Get_PatientPrev(PrevMatnoGE,Nsamples=TotalSamples)
G12_no53noGE<-Get_PatientPrev(PrevMatNoTP53noGE,Nsamples=TotalSamples)


cat(paste("Proportion patients pancan G12 no Gene Exp: ",G12_allnoGE$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Proportion patients pancan G12, no P53 no Gene Exp: ",G12_no53noGE$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)

cat(paste("Avg Targets PP pancan G12, no Gene Exp: ",G12_allnoGE$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Avg Targets PP pancan G12, no P53 no Gene Exp: ",G12_no53noGE$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)



#GROUP 2:
G2_all<-Get_PatientPrev(PrevMatG2,Nsamples=TotalSamples)
G2_no53<-Get_PatientPrev(PrevMatNoTP53G2,Nsamples=TotalSamples)


cat(paste("Proportion patients pancan G2: ",G2_all$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Proportion patients pancan G2, no P53: ",G2_no53$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)

cat(paste("Avg Targets PP pancan G2: ",G2_all$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Avg Targets PP pancan G2, no P53: ",G2_no53$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)

G2_allnoGE<-Get_PatientPrev(PrevMatnoGE2,Nsamples=TotalSamples)
G2_no53noGE<-Get_PatientPrev(PrevMatNoTP53noGE2,Nsamples=TotalSamples)


cat(paste("Proportion patients pancan G2 no Gene Exp: ",G2_allnoGE$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Proportion patients pancan G2, no P53 no Gene Exp: ",G2_no53noGE$PropPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)

cat(paste("Avg Targets PP pancan G2, no Gene Exp: ",G2_allnoGE$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)
cat(paste("Avg Targets PP pancan G2, no P53 no Gene Exp: ",G2_no53noGE$AvgTppPC),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)

#Group1:

G1_no53<-Get_PatientPrev(PrevMatNoTP53G1,Nsamples=TotalSamples)

G1_no53noGE<-Get_PatientPrev(PrevMatNoTP53noGE1,Nsamples=TotalSamples)
#Gene Expression only:
Gexp_G1<-Get_PatientPrev(GexpG1,Nsamples=TotalSamples)
Gexp_G12<-Get_PatientPrev(GexpG12,Nsamples=TotalSamples)

#All g 1 and 2 exc TP53 - G12 no TP53 no Gene expression (should give increase across groups 1 and 2 from Gene exp)
G1_no53$PropCT$Glioma<-0
G1_no53noGE$PropCT$Glioma<-0
G1_no53$PropCT$Kidney.Carcinoma<-0
G1_no53noGE$PropCT$Kidney.Carcinoma<-0

IncG12<-unlist(AllG12_no53$PropCT)-unlist(G12_no53noGE$PropCT[names(AllG12_no53$PropCT)])
IncNoGexp_G2<-(unlist(G12_no53noGE$PropCT)-unlist(G1_no53noGE$PropCT[names(G12_no53noGE$PropCT)]))
unlist(AllG12_no53$PropCT)-unlist(G1_no53$PropCT)-IncNoGexp_G2
AllG2diff<-unlist(AllG12_no53$PropCT)-unlist(G1_no53$PropCT[names(AllG_no53$PropCT)])
AllG2diffnoGE<-unlist(G12_no53noGE$PropCT)-unlist(G1_no53noGE$PropCT[names(G12_no53noGE$PropCT)])
#bar 1:
IncGexp_g1<-unlist(G1_no53$PropCT)-unlist(G1_no53noGE$PropCT[names(G1_no53$PropCT)])
#bar 2:
#group 1 and 2 all
IncGexp_g2<-IncG12-IncGexp_g1[names(unlist(IncG12))]

#get total number of unique group 2 targets with patient prevalence of associated marker:
inPatient<-lapply(PrevMatNoTP53noGE2,function(x) rownames(x)[rowSums(x)>0])
inPatientTarget<-sapply(names(inPatient),function(x) sapply(inPatient[[x]],function(y) strsplit(y,x)[[1]][1]))
cat(paste("Total number pancan G2, no P53 no Gene Exp: ",length(unique(unlist(inPatientTarget)))),file=paste0(outputdata,"/Prevalencelandscape.txt"),sep="\n",append=T)

#################################################################################################
#                                                                                               #
#   Plot target prevalence overview information on cancer type basis                            #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################
G1_allnoGE<-Get_PatientPrev(PrevMatnoGE1,Nsamples=TotalSamples)
G1_no53noGE<-Get_PatientPrev(PrevMatNoTP53noGE1,Nsamples=TotalSamples)

GdataCT<-data.table(cancer_type=names(AllG$PropCT),AllG=unlist(AllG$PropCT),Gno53noExp=unlist(AllG_no53noGE$PropCT),stringsAsFactors = FALSE)
GdataCT<-GdataCT[rowSums(GdataCT[,2:3])>0,]
GdataCT$cancer_type<-factor(GdataCT$cancer_type,levels=unlist(GdataCT[order(apply(GdataCT[,2:3],1,max),decreasing=F),"cancer_type"]))


PlotG2data<-melt(GdataCT,id="cancer_type")

pdf(paste0(outputdata,"/CancerTypePropsAll.pdf"),useDingbats = FALSE)
print(ggplot(PlotG2data,aes(x=cancer_type,y=value,fill=variable))+geom_bar(stat="identity",position="dodge")+coord_flip()+theme_bw()+ scale_fill_brewer(palette="Accent"))
dev.off()

G2dataCTnge<-data.table(cancer_type=names(G2_allnoGE$PropCT),AllG2=unlist(G2_all$PropCT),G2_no53noGE=unlist(G2_no53noGE$PropCT),
                        stringsAsFactors = FALSE)
G2dataCTnge<-G2dataCTnge[rowSums(G2dataCTnge[,2:3])>0,]
G2dataCTnge$cancer_type<-factor(G2dataCTnge$cancer_type,levels=unlist(G2dataCTnge[order(apply(G2dataCTnge[,2:3],1,max),decreasing=F),"cancer_type"]))


PlotG2datange<-melt(G2dataCTnge,id="cancer_type")

pdf(paste0(outputdata,"/CancerTypePropsG2.pdf"),useDingbats = FALSE)
print(ggplot(PlotG2datange,aes(x=cancer_type,y=value,fill=variable))+geom_bar(stat="identity",position="dodge")+coord_flip()+theme_bw()+ scale_fill_brewer(palette="Accent"))
dev.off()

GdataCTAll<-data.frame(cancer_type=names(AllG$PropCT),Prop=unlist(AllG$PropCT),Name="AllG"
                       ,stringsAsFactors = FALSE)
GdataCTAll<-rbind(GdataCTAll,data.frame(cancer_type=names(AllG_no53noGE$PropCT),Prop=unlist(AllG_no53noGE$PropCT),Name="All_noGExpno53"))
GdataCTAll<-rbind(GdataCTAll,data.frame(cancer_type=names(G2_no53noGE$PropCT),Prop=unlist(G2_no53noGE$PropCT),Name="G2_noGExpno53"))
sortGdata<-GdataCTAll[GdataCTAll$Name=="G2_noGExpno53",]
sortGdata<-sortGdata[order(sortGdata$Prop,decreasing=F),]
GdataCTAll$cancer_type<-factor(GdataCTAll$cancer_type,levels=c(setdiff(GdataCTAll$cancer_type,sortGdata$cancer_type),unique(unlist(sortGdata[,"cancer_type"]))))
GdataCTAll$Name<-factor(GdataCTAll$Name,levels=c("G2_noGExpno53","All_noGExpno53","AllG"))
pdf(paste0(outputdata,"/Figure5a.pdf"),useDingbats = FALSE)
print(ggplot(GdataCTAll[GdataCTAll$Name%in%c("AllG","All_noGExpno53","G2_noGExpno53"),],aes(x=cancer_type,y=Prop,fill=Name))+geom_bar(stat="identity",position="dodge")+coord_flip()+theme_bw()+ scale_fill_brewer(palette="Accent"))
dev.off()

GdataCTAll<-rbind(GdataCTAll,data.frame(cancer_type=names(G1_allnoGE$PropCT),Prop=unlist(G1_allnoGE$PropCT),Name="G1_noGExp"))
GdataCTAll<-rbind(GdataCTAll,data.frame(cancer_type=names(G1_no53$PropCT),Prop=unlist(G1_no53$PropCT),Name="G1_no53"))
GdataCTAll<-rbind(GdataCTAll,data.frame(cancer_type=names(G1_no53noGE$PropCT),Prop=unlist(G1_no53noGE$PropCT),Name="G1_noGExpno53"))
GdataCTAll<-rbind(GdataCTAll,data.frame(cancer_type=names(G12_no53noGE$PropCT),Prop=unlist(G12_no53noGE$PropCT),Name="G12_noGExpno53"))
GdataCTAll<-rbind(GdataCTAll,data.frame(cancer_type=names(G12_allnoGE$PropCT),Prop=unlist(G12_allnoGE$PropCT),Name="G12_noGExp"))
#find missing cancer types to add in:
setdiff(names(G1_no53noGE),names(G12_no53noGE))
useCT<-union(names(G1_no53noGE$PropCT),names(G12_no53noGE$PropCT))
G1_no53noGE$PropCT$Kidney.Carcinoma<-0
GdataCTAll<-rbind(GdataCTAll,data.frame(cancer_type=names(G12_allnoGE$PropCT),Prop=unlist(G12_no53noGE$PropCT)-unlist(G1_no53noGE$PropCT[names(G12_no53noGE$PropCT)]),Name="G12_noGExpInc"))

Psub1<-GdataCTAll[GdataCTAll$Name%in%c("G1_noGExpno53","G12_noGExpInc"),]
Psub1<-Psub1[order(Psub1$Prop,decreasing=F),]
Psub1$cancer_type<-factor(Psub1$cancer_type,levels=c(setdiff(GdataCTAll$cancer_type,sortGdata$cancer_type),unique(unlist(sortGdata[,"cancer_type"]))))
Psub1$Name<-factor(Psub1$Name,levels=c("G12_noGExpInc","G1_noGExpno53"))
pdf(paste0(outputdata,"/CancerTypePropsNoGE_1and2.pdf"),useDingbats = FALSE)
print(ggplot(Psub1,aes(x=cancer_type,y=Prop,fill=Name))+geom_bar(stat="identity",position="stack")+coord_flip()+theme_bw()+ scale_fill_brewer(palette="Paired"))
dev.off()

Psub2<-GdataCTAll[GdataCTAll$Name%in%c("G1_noGExpno53","G12_noGExpno53"),]
Psub2<-Psub2[order(Psub2$Prop,decreasing=F),]
Psub2$cancer_type<-factor(Psub2$cancer_type,levels=unique(unlist(Psub2[,"cancer_type"])))

pdf(paste0(outputdata,"/CancerTypePropsNoGE_1and2no53.pdf"),useDingbats = FALSE)
print(ggplot(Psub2,aes(x=cancer_type,y=Prop,fill=Name))+geom_bar(stat="identity",position="dodge")+coord_flip()+theme_bw()+ scale_fill_brewer(palette="Accent"))
dev.off()

GdataCTAll<-rbind(GdataCTAll,data.frame(cancer_type=names(AllG12_no53$PropCT),Prop=unlist(AllG12_no53$PropCT)-unlist(G1_no53$PropCT[names(AllG12_no53$PropCT)]),Name="G12_GExpInc"))

Psub3<-GdataCTAll[GdataCTAll$Name%in%c("G1_no53","G12_GExpInc"),]
Psub3<-Psub3[order(Psub3$Prop,decreasing=F),]
Psub3$cancer_type<-factor(Psub3$cancer_type,levels=c(setdiff(GdataCTAll$cancer_type,sortGdata$cancer_type),unique(unlist(sortGdata[,"cancer_type"]))))
Psub3$Name<-factor(Psub3$Name,levels=c("G12_GExpInc","G1_no53"))
pdf(paste0(outputdata,"/CancerTypePropsGE_1and2.pdf"),useDingbats = FALSE)
print(ggplot(Psub3,aes(x=cancer_type,y=Prop,fill=Name))+geom_bar(stat="identity",position="stack")+coord_flip()+theme_bw()+ scale_fill_brewer(palette="Pastel1"))
dev.off()

#################################################################################################
#                                                                                               #
#   Plot target prevalence information for each cancer type                                     #
#                                                                                               #                                                                                           
#                                                                                               #
#################################################################################################

DFAgroups<-unique(GdfSplit$SLgroup)
DFAgroups<-DFAgroups[c(1:3,6:10,12,17)]

pCol<-c("blue","purple","orange")
names(pCol)<-c("Mutation","Expression","CopyNumber")
for(i in ctypes){

  try(Plot_PrevalenceHM(PrevMat=PrevMat[[i]],GdfSplit,plotCol=pCol,outputname=paste0(outputdata,i),DFA=NULL))
 
}

