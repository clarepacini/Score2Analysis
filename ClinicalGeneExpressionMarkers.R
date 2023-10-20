
#################################################################################################
#                                                                                               #
#   Load annotations, tumor cluster assignments and gene markers                                #
#                                                                                               #
#                                                                                               #
#################################################################################################

ctypeMapSubtype<-read.csv(paste0(inputdata,"cTypeMapSubtype.csv"),header=F,stringsAsFactors = F)

load(paste0(inputdata,"allpartitionsC.Rdata"))

load(file=paste0(inputdata,"CellOut_metadata.rdata"))
#from RF:
load(file=paste0(inputdata,"RFbms2.Rdata"))
allmarkersRF<-RFbms2
#from DE:
load(file=paste0(inputdata,"DEgenes.Rdata"))
load(file=paste0(inputdata,"DEgenes2.Rdata"))
allmarkersDE<-c(DEgenes,DEgenes2)
#from NN:
allmarkersNN<-readRDS(file=paste0(inputdata,"GeneExpBMList_AllNN.Rds"))

#################################################################################################
#                                                                                               #
#   Combine markers from Neural Net, Random Forest & Differential Expression                    #
#                                                                                               #
#                                                                                               #
#################################################################################################

allmarkersSM<-list()
markers3<-list()
nsubtype<-list()
outnames<-c()
for(i in 1:length(allmarkersNN)){
  secondnames<-ctypeMapSubtype[match(names(allmarkersNN)[i],make.names(ctypeMapSubtype[,1])),7]
  RFnames<-ctypeMapSubtype[match(names(allmarkersNN)[i],make.names(ctypeMapSubtype[,1])),9]
  t1<-allmarkersRF[[RFnames]]
  t2<-allmarkersDE[grep(secondnames,names(allmarkersDE),ignore.case=T)]
  t3<-allpartitionsC[[secondnames]]
  nsubtype[[i]]<-length(unique(t3))
  compareList<-list()
  if(!is.null(t1)&length(t2)>0){
    
    compareList[[1]]<-allmarkersNN[i][[1]]
    compareList[[2]]<-t1
    compareList[[3]]<-unlist(lapply(t2,function(x) x$ID))
    compareList<-lapply(compareList,function(x) x[!is.na(x)])
    allmarkersSM[[names(allmarkersNN)[i]]]<-union(intersect(compareList[[1]],compareList[[3]]),intersect(compareList[[2]],compareList[[3]]))

  }else{
    allmarkersSM[[names(allmarkersNN)[i]]]<-allmarkersNN[[i]]
    
  }
  
}



saveRDS(allmarkersSM,file=paste0(outputdata,"GeneExpBMList.Rds"))

glistsize<-unlist(lapply(allmarkersSM,length))
markerdata<-cbind(unlist(allmarkersSM),unlist(sapply(1:length(allmarkersSM),function(x) rep(names(allmarkersSM)[x],glistsize[x]))))
colnames(markerdata)<-c("geneSymbol","CancerType")
write.table(markerdata,file=paste0(outputdata,"/SupplementaryTable6.tsv"),quote=F,sep="\t",row.names=F)


load(paste0(inputdata,"manifestPostCLGROUP.Rdata"))
useSubtypes<-unique(manifest$CL_GROUP)



#################################################################################################
#                                                                                               #
#   Load Progeny Pathway analysis of gene expression markers                                                 #
#                                                                                               #
#                                                                                               #
#################################################################################################

load(paste0(inputdata,"allpathway.Rdata"))
normPath<-apply(allpathway,2,function(x) (x-mean(x))/sd(x))
rownames(normPath)<-make.names(rownames(normPath))
subtypescore<-NULL
subtypetest<-NULL

for(i in 1:length(useSubtypes)){
  sel<-which(make.names(ctypeMapSubtype[,1])==useSubtypes[i])
  if(length(sel)>0){
    cmpid<-ctypeMapSubtype[sel[1],4]
    print(cmpid)
    #extract cell lines:
 
    
    
    sel2<-which(names(allpartitionsC)==make.names(ctypeMapSubtype[sel[1],7]))
    print(names(allpartitionsC)[sel2])
    try({
      partitionin<-allpartitionsC[sel2][[1]]
      npart<-length(unique(partitionin))
      
      subtypeP<-sapply(1:npart,function(x) apply(normPath[intersect(names(partitionin)[partitionin==x],rownames(normPath)),],2,median))
      #do wilcox test here for differences in pathway activities between patient clusters:
      combPT<-combn(1:npart,2)
      test1<-list()
      for(j in 1:ncol(combPT)){
        
        test1[[j]]<-sapply(1:14,function(x) wilcox.test(normPath[intersect(names(partitionin)[partitionin==combPT[1,j]],rownames(normPath)),x],
                                                        normPath[intersect(names(partitionin)[partitionin==combPT[2,j]],rownames(normPath)),x])$p.value)
      }
      temp<-do.call('rbind',test1)
      temp<-matrix(p.adjust(temp,"BH"),ncol=14)
      colnames(temp)<-colnames(normPath)
      rownames(temp)<-apply(combPT,2,function(x) paste(useSubtypes[i],paste(x,collapse="-")))
      subtypetest[[i]]<-temp
      
      colnames(subtypeP)<-paste(useSubtypes[i],1:npart,sep="_")
      if(is.null(subtypescore)){
        subtypescore<-subtypeP
      }else{
        subtypescore<-cbind(subtypescore,subtypeP)
      }
    })
    
    
  }
}
names(subtypetest)<-useSubtypes
allsubtypetest<-do.call('rbind',subtypetest)
write.table(allsubtypetest,file=paste0(outputdata,"/wilcox_progeny.tsv"),quote=FALSE,sep="\t",col.names=NA)


median(setdiff(unique(unlist(lapply(allmarkersSM,length))),0))


