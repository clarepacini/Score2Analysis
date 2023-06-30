AdamTissue<-function(Binary,annot,returnTissue=FALSE){
  CCRlabels<-data.frame(name=colnames(Binary),sidm=annot[match(colnames(Binary),annot$model_id),"tissue"],bid=annot[match(colnames(Binary),annot$BROAD_ID),"tissue"],stringsAsFactors = FALSE)
  
  CCRlabels[is.na(CCRlabels[,2]),2]<-0
  CCRlabels[is.na(CCRlabels[,3]),3]<-0
  utissues<-unique(annot$tissue)
  ntissues<-length(utissues)
  print(ntissues)
  coreCCR<-list()
  j=1
  tissuenames<-c()
  for(i in 1:ntissues){
    print(utissues[i])
    selCL<-CCRlabels[CCRlabels[,2]==utissues[i]|CCRlabels[,3]==utissues[i],1]
    selCL<-intersect(selCL,colnames(Binary))
    if(length(selCL)>9){
      coreCCR[[j]]<-ADAMwrapper(Binary[,selCL])
      j=j+1
      tissuenames<-c(tissuenames,utissues[i])
    }
  }
  names(coreCCR)<-tissuenames
  
  BinaryTissue<-matrix(0,nrow=nrow(Binary),ncol=length(coreCCR))
  rownames(BinaryTissue)<-rownames(Binary)
  colnames(BinaryTissue)<-names(coreCCR)
  for(i in 1:length(coreCCR)){
    BinaryTissue[coreCCR[[i]],i]<-1
  }
  PCcoreCCR<-ADAMwrapper(BinaryTissue)
  if(returnTissue){
    return(list(Tissue=BinaryTissue,Core=PCcoreCCR))
  }else{
    return(PCcoreCCR)}
}

contingencyPlot<-function(cmat,filename="",dir.Results){
  chisq<-chisq.test(cmat)
  round(chisq$residuals, 3)
  chisq$p.value
  print(paste("Global p-value",chisq$p.value))
  pvalsCT<-pnorm(-1*abs(chisq$stdres))
  if(nrow(cmat)>1){
    pdf(paste0(dir.Results,filename,".pdf"),useDingbats = FALSE)
    corrplot(chisq$residuals, is.cor = FALSE,p.mat=pvalsCT,sig.level=0.05,insig="blank",tl.col="black")
    dev.off()
    
    png(paste0(dir.Results,filename,".png"),width=5,height=5,units="in",res=300)
    corrplot(chisq$residuals, is.cor = FALSE,p.mat=pvalsCT,sig.level=0.05,insig="blank",tl.col="black")
    dev.off()
  }
  return(list(pvalTable=pvalsCT,globalPval=chisq$p.value))
  
}
get_SplitMarkers<-function(GLOBAL,splitVal="//"){
  GdfSplit<-GLOBAL
  cmarkers<-grep(splitVal,GdfSplit[,"MARKER"],fixed=T)
  smarkers<-grep(splitVal,GdfSplit[,"MARKER"],fixed=T,invert = T)
  
  if(length(cmarkers)>0){
    splitM<-list()
    count<-1
    for(k in cmarkers){
      splitM[[count]]<-msplitF(GdfSplit[k,],"MARKER",splitVal=splitVal)
      count<-count+1
    }
    
    if(is.list(splitM)){
      allCMarkers<-do.call(rbind,splitM)
      GdfSplit<-rbind(GdfSplit[smarkers,],allCMarkers)
    }else{
      GdfSplit<-GdfSplit[smarkers,]
      for(k in 1:ncol(splitM)){
        temp<-splitM[,k]
        GdfSplit<-rbind(GdfSplit,unique(temp))
      }
      
    }
  }
  for(i in 1:ncol(GdfSplit)){
    
    GdfSplit[,i]<-unlist(GdfSplit[,i])
  }
  if("PPI_min"%in%colnames(GdfSplit)&"PPI_distance"%in%colnames(GdfSplit)){
     GdfSplit[,"PPI_min"]<-GdfSplit[,"PPI_distance"]
  }
  GdfSplit$MARKER_TYPE<-gsub(" ","",GdfSplit$MARKER_TYPE)
  GdfSplit$MARKER<-gsub(" ","",GdfSplit$MARKER)
  GdfSplit$MARKER_TYPE<-tolower(GdfSplit$MARKER_TYPE)
  GdfSplit$MARKER_TYPE<-gsub("cna","cn",GdfSplit$MARKER_TYPE)
  
  return(GdfSplit)
}

get_DMA<-function(PriorityResults,paraloglist,LoFgenes,markerColumn
                       ,cancerDrivers,splitMarkers=TRUE){
  
  #need to separate markers by commas as well for same profile
  GLOBAL<-PriorityResults
  
  #Separate markers so each row is one target and one marker:
  
  if(splitMarkers){
    GLOBAL<-get_SplitMarkers(GLOBAL)}
  
  GLOBAL$markerName<-unlist(sapply(GLOBAL[,markerColumn],function(x) strsplit(x,"_",fixed=T)[[1]][1]))
  GLOBAL$ID<-paste0(GLOBAL$TARGET,GLOBAL$ctype)
  GLOBAL$SLgroup<-""
  
  #Flag dual markers - do not assign to DFA group
  sel<-grep("&",GLOBAL[,markerColumn])
  GLOBAL[sel,"SLgroup"]<-"Dual"
  
  #Split any markers with same binary pattern. Now that have marked dual markers
  if(splitMarkers){
    GLOBALSplit<-get_SplitMarkers(GLOBAL[GLOBAL$SLgroup!="Dual",],splitVal=",")
    GLOBAL<-rbind(GLOBALSplit,GLOBAL[GLOBAL$SLgroup=="Dual",])}
  GLOBAL$markerName<-unlist(sapply(GLOBAL[,markerColumn],function(x) strsplit(x,"_",fixed=T)[[1]][1]))
  
  
  #Self oncogenic addiction:
  TT1<-c("SelfAddiction_Mutation","SelfAddiction_Expression","SelfAddiction_CN","SelfAddiction_Variant","SelfAddiction_Protein")
  TargetTypeMatrix<-matrix(0,nrow=nrow(GLOBAL),ncol=5,dimnames=list(GLOBAL$ID,TT1))
  
  TargetLoc<-which(GLOBAL$PPI_min==0)
  checkfor<-c("_mut","_Expr","_CN","_var","_Prot")
  for(i in TargetLoc){
    
    markers<-unlist(GLOBAL[i,markerColumn])
    um<-unlist(sapply(unlist(strsplit(markers,"//",fixed=TRUE)),function(x) strsplit(x,",",fixed=TRUE)))
    if(markerColumn=="NetworkMARKER"){
      bmtype<-unlist(strsplit(GLOBAL[,"BMType"],"//",fixed=TRUE))
      um<-sapply(1:length(um),function(x) paste(um[x],bmtype[x],collapse="_"))
      um<-gsub("_variant","_var",um)
    }
    um<-gsub(" ","",um)
    um<-gsub("_..._var","_var",um)
    um<-gsub("_...._var","_var",um)
    #instead of this way do split by "_" and use bmtype to assign annotations/columns in TargetTypeMatrix
    if(length(grep("_var",markers))>0){
      um<-paste0(GLOBAL[i,"TARGET"],"_var")
    }
    checknames<-paste0(GLOBAL[i,"TARGET"],checkfor)
    temp<-sapply(checknames,function(x) is.element(x,um))
    TargetTypeMatrix[i,]<-temp+0
    GLOBAL[i,"SLgroup"]<-paste(TT1[temp],collapse=",")
    
    
  }
  
  #Paralogs:
  #need to split markers by // otherwise this wont work. 
  GLOBALP<-GLOBAL
  GLOBALP$paralogF<-sapply(GLOBALP$TARGET,function(x) which(unlist(lapply(paraloglist,function(y) x%in%y))))
  #GLOBAL$ms<-unlist(sapply(GLOBAL[,"MARKER"],function (y) strsplit(y,"_")[[1]][1]))
  GLOBALP$paralogM<-unlist(apply(GLOBALP,1,function(x) length(intersect(c(unlist(sapply(unlist(strsplit(unlist(x[markerColumn]),'//',fixed=T)),function (y) strsplit(y,"_"))),x["TARGET"]),unlist(paraloglist[unlist(x["paralogF"])])))>1))
  GLOBALP$paralogs<-unlist(apply(GLOBALP,1,function(x) paste(intersect(setdiff(unlist(sapply(unlist(strsplit(unlist(x[markerColumn]),'//',fixed=T)),function (y) strsplit(y,"_")[[1]])),x["TARGET"]),unlist(plist[unlist(x["paralogF"])])),collapse=",")))
  depSign<-unlist(GLOBALP[,"ASSOCIATION_EFFECT"])
  depSign<-gsub(" ","",depSign)
  #Decreased Dependency with gene expression or protein intensity marker:
  ParalogID<-which(GLOBALP$paralogM=="TRUE"&depSign=="DecreasedDep.")
  ExprCheck<-grep("_Expr|_Prot|_CN",GLOBALP[,markerColumn])
  ParalogID<-intersect(ParalogID,ExprCheck)
  #Add any Increased dependency on a LoF mutation or variant:
  ParalogID2<-which(GLOBALP$paralogM=="TRUE"&depSign=="IncreasedDep.")
  
  type<-sapply(unlist(GLOBALP[,markerColumn]),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
  type<-gsub(" ","",type)
  type<-type%in%c("mut","var")
  mLoFd<-GLOBALP[,"markerName"]%in%LoFgenes
  LoFCheck<-intersect(which(mLoFd&type),ParalogID2)
  ParalogID<-union(ParalogID,LoFCheck)
  Paralog<-rep(0,nrow(GLOBAL))
  Paralog[ParalogID]<-1
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Paralog)
  GLOBAL[ParalogID,"SLgroup"]<-sapply(ParalogID,function(x) 
    ifelse(GLOBAL[x,"SLgroup"]!="",paste(GLOBAL[x,"SLgroup"],"Paralog",sep=","),"Paralog"))
  
  
  sel<-1:nrow(GLOBAL)
  
  #Synthetic lethal with Tumour suppressor loss.
  LoFSLmut<-rep(0,nrow(GLOBAL))
  LoFSLother<-rep(0,nrow(GLOBAL))
  LoFActother<-rep(0,nrow(GLOBAL))
  
  for(i in sel){
    markers<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][1])
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    type<-gsub(" ","",type)
    type<-type%in%c("mut","var")
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    mLoFd<-markers%in%LoFgenes
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="IncreasedDep."
    if(is.na(type)){
      check=FALSE
    }else{
      check<-mLoFd&&dep&&type}
    
    if(sum(check)>0){
      LoFSLmut[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"LoF_mutation_SL",sep=","),"LoF_mutation_SL")
    }
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    type<-gsub(" ","",type)
    type<-type%in%c("Expr","CN","Prot")
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="DecreasedDep."
    if(is.na(type)){
      check<-FALSE
    }else{
      check<-mLoFd&&dep&&type}
    
    if(sum(check)>0){
      LoFSLother[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"LoF_other_SL",sep=","),"LoF_other_SL")
    }
    
    dep<-depSign=="IncreasedDep."
    if(is.na(type)){
      check<-FALSE
    }else{
      check<-mLoFd&&dep&&type}
    
    if(sum(check)>0){
      #LoFActother[i]<-1
      #GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"LoF_other_Addiction",sep=","),"LoF_other_Addiction")
    }
    
  }
  #TargetTypeMatrix<-cbind(TargetTypeMatrix,LoFSLmut,LoFSLother,LoFActother)
  TargetTypeMatrix<-cbind(TargetTypeMatrix,LoFSLmut,LoFSLother)
  #MSI
  MSI<-rep(0,nrow(GLOBAL))
  MSISub<-grep("MSI_Composite",GLOBAL[,markerColumn],ignore.case = T)
  MSI[MSISub]<-1
  TargetTypeMatrix<-cbind(TargetTypeMatrix,MSI)
  GLOBAL[MSISub,"SLgroup"]<-sapply(MSISub,function(x) ifelse(GLOBAL[x,"SLgroup"]=="","MSI",GLOBAL[x,"SLgroup"]))
  
  #Transcriptional subtypes
  TSub<-grep("Tsubtype",GLOBAL[,markerColumn])
  CMS<-grep("CMS",GLOBAL[,markerColumn])
  Breastsub<-grep("ScMod",GLOBAL[,markerColumn])
  Breastsub2<-grep('PAM50',GLOBAL[,markerColumn])
  celligner<-grep('celligner',GLOBAL[,markerColumn])
  composite<-grep('_Composite',GLOBAL[,markerColumn])
  MSI<-grep("MSI_Composite",GLOBAL[,markerColumn])
  composite<-setdiff(composite,MSI)
  TSub<-unique(c(TSub,CMS,Breastsub,Breastsub2,celligner,composite))
  Tsubtype<-rep(0,nrow(GLOBAL))
  Tsubtype[TSub]<-1
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Tsubtype)
  GLOBAL[TSub,"SLgroup"]<-sapply(TSub,function(x) ifelse(GLOBAL[x,"SLgroup"]=="","Transcriptional_Subtype",GLOBAL[x,"SLgroup"]))
  
  
  
  
  #Oncogenic addiction
  cDriver<-cbind(cancerDrivers[[1]],cancerDrivers[[2]])
  cDriverAct<-cDriver[cDriver[,2]=="Act",]
  
  cDriverAmb<-cDriver[cDriver[,2]=="ambiguous",]
  MarkerList<-sapply(GLOBAL[,markerColumn],function(x) sapply(strsplit(x,"//",fixed=TRUE),function(y) unlist(strsplit(y,"_",fixed=TRUE))))
  MarkerList<-lapply(MarkerList,function(x) ifelse(is.matrix(x),x[,1],x))
  
  #oncogenic activating
  
  sel<-1:nrow(GLOBAL)
  
  Target_cancerDriverAct<-rep(0,nrow(GLOBAL))
  for(i in sel){
    marker<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][1])
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    mAct<-marker%in%cDriverAct
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="IncreasedDep."
    check<-mAct&&dep
    
    if(sum(check)>0){
      Target_cancerDriverAct[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"OncogenicAddiction_Act",sep=","),"OncogenicAddiction_Act")
    }
  }
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Target_cancerDriverAct)
  
  
  sel<-1:nrow(GLOBAL)
  
  Target_MetabolicM<-rep(0,nrow(GLOBAL))
  
  for(i in sel){
    marker<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][1])
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    type<-gsub(" ","",type)
    type<-type%in%c("Met","met")
    
    check<-type
    if(is.na(check)){check<-FALSE}
    
    if(sum(check)>0){
      Target_MetabolicM[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"Metabolite",sep=","),"Metabolite")
    }
    
  }
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Target_MetabolicM)
  #activating addiction other
  
  
  Target_Activating<-rep(0,nrow(GLOBAL))
  sel<-1:nrow(GLOBAL)
  for(i in sel){
    
    marker<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][1])
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    type<-gsub(" ","",type)
    type<-type%in%c("Expr","Prot","CN")
    mAmb<-!(marker%in%c(cancerDrivers[[1]]))
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="IncreasedDep."
    check<-type&&dep&&mAmb
    if(is.na(check)){check<-FALSE}
    if(sum(check)>0){
      Target_Activating[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"AddictionND",sep=","),"AddictionND")
    }
    
  }
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Target_Activating)
  
  
 
  sel<-1:nrow(GLOBAL)
  SL<-rep(0,nrow(GLOBAL))
  #this is okay because we are keeping one dependency-feature pair per line
  for(i in sel){
    
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    type<-gsub(" ","",type)
    type<-type%in%c("Expr","CN","Prot")
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="DecreasedDep."
    check<-type&&dep
    
    if(sum(check)>0){
      SL[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"SL",sep=","),"SL")
    }
  }
  TargetTypeMatrix<-cbind(TargetTypeMatrix,SL)
  
  
  
  #sel<-which(GLOBAL$SLgroup%in%c("","OncogenicAddiction_Act","OncogenicAddiction_Amb"))
  
  
  
  
  GLOBAL<-as.data.frame(GLOBAL)
  
  for(i in 1:ncol(GLOBAL)){
    
    GLOBAL[,i]<-unlist(GLOBAL[,i])
  }
  
  return(list(DFAresults=GLOBAL,TargetTypeMatrix=TargetTypeMatrix))
}


get_DFAgroup<-function(PriorityResults,paraloglist,LoFgenes,markerColumn
                      ,cancerDrivers,splitMarkers=TRUE){
  
  #need to separate markers by commas as well for same profile
  GLOBAL<-PriorityResults
  
  #Separate markers so each row is one target and one marker:
  
  if(splitMarkers){
  GLOBAL<-get_SplitMarkers(GLOBAL)}

  GLOBAL$markerName<-unlist(sapply(GLOBAL[,markerColumn],function(x) strsplit(x,"_",fixed=T)[[1]][1]))
  GLOBAL$ID<-paste0(GLOBAL$TARGET,GLOBAL$ctype)
  GLOBAL$SLgroup<-""
  
  #Flag dual markers - do not assign to DFA group
  sel<-grep("&",GLOBAL[,markerColumn])
  GLOBAL[sel,"SLgroup"]<-"Dual"
  
  #Split any markers with same binary pattern. Now that have marked dual markers
  if(splitMarkers){
  GLOBALSplit<-get_SplitMarkers(GLOBAL[GLOBAL$SLgroup!="Dual",],splitVal=",")
  GLOBAL<-rbind(GLOBALSplit,GLOBAL[GLOBAL$SLgroup=="Dual",])}
  GLOBAL$markerName<-unlist(sapply(GLOBAL[,markerColumn],function(x) strsplit(x,"_",fixed=T)[[1]][1]))
  

  #Self oncogenic addiction:
  TT1<-c("SelfAddiction_Mutation","SelfAddiction_Expression","SelfAddiction_CN","SelfAddiction_Variant","SelfAddiction_Protein")
  TargetTypeMatrix<-matrix(0,nrow=nrow(GLOBAL),ncol=5,dimnames=list(GLOBAL$ID,TT1))
  
  TargetLoc<-which(GLOBAL$PPI_min==0)
  checkfor<-c("_mut","_Expr","_CN","_var","_Prot")
  for(i in TargetLoc){
    
    markers<-unlist(GLOBAL[i,markerColumn])
    um<-unlist(sapply(unlist(strsplit(markers,"//",fixed=TRUE)),function(x) strsplit(x,",",fixed=TRUE)))
    if(markerColumn=="NetworkMARKER"){
      bmtype<-unlist(strsplit(GLOBAL[,"BMType"],"//",fixed=TRUE))
      um<-sapply(1:length(um),function(x) paste(um[x],bmtype[x],collapse="_"))
      um<-gsub("_variant","_var",um)
    }
    um<-gsub(" ","",um)
    um<-gsub("_..._var","_var",um)
    um<-gsub("_...._var","_var",um)
    #instead of this way do split by "_" and use bmtype to assign annotations/columns in TargetTypeMatrix
    if(length(grep("_var",markers))>0){
      um<-paste0(GLOBAL[i,"TARGET"],"_var")
    }
    checknames<-paste0(GLOBAL[i,"TARGET"],checkfor)
    temp<-sapply(checknames,function(x) is.element(x,um))
    TargetTypeMatrix[i,]<-temp+0
    GLOBAL[i,"SLgroup"]<-paste(TT1[temp],collapse=",")
    
    
  }

  #Paralogs:
  #need to split markers by // otherwise this wont work. 
  GLOBALP<-GLOBAL
  GLOBALP$paralogF<-sapply(GLOBALP$TARGET,function(x) which(unlist(lapply(paraloglist,function(y) x%in%y))))
  #GLOBAL$ms<-unlist(sapply(GLOBAL[,"MARKER"],function (y) strsplit(y,"_")[[1]][1]))
  GLOBALP$paralogM<-unlist(apply(GLOBALP,1,function(x) length(intersect(c(unlist(sapply(unlist(strsplit(unlist(x[markerColumn]),'//',fixed=T)),function (y) strsplit(y,"_"))),x["TARGET"]),unlist(paraloglist[unlist(x["paralogF"])])))>1))
  GLOBALP$paralogs<-unlist(apply(GLOBALP,1,function(x) paste(intersect(setdiff(unlist(sapply(unlist(strsplit(unlist(x[markerColumn]),'//',fixed=T)),function (y) strsplit(y,"_")[[1]])),x["TARGET"]),unlist(plist[unlist(x["paralogF"])])),collapse=",")))
  depSign<-unlist(GLOBALP[,"ASSOCIATION_EFFECT"])
  depSign<-gsub(" ","",depSign)
  #Decreased Dependency with gene expression or protein intensity marker:
  ParalogID<-which(GLOBALP$paralogM=="TRUE"&depSign=="DecreasedDep.")
  ExprCheck<-grep("_Expr|_Prot|_CN",GLOBALP[,markerColumn])
  ParalogID<-intersect(ParalogID,ExprCheck)
  #Add any Increased dependency on a LoF mutation or variant:
  ParalogID2<-which(GLOBALP$paralogM=="TRUE"&depSign=="IncreasedDep.")
  
  type<-sapply(unlist(GLOBALP[,markerColumn]),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
  type<-gsub(" ","",type)
  type<-type%in%c("mut","var")
  mLoFd<-GLOBALP[,"markerName"]%in%LoFgenes
  LoFCheck<-intersect(which(mLoFd&type),ParalogID2)
  ParalogID<-union(ParalogID,LoFCheck)
  Paralog<-rep(0,nrow(GLOBAL))
  Paralog[ParalogID]<-1
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Paralog)
  GLOBAL[ParalogID,"SLgroup"]<-sapply(ParalogID,function(x) 
    ifelse(GLOBAL[x,"SLgroup"]!="",paste(GLOBAL[x,"SLgroup"],"Paralog",sep=","),"Paralog"))
  
  
  sel<-which(GLOBAL$SLgroup=="")
  
  #Synthetic lethal with Tumour suppressor loss.
  LoFSLmut<-rep(0,nrow(GLOBAL))
  LoFSLother<-rep(0,nrow(GLOBAL))
  LoFActother<-rep(0,nrow(GLOBAL))
  
  for(i in sel){
    markers<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][1])
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    type<-gsub(" ","",type)
    type<-type%in%c("mut","var")
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    mLoFd<-markers%in%LoFgenes
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="IncreasedDep."
    if(is.na(type)){
      check=FALSE
    }else{
      check<-mLoFd&&dep&&type}
    
    if(sum(check)>0){
      LoFSLmut[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"LoF_mutation_SL",sep=","),"LoF_mutation_SL")
    }
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    type<-gsub(" ","",type)
    type<-type%in%c("Expr","CN","Prot")
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="DecreasedDep."
    if(is.na(type)){
      check<-FALSE
    }else{
      check<-mLoFd&&dep&&type}
    
    if(sum(check)>0){
      LoFSLother[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"LoF_other_SL",sep=","),"LoF_other_SL")
    }
    
    dep<-depSign=="IncreasedDep."
    if(is.na(type)){
      check<-FALSE
    }else{
      check<-mLoFd&&dep&&type}
    
    if(sum(check)>0){
      #LoFActother[i]<-1
      #GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"LoF_other_Addiction",sep=","),"LoF_other_Addiction")
    }
    
  }
  #TargetTypeMatrix<-cbind(TargetTypeMatrix,LoFSLmut,LoFSLother,LoFActother)
  TargetTypeMatrix<-cbind(TargetTypeMatrix,LoFSLmut,LoFSLother)
  #MSI
  MSI<-rep(0,nrow(GLOBAL))
  MSISub<-grep("MSI_Composite",GLOBAL[,markerColumn],ignore.case = T)
  MSI[MSISub]<-1
  TargetTypeMatrix<-cbind(TargetTypeMatrix,MSI)
  GLOBAL[MSISub,"SLgroup"]<-sapply(MSISub,function(x) ifelse(GLOBAL[x,"SLgroup"]=="","MSI",GLOBAL[x,"SLgroup"]))
  
  #Transcriptional subtypes
  TSub<-grep("Tsubtype",GLOBAL[,markerColumn])
  CMS<-grep("CMS",GLOBAL[,markerColumn])
  Breastsub<-grep("ScMod",GLOBAL[,markerColumn])
  Breastsub2<-grep('PAM50',GLOBAL[,markerColumn])
  celligner<-grep('celligner',GLOBAL[,markerColumn])
  composite<-grep('_Composite',GLOBAL[,markerColumn])
  MSI<-grep("MSI_Composite",GLOBAL[,markerColumn])
  composite<-setdiff(composite,MSI)
  TSub<-unique(c(TSub,CMS,Breastsub,Breastsub2,celligner,composite))
  Tsubtype<-rep(0,nrow(GLOBAL))
  Tsubtype[TSub]<-1
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Tsubtype)
  GLOBAL[TSub,"SLgroup"]<-sapply(TSub,function(x) ifelse(GLOBAL[x,"SLgroup"]=="","Transcriptional_Subtype",GLOBAL[x,"SLgroup"]))
  
  
  
  
  #Oncogenic addiction
  cDriver<-cbind(cancerDrivers[[1]],cancerDrivers[[2]])
  cDriverAct<-cDriver[cDriver[,2]=="Act",]
  
  cDriverAmb<-cDriver[cDriver[,2]=="ambiguous",]
  MarkerList<-sapply(GLOBAL[,markerColumn],function(x) sapply(strsplit(x,"//",fixed=TRUE),function(y) unlist(strsplit(y,"_",fixed=TRUE))))
  MarkerList<-lapply(MarkerList,function(x) ifelse(is.matrix(x),x[,1],x))
  
  #oncogenic activating
  
  sel<-which(GLOBAL$SLgroup=="")
  
  Target_cancerDriverAct<-rep(0,nrow(GLOBAL))
  for(i in sel){
    marker<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][1])
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    mAct<-marker%in%cDriverAct
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="IncreasedDep."
    check<-mAct&&dep
    
    if(sum(check)>0){
      Target_cancerDriverAct[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"OncogenicAddiction_Act",sep=","),"OncogenicAddiction_Act")
    }
  }
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Target_cancerDriverAct)
  
  sel<-which(GLOBAL$SLgroup=="")
  
  Target_cancerDriverAmb<-rep(0,nrow(GLOBAL))
  for(i in sel){
    marker<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][1])
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    mAct<-marker%in%cDriverAmb
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="IncreasedDep."
    check<-mAct&&dep
    
    if(sum(check)>0){
      #Target_cancerDriverAmb[i]<-1
      #GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"OncogenicAddiction_Amb",sep=","),"OncogenicAddiction_Amb")
    }
  }
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Target_cancerDriverAmb)
  sel<-which(GLOBAL$SLgroup=="")
  
  Target_MetabolicM<-rep(0,nrow(GLOBAL))
  
  for(i in sel){
    marker<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][1])
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    type<-gsub(" ","",type)
    type<-type%in%c("Met","met")
    
    check<-type
    if(is.na(check)){check<-FALSE}
    
    if(sum(check)>0){
      #Target_MetabolicM[i]<-1
      #GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"Metabolite",sep=","),"Metabolite")
    }
    
  }
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Target_MetabolicM)
  #activating addiction other
  sel<-which(GLOBAL$SLgroup=="")
  
  Target_Activating<-rep(0,nrow(GLOBAL))
  
  for(i in sel){
    
    marker<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][1])
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    type<-gsub(" ","",type)
    type<-type%in%c("Expr","Prot","CN")
    mAmb<-!(marker%in%c(cancerDrivers[[1]]))
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="IncreasedDep."
    check<-mAmb&&dep
    if(is.na(check)){check<-FALSE}
    if(sum(check)>0){
      Target_Activating[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"Addiction",sep=","),"Addiction")
    }
    
  }
  TargetTypeMatrix<-cbind(TargetTypeMatrix,Target_Activating)
  

  sel<-which(GLOBAL$SLgroup=="")
  
  SL<-rep(0,nrow(GLOBAL))
  #this is okay because we are keeping one dependency-feature pair per line
  for(i in sel){
    
    type<-sapply(unlist(strsplit(unlist(GLOBAL[i,markerColumn]),"//",fixed=T)),function(z) strsplit(z,"_",fixed=T)[[1]][length(strsplit(z,"_",fixed=T)[[1]])])
    depSign<-unlist(strsplit(unlist(GLOBAL[i,"ASSOCIATION_EFFECT"]),"//",fixed=T))
    type<-gsub(" ","",type)
    type<-type%in%c("Expr","CN","Prot")
    depSign<-gsub(" ","",depSign)
    dep<-depSign=="DecreasedDep."
    check<-type&&dep
    
    if(sum(check)>0){
      SL[i]<-1
      GLOBAL[i,"SLgroup"]<-ifelse(GLOBAL[i,"SLgroup"]!="",paste(GLOBAL[i,"SLgroup"],"SL",sep=","),"SL")
    }
  }
  TargetTypeMatrix<-cbind(TargetTypeMatrix,SL)
  

  
  #sel<-which(GLOBAL$SLgroup%in%c("","OncogenicAddiction_Act","OncogenicAddiction_Amb"))
  

  
  #uncollapse duplicate values to avoid double counting:
  Gset<-apply(GLOBAL[,c("markerName","ID","SLgroup")],1,function(x) paste0(x,collapse="-"))
  
  NoDup<-GLOBAL[!Gset%in%Gset[duplicated(Gset)],]
  NoDttm<-TargetTypeMatrix[!Gset%in%Gset[duplicated(Gset)],]
  Dup<-NULL
  Dupttm<-NULL
  for(i in unique(Gset[duplicated(Gset)])){
    g1<-GLOBAL[which(Gset==i),]
    t1<-TargetTypeMatrix[which(Gset==i),]
    newRes<-g1[1,]
    newRes[,markerColumn]<-paste0(g1[,markerColumn],collapse=",")
    Dup<-rbind(newRes,Dup)
    Dupttm<-rbind(t1[1,],Dupttm)
    
  }
  GLOBAL<-rbind(NoDup,Dup)
  TargetTypeMatrix<-rbind(NoDttm,Dupttm)
  
 GLOBAL<-as.data.frame(GLOBAL)
  
  for(i in 1:ncol(GLOBAL)){
    
    GLOBAL[,i]<-unlist(GLOBAL[,i])
  }
  
  return(list(DFAresults=GLOBAL,TargetTypeMatrix=TargetTypeMatrix))
}

get_Description<-function(PriorityResults,markerColumn,genesets,outputname){
  GLOBAL<-PriorityResults
  blanks<-rep("",nrow(GLOBAL))
  GLOBAL<-cbind(GLOBAL,blanks)
  colnames(GLOBAL)[ncol(GLOBAL)]<-outputname
  TargetTypeMatrix<-NULL
  setnames<-names(genesets)
  TargetTypeMatrix<-matrix(0,nrow=nrow(GLOBAL),ncol=length(genesets))
  rownames(TargetTypeMatrix)<-GLOBAL[,"TARGET"]
  colnames(TargetTypeMatrix)<-setnames
  for(i in 1:length(genesets)){
    inclusion<-which(GLOBAL[,"TARGET"]%in%genesets[[i]])
    TargetTypeMatrix[inclusion,i]<-1
    GLOBAL[inclusion,outputname]<-sapply(inclusion,function(x) ifelse(GLOBAL[x,outputname]=="",setnames[i],GLOBAL[x,outputname]))
    
  }
  return(list(Res=GLOBAL,TTM=TargetTypeMatrix))
}

plot_Biomarkers<-function(Target,BM,BMtype,Input1,Input2=NULL,ctype,logFC,annot,BM2=NULL,celllinecol=NULL,pointcol=NULL,outputdir=NULL,shape=NULL){
  #plot specific biomarker/target results for different dependencies (Target) and biomarkers (BM) for different omics (BMtype)
  temp<-CLnameMapping(Input1,logFC,annotation=annot)
  ctype<-make.names(ctype)
  annot$cancer_type<-make.names(annot$cancer_type)
  PlotData<-NULL
  if(!is.null(celllinecol)){
    cellcol<-matrix(celllinecol,nrow=1,byrow=T)
    colnames(cellcol)<-names(celllinecol)
    if(is.null(names(celllinecol))){
      Error("Names missing for colour variable")
    }
    tempcol<-CLnameMapping(Input1,cellcol,annotation=annot)
    celllinecol<-tempcol$refdata
    
  }
  if(ctype[1]!="PANCAN"){
    CL<-make.names(annot[annot$cancer_type%in%ctype,"model_id"])}else{
      CL<-colnames(temp$inputdata)
    }
  if(is.null(Input2)){
    if(!is.null(celllinecol)){
      ccol<-rep(2,length(intersect(CL,colnames(temp$refdata))))
      names(ccol)<-intersect(CL,colnames(temp$refdata))
      ccol[intersect(names(celllinecol),names(ccol))]<-celllinecol[intersect(names(celllinecol),names(ccol))]
      PlotData<-data.frame(Target_fc=temp$refdata[Target,intersect(CL,colnames(temp$refdata))],
                           Feature=temp$inputdata[BM,intersect(CL,colnames(temp$refdata))],col=ccol[intersect(CL,colnames(temp$refdata))],
                           stringsAsFactors = FALSE)
    }else{
      PlotData<-data.frame(Target_fc=temp$refdata[Target,intersect(CL,colnames(temp$refdata))],
                           Feature=temp$inputdata[BM,intersect(CL,colnames(temp$refdata))],
                           stringsAsFactors = FALSE)
    }
    if(BMtype%in%c("Exp","CN","Prot","Met")){
      if(BMtype=="Exp"){markerUnit="log(TPM+1)" }
      if(BMtype=="CN"){markerUnit="log(CN/ploidy)+1"}
      if(BMtype=="Prot"){markerUnit="Protein intensity"}
      if(BMtype=="Met"){markerUnit="Metabolite intensity"}
      if(sum(PlotData$col,na.rm=T)>0){
        PlotData[is.na(PlotData[,"col"]),"col"]<-2
        colpalette<-c("lightblue","violet","orange","pink")
        print(plot(PlotData$Target_fc,PlotData$Feature,xlab=paste0("LogFC Target: ",Target),ylab=paste0(markerUnit," Marker: ",BM),main=ctype,col=colpalette[PlotData[,"col"]+1]))
        pdf(paste(outputdir,Target,ctype,".pdf",sep="-"))
          print(plot(PlotData$Target_fc,PlotData$Feature,xlab=paste0("LogFC Target: ",Target),ylab=paste0(markerUnit, " Marker: ",BM),main=ctype,col=colpalette[PlotData[,"col"]+1],cex=2))
        dev.off()
        png(paste(outputdir,Target,ctype,".png",sep="-"),res=300,pointsize = 12,unit="in",width=5,height=5)
          print(plot(PlotData$Target_fc,PlotData$Feature,xlab=paste0("LogFC Target: ",Target),ylab=paste0(markerUnit, " Marker: ",BM),main=ctype,col=colpalette[PlotData[,"col"]+1],cex=2))
        dev.off()
      }else{
        if(!is.null(pointcol)){
          print(plot(PlotData$Target_fc,PlotData$Feature,xlab=paste0("LogFC Target: ",Target),ylab=paste0(markerUnit," Marker: ",BM),main=ctype,bg=pointcol,col='grey',pch=21))          
          pdf(paste(outputdir,Target,ctype,".pdf",sep="-"))
            print(plot(PlotData$Target_fc,PlotData$Feature,xlab=paste0("LogFC Target: ",Target),ylab=paste0(markerUnit, " Marker: ",BM),main=ctype,bg=pointcol,col='grey',pch=21,cex=2))          
          dev.off()
          png(paste(outputdir,Target,ctype,".png",sep="-"),res=300,pointsize = 12,unit="in",width=5,height=5)
            print(plot(PlotData$Target_fc,PlotData$Feature,xlab=paste0("LogFC Target: ",Target),ylab=paste0(markerUnit, " Marker: ",BM),main=ctype,bg=pointcol,col='grey',pch=21,cex=2))          
          dev.off()
          
        }else{
          print(plot(PlotData$Target_fc,PlotData$Feature,xlab=paste0("LogFC Target: ",Target),ylab=paste0(markerUnit, " Marker: ",BM),main=ctype))
          pdf(paste(outputdir,Target,ctype,".pdf",sep="-"))
          print(plot(PlotData$Target_fc,PlotData$Feature,xlab=paste0("LogFC Target: ",Target),ylab=paste0(markerUnit, " Marker: ",BM),main=ctype,cex=1.2))
          dev.off()
          png(paste(outputdir,Target,ctype,".png",sep="-"),res=300,pointsize = 12,unit="in",width=5,height=5)
          print(plot(PlotData$Target_fc,PlotData$Feature,xlab=paste0("LogFC Target: ",Target),ylab=paste0(markerUnit, " Marker: ",BM),main=ctype,cex=1.2))
          dev.off()
          
        }
      }
      
      
    }
    if(BMtype=="Discrete"){
      Feature<-BM
      if(!is.null(pointcol)){
        bp<-ggplot(data=PlotData,aes(x=Feature,y=Target_fc,group=Feature))+geom_boxplot()+geom_point(aes(col=pointcol))+
          theme_bw()+
          xlab(paste("Feature",Feature, "Absent/Present"))+ylab("Log Fold Change")
      }else{
        bp<-ggplot(data=PlotData,aes(x=Feature,y=Target_fc,group=Feature))+geom_boxplot()+geom_point()+
          theme_bw()+
          xlab(paste("Feature",Feature, "Absent/Present"))+ylab("Log Fold Change")}
      print(bp)
      
      pdf(paste(outputdir,Target,ctype,".pdf",sep="-"),width=4,height=6,useDingbats = FALSE)
      print(bp)
      dev.off()
      png(paste(outputdir,Target,ctype,".png",sep="-"),res=300,pointsize = 12,unit="in",width=5,height=5)
      print(bp)
      dev.off()
    }
  }else{
    Input1<-temp$inputdata
    FCinput<-temp$refdata[Target,intersect(CL,colnames(temp$refdata))]
    temp2<-CLnameMapping(Input2,logFC,annotation=annot)
    Input2<-temp2$inputdata
    CL<-intersect(CL,colnames(Input2))
    CL<-intersect(CL,colnames(Input1))
    if(BMtype=="CompoundME"){
      cols <- rep("darkgray", length(FCinput))
      cols[which(clgroup == 2)] <- "darkgreen"
      cols[which(clgroup==4)]<-"darkred"
      #split once for interaction
      clgroup<-rep(1,length(FCinput))
      names(clgroup)<-CL
      
      
      t1<-Input1[BM,CL]
      t2<-Input2[BM2,CL]
      
      clgroup[t1]<-2
      
      clgroup[t1|t2]<-4
      beeswarm(FCinput ~ clgroup, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                          makeTransparent("darkgreen"),makeTransparent("darkred")), main=paste(mutationMarker,ms,dependency,ctype,sep=":"),pch = 21, col = c("gray",  "darkgreen",makeTransparent("darkred")), cex = 1.5, ylim = range(FCinput), las = 2, 
               labels = c("WT", "Feat1","Either Feature"), xlab = "", ylab = "",axes=F)
      par(new = TRUE)
      
      boxplot(FCinput ~ clgroup, col = NA, ylim = range(FCinput), frame.plot = FALSE, 
              xaxt = "n",  outline = FALSE,tcl=0.5,tck=-0.01)
      legend("bottom",legend=c("WT", "Feat1","Either Feature"),lty=1,col=c("gray","darkgreen",makeTransparent("darkred")))
      
      clgroup<-rep(1,length(FCinput))
      names(clgroup)<-CL
      clgroup[t2]<-2
      clgroup[t1|t2]<-4
      beeswarm(FCinput ~ clgroup, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                          makeTransparent("darkgreen"),makeTransparent("darkred")), main=paste(mutationMarker,ms,dependency,ctype,sep=":"),pch = 21, col = c("gray",  "darkgreen",makeTransparent("darkred")), cex = 1.5, ylim = range(FCinput), las = 2, 
               labels = c("WT", 'Feat2',"Either Feature"), xlab = "", ylab = "",axes=F)
      par(new = TRUE)
      
      boxplot(FCinput ~ clgroup, col = NA, ylim = range(FCinput), frame.plot = FALSE, 
              xaxt = "n",  outline = FALSE,tcl=0.5,tck=-0.01)
      legend("bottom",legend=c("WT", "Feat2","Either Feature"),lty=1,col=c("gray","darkgreen",makeTransparent("darkred")))
    }
    
    if(BMtype=="Compound"){
      clgroup<-rep(1,length(CL))
      #if have mutation:
      
      names(clgroup)<-CL
      clgroup<-clgroup+Input1[BM,CL]
      
      #second split for expression
      clgroup<-clgroup+Input2[BM2,CL]*5
      #beeswarm on four groups instead
      cols <- rep("darkgray", length(FCinput))
      cols[which(clgroup == 1)] <- "darkgreen"
      cols[which(clgroup==2)]<-"darkblue"
      cols[which(clgroup==6)]<-"darkred"
      FCinput<-FCinput[CL]
      if(!is.null(shape)){
        pointvals<-Input1[shape,CL]
        names(pointvals)<-CL
      
        pointvals[pointvals==0]<-21
        
        pointvals[pointvals==1]<-23
        
      }else{
        pointvals<-rep(21,length(CL))
        names(pointvals)<-CL
      }
      beeswarm(FCinput ~ clgroup, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                          makeTransparent("darkgreen"),makeTransparent("darkred"),makeTransparent("pink")), main=paste(BM,BM2,Target,ctype,sep=":"),pwpch = pointvals, col = c("gray",  "darkgreen",makeTransparent("darkred"),makeTransparent("pink")), cex = 1.5, ylim = range(FCinput), las = 2, 
               labels = c("WT", "WTexp+Mut","NoMut+Expr","Mut+Expr"), xlab = "", ylab = "",axes=F)
      par(new = TRUE)
      
      boxplot(FCinput ~ clgroup, col = NA, ylim = range(FCinput), frame.plot = FALSE, 
              xaxt = "n",  outline = FALSE,tcl=0.5,tck=-0.01)
      legend("top",legend=c("WT", "WTexp+Mut","NoMut+Expr","Mut+Expr"),lty=1,col=c("gray","darkgreen",makeTransparent("darkred"),makeTransparent("pink")))
      
      
      pdf(paste(outputdir,Target,ctype,".pdf",sep="-"))
      beeswarm(FCinput ~ clgroup, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                          makeTransparent("darkgreen"),makeTransparent("darkred"),makeTransparent("pink")), main=paste(BM,BM2,Target,ctype,sep=":"),pwpch = pointvals, col = c("gray",  "darkgreen",makeTransparent("darkred"),makeTransparent("pink")), cex = 1.5, ylim = range(FCinput), las = 2, 
               labels = c("WT", "WTexp+Mut","NoMut+Expr","Mut+Expr"), xlab = "", ylab = "",axes=F)
      par(new = TRUE)
      
      boxplot(FCinput ~ clgroup, col = NA, ylim = range(FCinput), frame.plot = FALSE, 
              xaxt = "n",  outline = FALSE,tcl=0.5,tck=-0.01)
      legend("top",legend=c("WT", "WTexp+Mut","NoMut+Expr","Mut+Expr"),lty=1,col=c("gray","darkgreen",makeTransparent("darkred"),makeTransparent("pink")))
      
      dev.off()
      png(paste(outputdir,Target,ctype,".png",sep="-"),res=300,pointsize = 12,unit="in",width=5,height=5)
      beeswarm(FCinput ~ clgroup, corral = "wrap", bg = c(makeTransparent("gray"), 
                                                          makeTransparent("darkgreen"),makeTransparent("darkred"),makeTransparent("pink")), main=paste(BM,BM2,Target,ctype,sep=":"),pch = 21, col = c("gray",  "darkgreen",makeTransparent("darkred"),makeTransparent("pink")), cex = 1.5, ylim = range(FCinput), las = 2, 
               labels = c("WT", "WTexp+Mut","NoMut+Expr","Mut+Expr"), xlab = "", ylab = "",axes=F)
      par(new = TRUE)
      
      boxplot(FCinput ~ clgroup, col = NA, ylim = range(FCinput), frame.plot = FALSE, 
              xaxt = "n",  outline = FALSE,tcl=0.5,tck=-0.01)
      legend("top",legend=c("WT", "WTexp+Mut","NoMut+Expr","Mut+Expr"),lty=1,col=c("gray","darkgreen",makeTransparent("darkred"),makeTransparent("pink")))
      
      
      dev.off()
    }
  }
  return(PlotData)
}

plotSingleBiomarker<-function(BMres,PPInet,vmap,BMCol='BM',TargetCol="ScoringTargets",TargetNetworkID="STRING_id",BMscoreCol="BMindScores"){
  #From BMPPIpriority output result, extract gene dependencies (Targets) associated with a biomarker and plot
  
  allv_df<-NULL
  BiomarkerTarget<-FALSE
  BM<-BMres[BMCol]
  BM<-strsplit(as.character(BM),'_',fixed=T)[[1]][1]
  Targets<-BMres[TargetCol]
  Targets<-unlist(strsplit(unlist(Targets),"//",fixed=T))
  BMsSID<-PPInet[which(PPInet$symbol==BM),TargetNetworkID]
  TargetIDs<-PPInet[which(PPInet$symbol%in%Targets),TargetNetworkID]
  BMrwr<-unlist(strsplit(unlist(BMres[BMscoreCol]),"//",fixed=T))
  Sids<-shortest_paths(PPIigraph,from=BMsSID,to=TargetIDs,output="both")
  epaths<-Sids$epath
  edges<-vector("list",length=length(epaths))
  
  if(BMsSID%in%TargetIDs){
    tid<-which(TargetIDs==BMsSID)
    BiomarkerTarget<-TRUE
    TargetIDs<-setdiff(TargetIDs,BMsSID)
    BMTrwr<-BMrwr[tid]
    BMrwr<-BMrwr[-tid]
    epaths<-epaths[-tid]
    
  }
  vertexdf<-data.frame(vnames=c(BMsSID,TargetIDs),size=5,color=c("purple",rep('orange',length(TargetIDs))),stringsAsFactors = FALSE)
  allv_df<-rbind(allv_df,vertexdf)
  vs<-unlist(Sids$vpath)
  interdf<-data.frame(vnames=setdiff(names(unlist(vs)),allv_df$vnames),size=3,color="grey",stringsAsFactors = FALSE)
  for(j in 1:length(epaths)){
    temp<-do.call(rbind,sapply(as_ids(epaths[[j]]),function(x) strsplit(x,"[|]")))
    temp<-data.frame(temp,stringsAsFactors = FALSE)
    temp$weight<-as.numeric(BMrwr[j])
    
    
    edges[[j]]<-temp
  }
  if(BiomarkerTarget){
    temp<-data.frame(X1=BMsSID,X2=BMsSID,weight=100)
    edges[[j+1]]<-temp
  }
  edge_df<-data.frame(do.call(rbind,edges),stringsAsFactors = FALSE)
  allv_df<-rbind(allv_df,interdf)
  allv_df<-unique(allv_df)
  allv_df$label.degree<-0
  selGraph<-graph_from_data_frame(edge_df,directed=FALSE,vertices=allv_df)
  
  #vmap[names(vmap)%in%allv_df[allv_df$color=='grey',]$vnames]<-""
  
  V(selGraph)$label<-vmap[V(selGraph)$name]
  
  E(selGraph)$width<-edge_attr(selGraph,"weight")/100
  V(selGraph)$color<-vertex_attr(selGraph,"color")
  V(selGraph)$label.degree<-vertex_attr(selGraph,"label.degree")
  V(selGraph)$label.dist<-2.5
  V(selGraph)$label.family<-"Helvetica"
  V(selGraph)$label.color<-"black"
  return(list("graph"=selGraph,edge_df=edge_df,vertex_df=allv_df))
}

priorityPlot2<-function(ctype,PAB,th,sigOnly=TRUE,IDENTIFY=FALSE,priorSet=NULL,CEX=2,inputPriority,compound=FALSE,
                        TissueTypeColors,PAI=NULL,Pvectors=NULL,allMarkers,priorityres,ymax=90,nlabels=15){
  
  
  TissueColors<-TissueTypeColors
  names(TissueColors)<-make.names(names(TissueColors))
  if(compound){
    load(paste0(inputPriority,'/35_PrioritizedHits/Compound/_PRIORITY_vectors.RData'))
    load(paste0(inputPriority,'/35_PrioritizedHits/Compound/_PRIORITISED_across_indications.RData')) 
  }else{
    if(is.null(PAI)){
      if(ctype!='PANCAN'){
        load(paste0(inputPriority,'/35_PrioritizedHits/_PRIORITY_vectors.RData'))
        load(paste0(inputPriority,'/35_PrioritizedHits/_PRIORITISED_across_indications.RData')) 
      }else{
        load(paste0(inputPriority,'/35_PrioritizedHits/_PRIORITY_vectors_PANCAN.RData'))
        load(paste0(inputPriority,'/35_PrioritizedHits/_PRIORITISED_across_indications_PANCAN.RData'))
      }
    }else{
      PRIORITISED_across_indications<-PAI
      PRIORITY_vectors<-Pvectors
    }
    
  }
  
  TissueColors_tr<-makeTransparent(TissueColors,100)
  
  if(ctype!='PANCAN'){
    indications<-PRIORITISED_across_indications[[ctype]]
  }else{
    indications<-PRIORITISED_across_indications
  }
  
  PCHSYM<-c(23,24,22,21)
  names(PCHSYM)<-names(indications)
  
  
  if(ctype!='PANCAN'){
    TRACTABLE<-PAB[[ctype]]    
  }else{
    TRACTABLE<-PAB
  }
  
  if(sigOnly){
    YL<-c(th-1,ymax)    
  }else{
    YL<-c(1,ymax)
  }   
  if(compound){
    plotmain=allMarkers[allMarkers$ANALYSIS==ctype,"FEATURE_parent"][1]
  }else{
    plotmain<-''}
  
  plot(0,0,xlim=c(0,11),ylim=YL,col=NA,frame.plot=FALSE,xaxt='n',
       main=plotmain,
       xlab='tractability bucket',ylab='target priority')
  
  XX<-NULL
  YY<-NULL
  NN<-NULL
  
  
  currentMarkers<-allMarkers[which(allMarkers$ANALYSIS==ctype),]    
  
  tmp<-lapply(TRACTABLE,function(x){
    ht<-names(which(x>=th))
    scores<-x[which(x>=th)]
  })
  res<-do.call('rbind',lapply(1:10,function(i){cbind(rep(i,length(tmp[[i]])),tmp[[i]])}))
  
  colnames(res)<-c('bucket','priority')
  res<-res[order(res[,'priority'],decreasing=TRUE),]
  print(head(res))
  if(is.null(priorSet)){
    priorSet<-rownames(res)[1:nlabels]
  }
  
  uu<-1:10
  
  
  for (j in 1:length(uu)){
    
    if(sigOnly){
      genes<-names(which(TRACTABLE[[j]]>=th))
    }else{
      genes<-names(TRACTABLE[[j]])
    }
    if(ctype=="PANCAN"){
      symbols<-PCHSYM[colSums(do.call(rbind,lapply(indications,function(x){is.element(genes,names(x))}))*1:3)]
      
    }else{
      symbols<-PCHSYM[colSums(do.call(rbind,lapply(indications,function(x){is.element(genes,names(x))}))*1:4)]
    }
    ngenes<-length(genes)
    
    if(ngenes>0){
      
      if(ngenes==1){
        xc<-(uu[j])
      }else{
        xc<-seq(uu[j]-0.4,uu[j]+0.4,0.8/(ngenes-1))    
      }
      xc<-jitter(xc)
      
      if(ctype=='PANCAN'){
        bgCol<-rgb(120,120,120,alpha = 80,maxColorValue = 255)
        bCol<-rgb(120,120,120,alpha = 255,maxColorValue = 255)
        bgCol<-TissueColors_tr[ctype]
        #bCol<-TissueColors[ctype]
        bCol<-"grey"
      }else{
        if(compound){
          inctype<-ctype
          inctype<-gsub("\\d","",inctype)
          bgCol<-TissueColors_tr[inctype]
          bCol<-TissueColors[inctype]
          
        }else{
          bgCol<-TissueColors_tr[ctype]
          #bCol<-TissueColors[ctype]
          bCol<-"grey"
        }
      }
      
      
      if(ctype!='PANCAN'){
        PP<-PRIORITY_vectors[[ctype]][genes]    
      }else{
        PP<-PRIORITY_vectors[genes]
      }
      PP<-jitter(PP)
      
      bgCol<-rep(bgCol,length(genes))
      names(bgCol)<-genes
      
      points(xc,sort(PP,decreasing=TRUE),bg=bgCol,col=bCol,
             pch=symbols[order(PP,decreasing=TRUE)],cex=CEX+0.5)
      
      ad<-names(sort(PP,decreasing=TRUE))
      if(ctype=="PANCAN"){
        indi<-do.call(rbind,lapply(indications[1:3],function(x){is.element(ad,names(x))}))+0
      }else{
        indi<-do.call(rbind,lapply(indications[1:4],function(x){is.element(ad,names(x))}))+0}
      colnames(indi)<-genes
      
      genes<-names(sort(PP,decreasing=TRUE))
      
      #if(ctype!='PANCAN'){
      MARKERclass<-NULL
      for (k in 1:length(genes)){
        
        tmpM<-sort(unlist(currentMarkers[which(currentMarkers$Depleted.Gene==genes[k]),'CLASS']))[1]
        if(length(tmpM)==0){tmpM<-'N/A'}
        MARKERclass<-c(MARKERclass,tmpM)
      }
      
      newSymbols<-c(8,3,4,NA,NA,NA)
      names(newSymbols)<-c('A','B','C','D','N/A')
      MARKERpoints<-newSymbols[MARKERclass]
      
      
      
      points(xc,sort(PP,decreasing=TRUE),col='black',
             pch=MARKERpoints,cex=max(CEX-2,0.75))
      
      XX<-c(XX,xc)
      YY<-c(YY,sort(PP,decreasing=TRUE))
      NN<-c(NN,names(sort(PP,decreasing=TRUE)))
    }
  }
  
  abline(h=th)
  
  
  abline(v=seq(0.5,10.5,1),col=makeTransparent('black'),lty=1)
  
  axis(side = 1,1:10)
  
  if(length(priorSet)){
    
    ii<-match(priorSet,NN)
    resSel<-priorityres[priorityres$TARGET%in%priorSet&priorityres$ctype==ctype,]
    bestMarkers<-apply(resSel,1,function(x) bestMarkerPPI(x["MARKER"],x["PPI_min"],x["PPI_distance"]))
    names(bestMarkers)<-resSel$TARGET
    #text(XX[ii]+1,YY[ii]+0.5,paste(NN[ii],bestMarkers[NN[ii]],sep=';'),cex=0.5)
    text(XX[ii]+0.4,YY[ii]+0.5,NN[ii])
  }
  
  if(IDENTIFY){
    identify(XX,YY,NN)
  }
  
  
  
  return(res)
}

bestMarkerPPI<-function(markers,ppi_min,ppi_dist){
  mlist<-unlist(strsplit(markers,"//",fixed=TRUE))
  ppilist<-unlist(strsplit(ppi_dist,"//",fixed=TRUE))
  if(!ppi_min%in%c("Inf","N/A")){
    suppressWarnings(muse<-mlist[which(as.numeric(ppilist)<=as.numeric(ppi_min))])
  }else{
    muse<-mlist
  }
  msi<-grep("MSI_Composite",mlist)
  if(length(msi)>0){muse<-c(muse,"MSI_Composite")}
  comp<-grep("Tsubtype",mlist)
  if(length(comp)>0){
    muse<-c(muse,mlist[comp])
  }
  comp<-grep("Celligner",mlist)
  if(length(comp)>0){
    muse<-c(muse,mlist[comp])
  }
  comp<-grep("Progeny",mlist)
  if(length(comp)>0){
    muse<-c(muse,mlist[comp])
  }
  comp<-grep("_met",mlist)
  if(length(comp)>0){
    muse<-c(muse,mlist[comp])
  }
  muse<-unique(muse)
  return(paste(muse,collapse="//"))
}
makeTransparent<-function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
superPriorityPlot<-function(TOTRES,allMarkers,plotname,TissueColors,textThresh=50,ylims=c(30,90),shape=NULL,indi=NULL,pointsize=2,plotsuffix=NULL){
  if(!is.null(shape)){
    TOTRES$shape<-23
    TOTRES[TOTRES$MARKERCLASS%in%c("A","B","C","D"),"shape"]<-21
    
  }else{
    TOTRES$shape<-21
  }
  TOTRES$symbol<-NA
  if(!is.null(indi)){
    
    TOTRES[TOTRES$RWRscore=="100","symbol"]<-8
    #TOTRES[TOTRES$RWRscore=="75","symbol"]<-4
  }
  if(is.null(plotsuffix)){
    pdf(paste0(plotname,'SuperPriorityPlot.pdf'),20,30,useDingbats = FALSE)
  }else{
    pdf(paste0(plotname,plotsuffix),20,30,useDingbats = FALSE)
  }
  layout(t(matrix(c(rep(1,6),rep(2:4,5)),3,7)))
  uu<-unique(TOTRES$BUCKET)
  par(mar=c(4,4,4,0))
  plot(0,0,xlim=c(min(uu)-1,max(uu)+1),ylim=ylims,col=NA,frame.plot=FALSE,xaxt='n',cex.axis=3,
       main='',
       xlab='tractability bucket',ylab='')
  TissueColors_tr<-makeTransparent(TissueColors,100)
  
  
  
  ALLX<-NULL
  ALLY<-NULL
  ALLLAB<-NULL
  for (i in 1:length(uu)){
    
    id<-which(TOTRES$BUCKET==uu[i] & TOTRES$ctype!='PANCAN')
    id<-which(TOTRES$BUCKET==uu[i])
    id<-id[order(TOTRES$PRIORITYL3[id],decreasing = TRUE)]
    
    genes<-TOTRES$TARGET[id]
    
    ngenes<-length(genes)
    
    if(ngenes==1){
      xc<-(uu[i])
    }else{
      xc<-seq(uu[i]-0.4,uu[i]+0.4,0.8/(ngenes-1))
    }
    
    points(xc,TOTRES$PRIORITYL3[id],pch=TOTRES$shape[id],cex=4,
           bg=TissueColors_tr[TOTRES$ctype[id]],col=TissueColors[TOTRES$ctype[id]])
    if(sum(!is.na(TOTRES$symbol[id]))>0){
      idcheck<-id[which(!is.na(TOTRES$symbol[id]))]
      points(xc[which(!is.na(TOTRES$symbol[id]))],TOTRES$PRIORITYL3[idcheck],pch=TOTRES$symbol[idcheck],cex=2)
    }
    #text(xc,TOTRES$PRIORITY[id],genes,pos = 4,cex=0.2)
    
    ALLX<-c(ALLX,xc)
    ALLY<-c(ALLY,TOTRES$PRIORITYL3[id])
    ALLLAB<-c(ALLLAB,genes)
  }
  
  #identify(ALLX,ALLY,ALLLAB)
  
  abline(h=th)
  abline(h=priority_threshold,col='darkgray')
  
  #abline(v=1:10,col=makeTransparent('black'),lty=2)
  abline(v=seq(0.5,10.5,1),col=makeTransparent('black'),lty=1)
  
  axis(side = 1,1:10,cex.axis=3)
  
  
  id<-which(TOTRES$BUCKET<4)
  targets<-TOTRES$TARGET[id]
  targets<-unique(unlist(allMarkers[which(is.element(allMarkers$Depleted.Gene,targets) & is.element(allMarkers$CLASS,c('A','B','C'))),3]))
  currentTOTRES<-TOTRES[id,]
  currentTOTRES<-currentTOTRES[which(is.element(currentTOTRES$TARGET,targets)),]
  currentTOTRES<-currentTOTRES[order(currentTOTRES$PRIORITYL3,decreasing = TRUE),]
  
  
  
  
  id<-which(TOTRES$BUCKET>3 & TOTRES$BUCKET<7)
  targets<-TOTRES$TARGET[id]
  targets<-unique(unlist(allMarkers[which(is.element(allMarkers$Depleted.Gene,targets) & is.element(allMarkers$CLASS,c('A','B','C'))),3]))
  currentTOTRES<-TOTRES[id,]
  currentTOTRES<-currentTOTRES[which(is.element(currentTOTRES$TARGET,targets)),]
  currentTOTRES<-currentTOTRES[order(currentTOTRES$PRIORITYL3,decreasing = TRUE),]
  
  
  
  ALLX<-NULL
  ALLY<-NULL
  ALLLAB<-NULL
  for (i in 1:length(uu)){
    
    id<-which(TOTRES$BUCKET==uu[i] & TOTRES$ctype!='PANCAN')
    id<-which(TOTRES$BUCKET==uu[i])
    id<-id[order(TOTRES$PRIORITYL3[id],decreasing = TRUE)]
    
    genes<-TOTRES$TARGET[id]
    
    ngenes<-length(genes)
    
    if(ngenes==1){
      xc<-(uu[i])
    }else{
      xc<-seq(uu[i]-0.4,uu[i]+0.4,0.8/(ngenes-1))
    }
    
    points(xc,TOTRES$PRIORITYL3[id],pch=TOTRES$shape[id],cex=1,
           bg=TissueColors_tr[TOTRES$ctype[id]],col=TissueColors[TOTRES$ctype[id]])
    if(sum(!is.na(TOTRES$symbol[id]))>0){
      idcheck<-id[which(!is.na(TOTRES$symbol[id]))]
      points(xc[which(!is.na(TOTRES$symbol[id]))],TOTRES$PRIORITYL3[idcheck],pch=TOTRES$symbol[idcheck],cex=1)
    }
    nID<-length(id)
    for(k in 1:nID){
      if(TOTRES$PRIORITYL3[id[k]]>textThresh){
        text(xc[k],TOTRES$PRIORITYL3[id[k]],genes[k],pos = 4,cex=pointsize)
        if(i==1){
          print(TOTRES[id[k],])
        }
      }
    }
    ALLX<-c(ALLX,xc)
    ALLY<-c(ALLY,TOTRES$PRIORITYL3[id])
    ALLLAB<-c(ALLLAB,genes)
  }
  plot(0,0,xlim=c(0,4),ylim=c(1,length(TissueColors)))
  for(i in 1:length(TissueColors)){
    points(1,i,pch=19,col=TissueColors_tr[i],cex=2)
    text(2,i,names(TissueColors)[i])
  }
  dev.off()
}

msplitF<-function(markerrow,featCol="FEATURE",typeCol="MARKER_TYPE",assoCol="ASSOCIATION_EFFECT",IQRcol="IQR",SLcol=NULL,PPIcol="PPI_distance",splitVal="//"){
  #need to account for same feature prevalence e.g. where separated with a comma.
  mkr<-markerrow[featCol]
  ftype<-NULL
  splitSL<-NULL
  splitC1<-unlist(sapply(unlist(mkr),function(x) strsplit(x,splitVal,fixed=T)))
  splitC<-unlist(sapply(splitC1,function(x) strsplit(x,",",fixed=T)))
  nvals<-unlist(lapply(sapply(splitC1,function(x) strsplit(x,",",fixed=T)),length))
  CoOccID<-grep(",",splitC1)
  if(!is.null(typeCol)){
    splitM<-unlist(sapply(splitC,function(x) ifelse(length(strsplit(x,"_",fixed=T)[[1]])==2,strsplit(x,"_",fixed=T)[[1]][2],ifelse(strsplit(x,"_",fixed=T)[[1]][3]%in%c("var","Composite"),strsplit(x,"_",fixed=T)[[1]][3],"dual"))))}
  
  if(!is.null(SLcol)){
    splitSL<-unlist(sapply(unlist(markerrow[SLcol]),function(x) strsplit(x,",",fixed=T)))
  }
  nmarkers<-length(splitC)
  markers<-splitC
  
  markers<-gsub(" ","",markers, fixed = T)
  if(!is.null(assoCol)){
    if(length(splitC1)==length(splitC)){
      assocC<-unlist(sapply(unlist(markerrow[assoCol]),function(x) strsplit(x,splitVal,fixed=T)))
    }else{
      temp<-unlist(sapply(unlist(markerrow[assoCol]),function(x) strsplit(x,splitVal,fixed=T)))
      assocC<-c()
      for(k in 1:length(temp)){
        if(k%in%CoOccID){
          assocC<-c(assocC,rep(temp[k],nvals[k]))
        }else{
          assocC<-c(assocC,temp[k])
        }
      }
    }
    
  }
  if(!is.null(PPIcol)){
    if(length(splitC1)==length(splitC)){
      PPIval<-unlist(sapply(unlist(markerrow[PPIcol]),function(x) strsplit(x,splitVal,fixed=T)))
    }else{
      temp<-unlist(sapply(unlist(markerrow[PPIcol]),function(x) strsplit(x,splitVal,fixed=T)))
      PPIval<-c()
      for(k in 1:length(temp)){
        if(k%in%CoOccID){
          PPIval<-c(PPIval,rep(temp[k],nvals[k]))
        }else{
          PPIval<-c(PPIval,temp[k])
        }
      }
    }
    
  }
  if(IQRcol%in%names(markerrow)){
    if(length(splitC1)==length(splitC)){
      IQR<-unlist(sapply(unlist(markerrow[IQRcol]),function(x) strsplit(x,splitVal,fixed=T)))
    }else{
      temp<-unlist(sapply(unlist(markerrow[IQRcol]),function(x) strsplit(x,splitVal,fixed=T)))
      IQR<-c()
      for(k in 1:length(temp)){
        if(k%in%CoOccID){
          IQR<-c(IQR,rep(temp[k],nvals[k]))
        }else{
          IQR<-c(IQR,temp[k])
        }
      }
    }
    
  }
  
  if(nmarkers>1){
    newRes<-matrix(rep(markerrow,nmarkers),byrow = T,ncol=length(markerrow))
    colnames(newRes)<-names(markerrow)
    newRes[,featCol]<-splitC
    if(!is.null(typeCol)){
      
      newRes[,typeCol]<-splitM
      
    }
    if(!is.null(assoCol)){
      newRes[,assoCol]<-assocC
    }
    if(!is.null(PPIcol)){
      newRes[,PPIcol]<-PPIval
    }

    if(!is.null(SLcol)&(length(splitSL)==nrow(newRes))){
      newRes[,SLcol]<-splitSL
    }
    if(IQRcol%in%names(markerrow)){
      newRes[,IQRcol]<-IQR  
      
    }
  }else{
    
    newRes<-matrix(markerrow,nrow=1)
    
    
    
  }
  return(newRes)
}

Get_PrevalenceMatrix<-function(PriorityTargets,MutMat,CNgain,CNloss,ExprMat=NULL,ExprUp=NULL,ExprDown=NULL,Mcol="MARKER_TYPE",Score="PRIORITYL3"){
  PT<-PriorityTargets[order(PriorityTargets[,Score],decreasing=T),]
  PT$M<-unlist(sapply(PT$MARKER,function(x) gsub(" ","",strsplit(x,"_",fixed=T)[[1]][1])))
  if(!is.null(ExprMat)){
    Psamples<-intersect(colnames(MutMat),intersect(intersect(colnames(CNgain),colnames(CNloss)),colnames(ExprMat)))
  }else{
    if(!is.null(ExprUp)){
      Psamples<-intersect(colnames(MutMat),intersect(intersect(colnames(CNgain),colnames(CNloss)),colnames(ExprUp)))
    }else{
      Psamples<-intersect(colnames(MutMat),intersect(colnames(CNgain),colnames(CNloss)))
      
    }
  }
  if(!is.null(MutMat)){
  MutMat<-MutMat[,Psamples]}
  if(!is.null(ExprMat)){
    Emat<-dimnames(ExprMat)
    ExprMat<-ExprMat[,Psamples]
    if(!is.matrix(ExprMat)){
      ExprMat<-matrix(ExprMat,nrow=1,ncol=length(Psamples))
      rownames(ExprMat)<-Emat[[1]]
      colnames(ExprMat)<-Psamples
    }
  }
  PriorityMatrix<-matrix(0,nrow=nrow(PT),ncol=length(Psamples))
  rownames(PriorityMatrix)<-PT$ID
  colnames(PriorityMatrix)<-Psamples
  for(i in 1:nrow(PT)){
    mtype<-PT[i,Mcol]
    if(mtype%in%c("mut")){
      #check direction of dependency
      if(PT[i,"M"]%in%rownames(MutMat)){
        if(PT[i,"ASSOCIATION_EFFECT"]%in%c("Increased Dep.","IncreasedDep.")){
          
          sel<-MutMat[PT[i,"M"],]
          PriorityMatrix[i,]<-sel[Psamples]
        }else{
          sel<-MutMat[PT[i,"M"],]
          PriorityMatrix[i,]<-(sel[Psamples]==0)+0
        }
      }else{
        PriorityMatrix[i,]<-0
      }
    }
    if(mtype%in%c("cna","cn")){
      #check direction of dependency
      if(PT[i,"M"]%in%rownames(CNloss)){
        if(PT[i,"ASSOCIATION_EFFECT"]%in%c("Decreased Dep.","DecreasedDep.")){
          sel<-CNloss[PT[i,"M"],]
          PriorityMatrix[i,]<-sel[Psamples]
        }else{
          sel<-CNgain[PT[i,"M"],]
          PriorityMatrix[i,]<-sel[Psamples]
        }
      }else{
        PriorityMatrix[i,]<-0
      }
    }
    if(mtype%in%c("expr")){
      if(!is.null(ExprMat)){
        gID<-paste(PT[i,"TARGET"],PT[i,"M"],sep="-")
        sel<-ExprMat[gID,]
        PriorityMatrix[i,]<-sel[Psamples]
      }else{
        if(PT[i,"M"]%in%rownames(ExprUp)){
          if(PT[i,"ASSOCIATION_EFFECT"]%in%c("Increased Dep.","IncreasedDep.")){
            
            sel<-ExprUp[PT[i,"M"],]
            PriorityMatrix[i,]<-sel[Psamples]
          }else{
            sel<-ExprDown[PT[i,"M"],]
            PriorityMatrix[i,]<-sel[Psamples]
          }
        }else{
          PriorityMatrix[i,]<-0
        }
      }
    }
  }
  
  return(PriorityMatrix)
}

Plot_PrevalenceHM<-function(PrevMat,PriorityRes,plotCol,outputname=NULL,DFA=NULL,PSorder=FALSE){
  inputMat<-PrevMat[rowSums(PrevMat,na.rm=T)>0,]
  inputMat<-inputMat[order(rowSums(inputMat),decreasing=T),]
  PRes<-PriorityRes[PriorityRes$ID%in%rownames(inputMat),]
  ctype<-PRes$ctype
  mat_list<-vector("list",length=3)
  Einput<-NULL
  Cinput<-NULL
  Minput<-NULL
  nTargets<-unique(paste(PRes$TARGET,PRes$M,sep=":"))
  
  #Expression markers:
  Emat<-matrix(0,nrow=length(nTargets),ncol=ncol(inputMat))
  rownames(Emat)<-nTargets
  colnames(Emat)<-colnames(inputMat)
  emarkers<-paste0(PRes$TARGET,PRes$ctype,PRes$M,"expr")
  #patient expression marker subset, Einput:
  usem<-intersect(emarkers,rownames(inputMat))
  if(length(usem)>1){
    Einput<-inputMat[usem,]
  }
  if(length(usem)==1){
    Einput<-matrix(inputMat[usem,],nrow=1,ncol=ncol(inputMat),byrow=TRUE)
    rownames(Einput)<-usem
    colnames(Einput)<-colnames(inputMat)
  }
  if(!is.null(Einput)){
    #targets relating to the expression markers:
    inT<-unlist(sapply(rownames(Einput),function(x) strsplit(x,ctype,fixed=T)[[1]][1]))
    exprMarker<-rownames(Einput)
    exprMarker<-unlist(sapply(exprMarker,function(x) strsplit(x,ctype,fixed=T)[[1]][2]))
    exprMarker<-unlist(sapply(exprMarker,function(x) strsplit(x,"expr",fixed=T)[[1]][1]))
    Emat[paste(inT,exprMarker,sep=":"),]<-Einput
  }
  #Copy number markers:
  CNmat<-matrix(0,nrow=length(nTargets),ncol=ncol(inputMat))
  rownames(CNmat)<-nTargets
  colnames(CNmat)<-colnames(inputMat)
  cnmarkers<-paste0(PRes$TARGET,PRes$ctype,PRes$M,"cn")
  
  usem<-intersect(cnmarkers,rownames(inputMat))
  if(length(usem)>1){
    Cinput<-inputMat[usem,]
  }
  if(length(usem)==1){
    Cinput<-matrix(inputMat[usem,],nrow=1,ncol=ncol(inputMat),byrow=TRUE)
    rownames(Cinput)<-usem
    colnames(Cinput)<-colnames(inputMat)
  }
  if(!is.null(Cinput)){
    inT<-unlist(sapply(rownames(Cinput),function(x) strsplit(x,ctype,fixed=T)[[1]][1]))
    cnMarker<-rownames(Cinput)
    cnMarker<-unlist(sapply(cnMarker,function(x) strsplit(x,ctype,fixed=T)[[1]][2]))
    cnMarker<-unlist(sapply(cnMarker,function(x) strsplit(x,"cn",fixed=T)[[1]][1]))
    CNmat[paste(inT,cnMarker,sep=":"),]<-Cinput
  }
  #Mutation markers:
  Mutmat<-matrix(0,nrow=length(nTargets),ncol=ncol(inputMat))
  rownames(Mutmat)<-nTargets
  colnames(Mutmat)<-colnames(inputMat)
  mmarkers<-paste0(PRes$TARGET,PRes$ctype,PRes$M,"mut")
  
  usem<-intersect(mmarkers,rownames(inputMat))
  if(length(usem)>1){
    Minput<-inputMat[usem,]
  }
  if(length(usem)==1){
    Minput<-matrix(inputMat[usem,],nrow=1,ncol=ncol(inputMat),byrow=TRUE)
    rownames(Minput)<-usem
    colnames(Minput)<-colnames(inputMat)
  }
  if(!is.null(Minput)){
    inT<-unlist(sapply(rownames(Minput),function(x) strsplit(x,ctype,fixed=T)[[1]][1]))
    mMarker<-rownames(Minput)
    mMarker<-unlist(sapply(mMarker,function(x) strsplit(x,ctype,fixed=T)[[1]][2]))
    mMarker<-unlist(sapply(mMarker,function(x) strsplit(x,"mut",fixed=T)[[1]][1]))
    Mutmat[paste(inT,mMarker,sep=":"),]<-Minput
  }
  
  
  mat_list = list(Mutmat,Emat,CNmat)
  names(mat_list)<-c("Mutation","Expression","CopyNumber")
  rownames(mat_list$Mutation) = rownames(mat_list$Expression) = rownames(mat_list$CopyNumber) = nTargets
  colnames(mat_list$Mutation) = colnames(mat_list$Expression) = colnames(mat_list$CopyNumber) = colnames(PrevMat)
  
  col_funPS = colorRamp2(c(0, 65, 100), c("yellow", "orange", "red"))
  tractCols = brewer.pal(7,"Blues")
  names(tractCols)<-as.character(c(1:7))
  
  PScore<-PRes[match(nTargets,paste(PRes$TARGET,PRes$M,sep=":")),"PRIORITYL3"]
  names(PScore)<-nTargets
  PScore<-as.numeric(PScore)
  Tgroups<-as.character(PRes[match(nTargets,paste(PRes$TARGET,PRes$M,sep=":")),"GROUP"])
  names(Tgroups)<-nTargets
  if(!is.null(DFA)){
    SLcols<-brewer.pal(length(DFA)+1, "Set3")
    names(SLcols)<-c(DFA,"None")
    PRes[!PRes$SLgroup%in%DFA,"SLgroup"]<-'None'
    SLgroups<-PRes[match(nTargets,paste(PRes$TARGET,PRes$M,sep=":")),"SLgroup"]
    names(SLgroups)<-nTargets
    
    rowAnnot<-rowAnnotation(PS=PScore,TG=Tgroups,DFA=SLgroups,col=list(PS=col_funPS,TG=tractCols,DFA=SLcols))
  }else{
    rowAnnot<-rowAnnotation(PS=PScore,TG=Tgroups,col=list(PS=col_funPS,TG=tractCols))
    
  }
  if(PSorder){
    rorder<-order(PScore,decreasing=T)
  }else{
    rorder<-NULL
  }
  if(!is.null(outputname)){
    pdf(paste0(outputname,".pdf"),useDingbats = FALSE,width=4.5,height=max(round(nrow(inputMat)*0.25),2.5))
      print(oncoPrint(mat_list,
                    left_annotation = rowAnnot,row_order=rorder,
                    alter_fun = list(
                      background = function(...) NULL,remove_empty_rows = FALSE,
                      
                      Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                gp = gpar(fill = plotCol["Mutation"], col = NA)),
                      Expression = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                                  gp = gpar(fill = plotCol["Expression"], col = NA)),
                      CopyNumber = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, 
                                                                  gp = gpar(fill = plotCol["CopyNumber"], col = NA))
                    ), col = plotCol,top_annotation = NULL,right_annotation = NULL)
      )
    dev.off()
    png(paste0(outputname,".png"),width=8,height=4,units="in",res=300,pointsize = 12)
      print(oncoPrint(mat_list,
                    left_annotation = rowAnnot,row_order=rorder,
                    alter_fun = list(
                      background = function(...) NULL,remove_empty_rows = FALSE,
                      
                      Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                gp = gpar(fill = plotCol["Mutation"], col = NA)),
                      Expression = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                                  gp = gpar(fill = plotCol["Expression"], col = NA)),
                      CopyNumber = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, 
                                                                  gp = gpar(fill = plotCol["CopyNumber"], col = NA))
                    ), col = plotCol,top_annotation=NULL)
      )
    dev.off()
  }else{
    print(oncoPrint(mat_list,
                    alter_fun = list(
                      Mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                                gp = gpar(fill = plotCol["Mutation"], col = NA)),
                      Expression = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                                  gp = gpar(fill = plotCol["Expression"], col = NA)),
                      CopyNumber = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2, 
                                                                  gp = gpar(fill = plotCol["CopyNumber"], col = NA))
                    ), col = plotCol,top_annotation=NULL)
    )
  }
}
