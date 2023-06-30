
load(paste0(inputdata,"EXPpcBC.Rdata"))

load(paste0(inputdata,"07_EssMatrix_qnorm_corrected_logFCs.RData"))
load(paste0(inputdata,"07_EssMatrix_bDepletionsB2.Rdata"))
load(paste0(outputdata,"DFAdataPANCAN.Rdata"))
GLOBALP<-DFAdata[[1]]
#want to know how many of the paralog markers are also depleted in the same cancer type. 
ml<-read.csv(paste0(inputdata,"model_list_latest.csv"),header=T,stringsAsFactors=F)
ml2<-ml
ml2$model_id<-make.names(ml2$BROAD_ID)
ml<-rbind(ml,ml2)
ml$cancer_type<-make.names(ml$cancer_type)

temp<-CLnameMapping(EXPpcBC,qnorm_corrected_logFCs,annotation=ml)

EXP<-temp$inputdata
CData<-temp$refdata

temp<-CLnameMapping(EXPpcBC,bDepletionsB2,annotation=ml)
BData<-temp$refdata
#now need to add columns min and max expression of target
#data is log (TPM+1,2) so val less than 1 means less than 1 tpm.
ParalogPriority<-GLOBALP[grep("Paralog",GLOBALP$SLgroup),]
ParalogPriority<-ParalogPriority[ParalogPriority$MARKER_TYPE=="expr",]

ParalogPriority$target_minExp<-sapply(ParalogPriority$TARGET,function(x) min(EXP[x,]))
ParalogPriority$target_m_cor<-apply(ParalogPriority,1,function(x) ifelse(x["markerName"]%in%rownames(EXP),cor(CData[unlist(x["TARGET"]),],EXP[unlist(x["markerName"]),]),NA))
ParalogPriority$marker_t_cor<-apply(ParalogPriority,1,function(x) ifelse(x["markerName"]%in%rownames(CData),cor(CData[unlist(x["markerName"]),],EXP[unlist(x["TARGET"]),]),NA))
ParalogPriority$marker_nDep<-apply(ParalogPriority,1,function(x) ifelse(x["markerName"]%in%rownames(BData),sum(BData[unlist(x["markerName"]),],na.rm=T),NA))


paralogThresh<-quantile(ParalogPriority$target_m_cor,0.25,na.rm=T)
cat(paste("Paralog correlation threshold:",paralogThresh),file=paste0(outputdata,"/ParalogAnalysis.txt"),sep="\n",append=F)

PP<-ParalogPriority[ParalogPriority$marker_nDep>20,]
cat(paste("Number of paralog pairs showing reciprocal buffering:",sum(PP$marker_t_cor>paralogThresh|PP$marker_t_cor>PP$target_m_cor,na.rm=T)
          ),file=paste0(outputdata,"/ParalogAnalysis.txt"),sep="\n",append=T)

cat(paste("Number of paralog pairs no buffering, with partner expressed and marker depleted in at least 20 cell line pancan:",
          sum(PP$marker_t_cor<paralogThresh&PP$target_minExp>1&PP$marker_t_cor<PP$target_m_cor,na.rm=T)
),file=paste0(outputdata,"/ParalogAnalysis.txt"),sep="\n",append=T)

cat(paste("Number of paralog pairs with partner depleted:",sum(ParalogPriority$marker_nDep>20,na.rm=T)
),file=paste0(outputdata,"/ParalogAnalysis.txt"),sep="\n",append=T)



pdata<-data.frame(model=colnames(CData),var1=CData["CDK6",],var2=EXP["CDK6",colnames(CData)],ctype=ml[match(colnames(CData),ml$model_id),"cancer_type"],stringsAsFactors = FALSE)
pdf(paste0(outputdata,"Pancan_CDK6.pdf"),useDingbats = FALSE)
print(ggplot(data=pdata, aes(x=var1,y=var2,col=ctype))+geom_point()+theme(
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "grey"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "grey"),legend.position = "none"))
 dev.off()   
 
 pdata<-data.frame(model=colnames(CData),var1=CData["CDK6",],var2=EXP["CDK4",colnames(CData)],ctype=ml[match(colnames(CData),ml$model_id),"cancer_type"],stringsAsFactors = FALSE)
 pdf(paste0(outputdata,"Pancan_CDK6_CDK4expr.pdf"),useDingbats = FALSE)
 print(ggplot(data=pdata, aes(x=var1,y=var2,col=ctype))+geom_point()+theme(
   panel.background = element_rect(fill = "white",
                                   colour = "white",
                                   size = 0.5, linetype = "solid"),
   panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                   colour = "grey"), 
   panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                   colour = "grey"),legend.position = "none"))
 dev.off()   
 
 pdata<-data.frame(model=colnames(CData),var1=CData["CDK4",],var2=EXP["CDK6",colnames(CData)],ctype=ml[match(colnames(CData),ml$model_id),"cancer_type"],stringsAsFactors = FALSE)
 pdf(paste0(outputdata,"Pancan_CDK4_CDK6expr.pdf"),useDingbats = FALSE)
 print(ggplot(data=pdata, aes(x=var1,y=var2,col=ctype))+geom_point()+theme(
   panel.background = element_rect(fill = "white",
                                   colour = "white",
                                   size = 0.5, linetype = "solid"),
   panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                   colour = "grey"), 
   panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                   colour = "grey"),legend.position = "none"))
 dev.off()    
 



