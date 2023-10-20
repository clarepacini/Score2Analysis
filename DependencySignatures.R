
#################################################################################################
#                                                                                               #
#   Analysis of NMF results using selective dependencies                                        #
#                                                                                               #
#                                                                                               #
#################################################################################################



load(paste0(inputdata,"pcgene.rdata"))
hallmark<-read.gmt(paste0(inputdata,"hallmark.gmt"))
kegg<-read.gmt(paste0(inputdata,"kegg.gmt"))
#using normLRT:
load(paste0(outputdata,"/normLRTresnq.Rdata"))
LRT<-normLRTres$LRT
summary(LRT)
top<-2000
topLRTgenes<-names(sort(LRT,decreasing=T))[1:top]
load(paste0(inputdata,'10_PANCANCER_coreFitness_genes.RData'))
load(paste0(inputdata,'10_PANCANCER_CommonEssentialsAUC.RData'))
load(paste0(inputdata,"07_EssMatrix_bDepletionsB2.Rdata"))
bDepletions<-bDepletionsB2
cmp<-read.csv(paste0(inputdata,'model_list_latest.csv'),header=T,stringsAsFactors = FALSE)

cmp2<-cmp
cmp2$model_id<-cmp2$BROAD_ID
cmp<-rbind(cmp,cmp2)
load(paste0(outputdata,"/inputdataNMFnq.Rdata"))
load(paste0(inputdata,"/nmfRes2_10nq.Rdata"))

pdf(paste0(outputdata,"/cophenetic_DepSig.pdf"),useDingbats = FALSE)
plot(nmfRes2_10$measures$rank,nmfRes2_10$measures$cophenetic,type="l")
dev.off()
pdf(paste0(outputdata,"/silhouette_DepSig.pdf"),useDingbats = FALSE)
plot(nmfRes2_10$measures$rank,nmfRes2_10$measures$silhouette.consensus,type="l")
dev.off()

median(nmfRes2_10$measures$cophenetic)
median(nmfRes2_10$measures$silhouette.consensus)
urank<-3
set.seed(123)
nmfRes<-nmf(inputdataNMF[setdiff(topLRTgenes,union(PanCancerCoreFitnessGenes,CommonEssentialsAUC$cfgenes)),],urank,seed=123)
W<-nmfRes@fit@W
H<-nmfRes@fit@H


softThresh<-median(apply(W,2,function(x) quantile(x,0.7)))

Gsets<-apply(W,2,function(x) rownames(W)[which(x>softThresh)])
enrichResNMF<-lapply(Gsets,function(x) gost(x,user_threshold = 1e-03))

Winc<-W>softThresh+0

KEGGenrichNMF<-lapply(enrichResNMF,function(x) x$result[x$result$source=="KEGG",])
ReactomeEnrichNMF<-lapply(enrichResNMF,function(x) x$result[x$result$source=="REAC",])
GOBPEnrichNMF<-lapply(enrichResNMF,function(x) x$result[x$result$source=="GO:BP",])
rownames(H)<-c("Metabolic","Apoptosis","EMT")
save(H,file=paste0(outputdata,"/H.Rdata"))

#Supplementary Table 3:
for(i in 1:length(KEGGenrichNMF)){
  write.table(as.matrix(KEGGenrichNMF[[i]]),file=paste0(outputdata,"/KEGGenrichSignatureTop2k_",i,".txt"),quote=F,sep='\t',row.names=FALSE)
  write.table(as.matrix(ReactomeEnrichNMF[[i]]),file=paste0(outputdata,"/ReactomeenrichSignatureTop2k_",i,".txt"),quote=F,sep='\t',row.names=FALSE)
  write.table(as.matrix(GOBPEnrichNMF[[i]]),file=paste0(outputdata,"/GOBPenrichSignatureTop2k_",i,".txt"),quote=F,sep='\t',row.names=FALSE)
  
}





