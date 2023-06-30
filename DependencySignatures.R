
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

#Supplementary Table 3:
for(i in 1:length(KEGGenrichNMF)){
  write.table(as.matrix(KEGGenrichNMF[[i]]),file=paste0(outputdata,"/KEGGenrichSignatureTop2k_",i,".txt"),quote=F,sep='\t',row.names=FALSE)
  write.table(as.matrix(ReactomeEnrichNMF[[i]]),file=paste0(outputdata,"/ReactomeenrichSignatureTop2k_",i,".txt"),quote=F,sep='\t',row.names=FALSE)
  write.table(as.matrix(GOBPEnrichNMF[[i]]),file=paste0(outputdata,"/GOBPenrichSignatureTop2k_",i,".txt"),quote=F,sep='\t',row.names=FALSE)
  
}




GOBPfilt<-GOBPEnrichNMF
for(i in 1:length(GOBPEnrichNMF)){
  GOBPfilt[[i]]<-GOBPfilt[[i]][GOBPfilt[[i]][,"recall"]>0.1,]
}

GOBPfilt2<-GOBPfilt
idx<-c(2,3,1)
for(i in 1:length(GOBPfilt)){
  temp<-GOBPfilt[[idx[i]]]
  temp<-temp[order(temp[,"p_value"],decreasing=F),]
  nval<-min(nrow(temp),15)
  GOBPfilt2[[i]]<-temp[1:nval,]
}

allKEGG<-unique(unlist(lapply(KEGGenrichNMF,function(x) x$term_name)))
allREAC<-unique(unlist(lapply(ReactomeEnrichNMF,function(x) x$term_name)))
allGOBP<-unique(unlist(lapply(GOBPfilt,function(x) x$term_name)))
allGOBP2<-unique(unlist(lapply(GOBPfilt2,function(x) x$term_name)))

KEGGmat<-matrix(NA,nrow=length(allKEGG),ncol=nrow(H),dimnames=list(allKEGG,paste0("Sig",1:urank)))

Reacmat<-matrix(NA,nrow=length(allREAC),ncol=nrow(H),dimnames=list(allREAC,paste0("Sig",1:urank)))

GOBPmat<-matrix(NA,nrow=length(allGOBP),ncol=nrow(H),dimnames=list(allGOBP,paste0("Sig",1:urank)))
GOBPmat2<-matrix(NA,nrow=length(allGOBP2),ncol=nrow(H),dimnames=list(allGOBP2,paste0("Sig",1:urank)))

for(i in 1:nrow(H)){
  KEGGmat[KEGGenrichNMF[[i]][,"term_name"],i]<-(-log(KEGGenrichNMF[[i]][,"p_value"]))
  Reacmat[ReactomeEnrichNMF[[i]][,"term_name"],i]<-(-log(ReactomeEnrichNMF[[i]][,"p_value"]))
  
  GOBPmat[GOBPfilt[[i]][,"term_name"],i]<-(-log(GOBPfilt[[i]][,"p_value"]))
  GOBPmat2[GOBPfilt2[[i]][,"term_name"],i]<-(-log(GOBPfilt2[[i]][,"p_value"]))
}
outputdataU<-paste0(outputdata,"/NMFnqG",urank,"/")
if(!dir.exists(outputdataU)){dir.create(outputdataU)}

pdf(paste0(outputdataU,"/HeatmapGOBP_SigEnrichPval_filt15.pdf"),useDingbats = FALSE)
print(pheatmap(GOBPmat2,fontsize_row=6,cluster_cols=FALSE,cluster_rows=FALSE,treeheight_col = 0,show_colnames=F,treeheight_row = 0))
dev.off()

pdf(paste0(outputdataU,"/HeatmapKEGG_SigEnrichPval.pdf"),useDingbats = FALSE)
print(pheatmap(KEGGmat,cluster_cols=FALSE,cluster_rows=FALSE,treeheight_col = 0,show_colnames=F,treeheight_row = 0))
dev.off()

png(paste0(outputdataU,"/HeatmapKEGG_SigEnrichPval.png"),width=5,height=5,units="in",res=300,pointsize = 12)
print(pheatmap(KEGGmat,cluster_cols=FALSE,cluster_rows=FALSE,treeheight_col = 0,show_colnames=F,treeheight_row = 0))
dev.off()

pdf(paste0(outputdataU,"/HeatmapReactome_SigEnrichPval.pdf"),useDingbats = FALSE)
print(pheatmap(Reacmat,fontsize_row=6,cluster_cols=FALSE,cluster_rows=FALSE,treeheight_col = 0,show_colnames=F,treeheight_row = 0))
dev.off()

png(paste0(outputdataU,"/HeatmapReactome_SigEnrichPval.png"),width=5,height=5,units="in",res=300,pointsize = 12)
print(pheatmap(Reacmat,fontsize_row=4,cluster_cols=FALSE,cluster_rows=FALSE,treeheight_col = 0,show_colnames=F,treeheight_row = 0))
dev.off()

pdf(paste0(outputdataU,"/HeatmapGOBP_SigEnrichPval.pdf"),useDingbats = FALSE)
print(pheatmap(GOBPmat,fontsize_row=6,cluster_cols=FALSE,cluster_rows=FALSE,treeheight_col = 0,show_colnames=F,treeheight_row = 0))
dev.off()

png(paste0(outputdataU,"/HeatmapGOBP_SigEnrichPval.png"),width=5,height=5,units="in",res=300,pointsize = 12)
print(pheatmap(GOBPmat,cluster_cols=FALSE,cluster_rows=FALSE,treeheight_col = 0,show_colnames=F,treeheight_row = 0))
dev.off()




annotSigs<-t(H)
colnames(annotSigs)<-paste0("Signature",1:3)
BinaryDep<-bDepletions[setdiff(topLRTgenes,union(PanCancerCoreFitnessGenes,CommonEssentialsAUC$cfgenes)),]



annotCL<-data.frame(CT=cmp[match(colnames(bDepletions),cmp[,"model_id"]),"cancer_type"],
                    Signature1=annotSigs[,1],Signature2=annotSigs[,2],Signature3=annotSigs[,3],stringsAsFactors = FALSE)

rownames(annotCL)<-colnames(bDepletions)



pdf(paste0(outputdata,"HeatmapDependencies.pdf"),width=3,height=5)
print(pheatmap(BinaryDep,show_rownames = FALSE,show_colnames = FALSE,treeheight_row=0,treeheight_col=0,legend=FALSE,annotation_col=annotCL,annotation_legend = FALSE,color=colorRampPalette(c("white", "red","navy"))(50)))
dev.off()
png(paste0(outputdata,"BinaryHeatmapDependencies.png"),width=3,height=5,res=300,pointsize = 12,units="in")
print(pheatmap(BinaryDep,show_rownames = FALSE,show_colnames = FALSE,treeheight_row=0,treeheight_col=0,annotation_col=annotCL,legend=FALSE,annotation_legend = FALSE,color=colorRampPalette(c("white", "red","navy"))(50)))
dev.off()



