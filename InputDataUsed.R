

#################################################################################################
#                                                                                               #
#   Load CRISPR depletion and get coverage of models with omic data sets                        #
#                                                                                               #
#                                                                                               #
#################################################################################################

#####PAN CANCER ######

MASTER_LIST<-read.csv(paste0(inputdata,'model_list_latest.csv'),header=T,stringsAsFactors = FALSE)
M2<-MASTER_LIST
M2$model_id<-make.names(M2$BROAD_ID)
Annot<-rbind(MASTER_LIST,M2)
load(paste0(inputdata,"/07_EssMatrix_bDepletionsB2.rdata"))
bDepletions<-bDepletionsB2

#Binary features
TOTALBEM<-readRDS(file=paste0(inputdata,"ConsBEM.Rds"))
Cdrivers<-rownames(TOTALBEM)[grep("_mut",rownames(TOTALBEM))]
Cdrivers<-unlist(sapply(Cdrivers,function(x) strsplit(x,"_mut",fixed=T)[[1]][1]))
write.table(Cdrivers,file=paste0(outputdata,"/SupplementaryTable4.tsv"),quote=F,sep="\t")
temp<-CLnameMapping(bDepletions,TOTALBEM,"col",annotation=MASTER_LIST)
UseBinary<-temp$inputdata


#Gene Expression data
load(file=paste0(inputdata,"EXPpcBC.Rdata"))
EXP<-EXPpcBC
temp<-CLnameMapping(bDepletions,EXP,"col",annotation=MASTER_LIST)
UseExpr<-temp$inputdata

#Copy number data
CN<-readRDS(file=paste0(inputdata,"Combined_CNdriver_29.7.22.Rds"))
temp<-CLnameMapping(bDepletions,CN,"col",annotation=MASTER_LIST)
UseCN<-temp$inputdata

#Metabolite data
Met<-readRDS(file=paste0(inputdata,"MetaboliteData.Rds"))
temp<-CLnameMapping(bDepletions,Met,"col",annotation=MASTER_LIST)
UseMet<-temp$inputdata

#Proteomics data
Prot<-readRDS(file=paste0(inputdata,"ProteomicMarkers.Rds"))
temp<-CLnameMapping(bDepletions,Prot,"col",annotation=MASTER_LIST)
UseProt<-temp$inputdata

CRISPR<-make.names(colnames(bDepletions))
Binary<-colnames(UseBinary)
Expr<-colnames(UseExpr)
CN<-colnames(UseCN)
Met<-colnames(UseMet)
Prot<-colnames(UseProt)
plotdata<-rbind(cbind(model=CRISPR,dataset="CRISPR",cancerType=Annot[match(CRISPR,Annot$model_id),"cancer_type"]),
                cbind(model=Binary,dataset="Mutation",cancerType=Annot[match(Binary,Annot$model_id),"cancer_type"]),
                cbind(model=Expr,dataset="Expression",cancerType=Annot[match(Expr,Annot$model_id),"cancer_type"]),
                cbind(model=CN,dataset="Copy Number",cancerType=Annot[match(CN,Annot$model_id),"cancer_type"]),
                cbind(model=Met,dataset="Metabolite",cancerType=Annot[match(Met,Annot$model_id),"cancer_type"]),
                cbind(model=Prot,dataset="Proteomic",cancerType=Annot[match(Prot,Annot$model_id),"cancer_type"]))


pdf(paste0(outputdata,"/barchartDataInput.pdf"),useDingbats = FALSE)
print(ggplot(data=as.data.frame(plotdata),aes(x=cancerType,fill=dataset))+geom_bar(position="dodge",stat="count")+
        theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

png(paste0(outputdata,"/barchartDataInput.png"))
print(ggplot(data=as.data.frame(plotdata),aes(x=cancerType,fill=dataset))+geom_bar(position="dodge",stat="count")+
        theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

#generate info for supp table:
datamat<-matrix(0,nrow=length(unique(plotdata[,"model"])),ncol=6)
rownames(datamat)<-unique(plotdata[,"model"])
colnames(datamat)<-c("CRISPR","Mutation","Gene Expression","Copy Number","Metabolite","Proteomic")
datamat[CRISPR,"CRISPR"]<-1
datamat[Binary,"Mutation"]<-1
datamat[Expr,"Gene Expression"]<-1
datamat[CN,"Copy Number"]<-1
datamat[Met,"Metabolite"]<-1
datamat[Prot,"Proteomic"]<-1


write.table(datamat,file=paste0(outputdata,"/SupplementaryTable5.tsv"),quote=F,sep="\t")

pdata<-data.frame(dataset=c("CRISPR","Binary","Expression","Copy Number","Metabolite","Proteomic"),
                            counts=c(sum(plotdata[,"dataset"]=="CRISPR"),
                  sum(plotdata[,"dataset"]=="Mutation"),
                  sum(plotdata[,"dataset"]=="Expression"),
                  sum(plotdata[,"dataset"]=="Copy Number"),
                  sum(plotdata[,"dataset"]=="Metabolite"),
                  sum(plotdata[,"dataset"]=="Proteomic")),
                  stringsAsFactors = FALSE)
pdata$dataset<-factor(pdata$dataset,levels=c("Proteomic","Metabolite","Copy Number","Binary","Expression","CRISPR"))
pdf(paste0(outputdata,"/barchartDataInputPC.pdf"),useDingbats = FALSE)
print(ggplot(data=as.data.frame(pdata),aes(x=dataset,fill=dataset,y=counts))+geom_bar(position="dodge",stat="identity")+
        theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_flip())
dev.off()

png(paste0(outputdata,"/barchartDataInputPC.png"))
print(ggplot(data=as.data.frame(pdata),aes(x=dataset,fill=dataset,y=counts))+geom_bar(position="dodge",stat="identity")+
        theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_flip())
dev.off()


load(paste0(inputdata,'manifestCombined.Rdata'))


cmpcombined<-Annot
manifest$TISSUE<-make.names(cmpcombined[match(make.names(manifest$INSTITUTE_ID),make.names(cmpcombined$model_id)),"tissue"])
manifest$CANCER_TYPE<-make.names(cmpcombined[match(make.names(manifest$INSTITUTE_ID),make.names(cmpcombined$model_id)),"cancer_type"])

##### CANCER TYPE SPECIFIC COVERAGE ######
source("./SelectCLgroups.R")

pdata<-plotdata[plotdata[,1]%in%make.names(manifest$INSTITUTE_ID),]
pdata<-pdata[!pdata[,"cancerType"]=="Other Solid Cancers",]
pdata<-as.data.frame(pdata)
pdata[,"dataset"]<-factor(pdata[,"dataset"],levels=c("Proteomic","Metabolite","Copy Number","Mutation","Expression","CRISPR"))
nct<-table(pdata$cancerType)
nct<-as.vector(nct)
names(nct)<-unique(sort(pdata$cancerType))
pdata[,"cancerType"]<-factor(pdata[,"cancerType"],levels=names(nct)[order(nct,decreasing=FALSE)])
save(pdata,file=paste0(outputdata,"/pdata.Rdata"))
gout<-ggplot(data=as.data.frame(pdata),aes(x=cancerType,fill=dataset))+geom_bar(position="dodge",stat="count")+theme_bw()+coord_flip()
ggsave(filename="barchartDataInputUseCT.pdf", plot=gout, device="pdf", path=outputdata)

ggsave(filename="barchartDataInputUseCT.png", plot=gout, device="png", path=outputdata)

load(paste0(inputdata,"TissueColours.Rdata"))
save(TissueColours,file=paste0(outputdata,'TissueColours.Rdata'))
manifestColors<-manifest[,c("TISSUE","CL_GROUP")]
manifestColors<-unique(manifestColors)
names(TissueColours)<-make.names(names(TissueColours))
manifestColors$TissueColor<-TissueColours[manifestColors$TISSUE]
manifestColors<-manifestColors[manifestColors$CL_GROUP!="Other.Solid.Cancers",]
tissuenames<-unique(manifestColors$TISSUE)
tissuenames<-sort(tissuenames)
numberCT<-as.vector(table(manifestColors$TISSUE))
names(numberCT)<-tissuenames

manifestColors$CTcolor<-NA

for(i in tissuenames){
  if(numberCT[i]==1){
    manifestColors[which(manifestColors$TISSUE==i),"CTcolor"]<-manifestColors[which(manifestColors$TISSUE==i),"TissueColor"]
  }else{
    numberSubcols<-numberCT[i]
    origcol<-manifestColors[which(manifestColors$TISSUE==i),"TissueColor"][1]
    transvalues<-seq(50,100,length.out=numberSubcols)
    newCols<-sapply(transvalues,function(x) makeTransparent(origcol,x))
    names(newCols)<-manifestColors[which(manifestColors$TISSUE==i),"CL_GROUP"]
    for(j in 1:length(newCols)){
      manifestColors[which(manifestColors$CL_GROUP==names(newCols)[j]),"CTcolor"]<-newCols[j]
    }
  }
}
rownames(manifestColors)<-NULL

CancerTypeColors<-manifestColors$CTcolor
names(CancerTypeColors)<-manifestColors$CL_GROUP
CancerTypeColors<-c(CancerTypeColors,"Other.Solid.Cancers"="#9EC552")
save(CancerTypeColors,file=paste0(outputdata,"CancerTypeColors.Rdata"))

mypalette<-brewer.pal(12,"Paired")
OmicColors<-mypalette[c(1:5,7:9)]
names(OmicColors)<-c("Genomic","Copy Number","Expression","Metabolite","Proteomic","Epigenetic","Composite","Other")

save(OmicColors,file=paste0(outputdata,"OmicColors.Rdata"))
