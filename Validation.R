
clines<-c("SW626","HT29","SW837","MDST8","HCT116","KM12","RKO","LS180")
ids<-c("SIDM01168","SIDM00136","SIDM00833","SIDM00527","SIDM00783","SIDM00150","SIDM01090","SIDM00680")
ids2<-c("ACH.001399","ACH.000552","ACH.000421","ACH.000935","ACH.000971","ACH.000969","ACH.000943","ACH.000957")

valData<-read.csv(paste0(inputdata,"ColoRes.csv"),stringsAsFactors = FALSE,header=T)



####plot for WRN/BRCA2:
idmatch2<-cbind(clines,ids2)
DoutT<-valData
DoutT$cellline<-idmatch2[match(DoutT$Annotation,idmatch2[,1]),2]
DoutT$MSI<-0
DoutT$KRASmut<-0
DoutT[DoutT$cellline%in%c("ACH.000971","ACH.000969","ACH.000943","ACH.000957"),"MSI"]<-1
DoutT[DoutT$cellline%in%c("ACH.001399","ACH.000421","ACH.000971","ACH.000957"),"KRASmut"]<-1
DoutT$biomarker<-"0"
DoutT[DoutT$X=="KRAS"&DoutT$KRASmut==1,"biomarker"]<-"1"
DoutT[DoutT$X%in%c("BRCA2","WRN")&DoutT$MSI==1,"biomarker"]<-"2"
gsel<-c("WRN","BRCA2")
Dplot<-DoutT[DoutT$X%in%gsel,]
Dplot$X.4<-as.numeric(Dplot$X.4)

pdf(paste0(outputdata,"Individ_MSITargets2.pdf"),useDingbats = FALSE,width=6,height=3)
bp<-ggplot(data=Dplot,aes(x=biomarker,y=X.4,group=biomarker))+geom_boxplot()+geom_point()+
  facet_grid(.~X)+theme_bw()+
  xlab("MSI status")+ylab("Validation Score")+theme(legend.position = "none")+coord_flip()
print(bp)
dev.off()

