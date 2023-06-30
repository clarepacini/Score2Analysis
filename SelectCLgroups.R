manifest$CL_GROUP<-NA

#load the information from cell model passports for consistency and set tissue and cancer type

manifest$CANCER_TYPE<-make.names(manifest$CANCER_TYPE)
submanifest<-manifest[,c("TISSUE","INSTITUTE_ID","CANCER_TYPE","CL_GROUP")]
submanifest<-unique(submanifest)
tissueTypes<-unique(submanifest$TISSUE)
NoCancerType<-which(submanifest$CANCER_TYPE%in%c(""))
NoCancerType<-c(NoCancerType,which(is.na(submanifest$CANCER_TYPE)))
submanifest<-submanifest[!rownames(submanifest)%in%NoCancerType,]

CTNumbers<-table(submanifest$CANCER_TYPE)
CTTypes<-names(CTNumbers)[CTNumbers>9]
manifest<-manifest[manifest$CANCER_TYPE%in%CTTypes,]
manifest$CL_GROUP<-manifest$CANCER_TYPE
