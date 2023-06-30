source("./Functions.R")
source("./InstallPackages.R")
source("./CellignerFunctions.R")


outputdata<-"../outputdata/"
if(!dir.exists(outputdata)){dir.create(outputdata)}
inputdata<-"../InputData/"

#Overview of dependency data:
source("./DependencyResults.R")

#Dependency signature analysis:
source("./DependencySignatures.R")

#Input data:
source("./InputDataUsed.R")

#for the Celligner Pan can clusters &figures:
source("./GExpClusters.R")
source("./Celligner_Checks.R")

source("./ClinicalGeneExpressionMarkers.R")

#DFA landscape analysis
source("./DMA_Analysis.R")


#RWR analysis
source("./RWR_landscape.R")

#Priority Target Results
source("./CancerType_PriorityResults.R")

source("./PANCAN_PriorityResults.R")

#processing validation data:
source("./Validation.R")


#Paralog analysis pancancer for non-reciprocal buffering:
source("./Paralog_Analysis.R")

#Assigning patient prevalence to priority results:
source('./Target_Prevalence.R')


