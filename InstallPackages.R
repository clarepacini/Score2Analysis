options(repos = c("https://www.stats.bris.ac.uk/R/"))
options(repos = c("https://cran.rstudio.com/"))
cran_packages <- c('ggplot2','ggrepel','tidyverse','RColorBrewer','NMF','gprofiler2','pheatmap','beeswarm',
                   'reshape2', 'dplyr', 'data.table', 'Seurat','devtools','survival','GGally','stringr','heatmap3',
                   'igraph','magrittr','Matrix','corrplot','mclust','diptest','VennDiagram','circlize','remotes','pwrss')
new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

bioconductor_packages <- c('qusage','TCGAutils','rtracklayer','ComplexHeatmap','progeny','maftools')
new_bioconductor_packages <- bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
if(length(new_bioconductor_packages)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(new_bioconductor_packages)
}
library(remotes)
if(!"CRISPRcleanR"%in%installed.packages()[,"Package"]){
  remotes::install_github("francescojm/CRISPRcleanR")}

if(!"ScorePackage"%in%installed.packages()[,"Package"]){
  remotes::install_github("DepMap-Analytics/ScorePackage")}

library(ScorePackage)
library(CRISPRcleanR)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(RColorBrewer)
library(NMF)
library(gprofiler2)
library(pheatmap)
library(beeswarm)
library(reshape2)
library(dplyr)
library(data.table)
library(Seurat)
library(survival)
library(GGally)
library(stringr)
library(igraph)
library(magrittr)
library(Matrix)
library(corrplot)
library(mclust)
library(diptest)
library(VennDiagram)
library(circlize)
library(qusage)
library(TCGAutils)
library(rtracklayer)
library(ComplexHeatmap)
library(progeny)
library(maftools)
library(heatmap3)
library(pwrss)



