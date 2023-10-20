# Score2Analysis
Scripts for producing figures included in Second-generation map of genetic vulnerabilities in cancer cells

The main script Master_Script.R should be used to run the analyses. The script assumes you run the analysis from the directory
containing this GitHub repository. The inputdata variable needs to be set to the directory where you have downloaded all
required inputs. The required inputs are available at this Figshare repository: DOI 10.6084/m9.figshare.24411922

Workflow in R:

inputdata<-"/path/to/the/folder/containing/figshare/data"
source("Master_Script.R")
