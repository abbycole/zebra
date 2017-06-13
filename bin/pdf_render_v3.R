#Load packages required for the project

library(rmarkdown)
library(knitr)
library(ggplot2)
library(reshape2)
library(tibble)
library(data.table)
library(vegan)
library(ape)
library(RColorBrewer)
library(dplyr)

####Pre Processing####

#set to working directory to project zebra directory
setwd("C:/Users/arzan/Desktop/Bioinformatics/MCT/GIT/zebra")
#### Read in mapping file ####
map <- read.table("raw/new_map_with_treatment.txt", sep = "\t", header = TRUE, comment = "")

#### Read in taxa table #### 
taxa <- read.delim("raw/youbet.txt", row =1)

#### Summarizing ####
split <- strsplit(rownames(taxa),";")                               # Split and rejoin on lv7
taxaStrings <- sapply(split,function(x) paste(x[1:7],collapse=";")) # level 7 is species, 8 is strain
taxa <- rowsum(taxa,taxaStrings)                                    # Collapse by taxonomy name




 #### Manicure the samplenames, grab latest (mct study only and only when re-running things from a new taxa table)
   #colnames(taxa) = gsub(".S[0-9]+.R1.001","",colnames(taxa));      # Clean old plate IDs
   #taxa = taxa[,order(colnames(taxa))];              # Sort nicely by sample ID
   #taxa = taxa[,-(grep("L0",colnames(taxa))-1)];     # Keep new runs only
   #colnames(taxa) = gsub(".S[0-9]+.L001.R1.001","",colnames(taxa)); # Clean new plate IDs



#### Massaging ####
taxa = taxa[,colSums(taxa)>=20000]                  #use taxa var to create OTU table (counts)
staxa = taxa                                         #staxa is utilized for relative abundances, taxa for counts
staxa = sweep(staxa,2,colSums(staxa),'/');           # Normalize
staxa = staxa[order(rowMeans(staxa),decreasing=F),]; # Sort by avg. abundance
staxa = staxa[rowMeans(staxa) >= 0.0001,];            # Drop rare taxa (abundance)
staxa = staxa[rowSums(staxa > 0) > 10,];             # Drop rare taxa (prevalence)


###Instantiate Plots###

for (id in unique(map$UserName)){
  
  submap <- map[map$UserName == id,] #create submap variable we can use to iterate through 
  subtaxa <- staxa[(colnames(staxa) %in% submap[,"X.SampleID"])] #create subtaxa variable - looks at each subject individually
  subtaxasp <- setDT(subtaxa, keep.rownames=TRUE) 
  subtaxaalpha <- taxa[(colnames(taxa) %in% submap[,"X.SampleID"])] #subtaxa variable used for alpha diversity
  betataxa<- staxa[(colnames(staxa) %in% map[,"X.SampleID"])] #Create table with just taxa and subjects
  betataxa <- t(betataxa) #transpose 
  render(input = "lib/MCTS_pdf_mcb_v2.Rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         output_dir = "output/B") 

}










