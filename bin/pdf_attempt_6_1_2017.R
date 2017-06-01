#Load packages required for the project

library(rmarkdown)
library(knitr)
library(ggplot2)
library(reshape2)
library(tibble)
library(data.table)
library(vegan)

#set to working directory to project zebra directory
setwd("C:/Users/arzan/Desktop/Bioinformatics/MCT/zebra-master/zebra-master")

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
taxa = taxa[,colSums(taxa)>=20000]
staxa = taxa                                         #staxa is utilized for relative abundances, taxa for counts
staxa = sweep(staxa,2,colSums(staxa),'/');           # Normalize
staxa = staxa[order(rowMeans(staxa),decreasing=F),]; # Sort by avg. abundance
staxa = staxa[rowMeans(staxa) >= 0.0001,];            # Drop rare taxa (abundance)
staxa = staxa[rowSums(staxa > 0) > 10,];             # Drop rare taxa (prevalence)


### Get summary infomration for the pdf doc specific to nutrition
for (id in unique(map$UserName)){
  submap <- map[map$UserName == id,]
  subtaxa <- staxa[(colnames(staxa) %in% submap[,"X.SampleID"])]
  #render(input = "lib/MCTS_pdf_mcb.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
  #output_dir = "output/5.31.test2") 
}


###Instantiate Plots###


for (id in unique(map$UserName)){
  submap <- map[map$UserName == id,]
  subtaxa <- staxa[(colnames(staxa) %in% submap[,"X.SampleID"])]
  subtaxasp <- setDT(subtaxa, keep.rownames=TRUE) 
  
  #melt subtaxasp to get our dataframe in the long format for future usage
  meltdf<- melt(subtaxasp)
  
  #merge to get access to Day var
  mergedf<- merge(x=meltdf, y=map, by.x = "variable", by.y= "X.SampleID", all.x=TRUE)
  
  #convert our dataframe species (rn) column to a character
  mergedf$rn <- as.character(mergedf$rn) 
  
  #series of gsub commands meant to neaten and clarify legend content
  mergedf$rn <- gsub(".*s__", "", mergedf$rn)
  mergedf$rn <- gsub("\\[", "",mergedf$rn)
  mergedf$rn <- gsub("\\]", "",mergedf$rn)
  mergedf$rn <- gsub("_", " ",mergedf$rn)
  
  
  #create <10% abundance category
  mergedf$rn[mergedf$value < 0.1] <- "<10% abundance"
  
 
  
 #render(input = "lib/MCTS_pdf_mcb.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         #output_dir = "output/5.31.test3") 
}


###alpha diversity###

taxa = taxa[,colSums(taxa)>=20000] #use taxa var to create OTU table (counts)

#for loop to yield subject alpha diversity to use for graphing
for (id in unique(map$UserName)){
  submapalpha <- map[map$UserName == id,] 
  subtaxaalpha <- taxa[(colnames(taxa) %in% submapalpha[,"X.SampleID"])]
  
  #transpose
   subtaxaalpha <- t(subtaxaalpha)
 
  #calculate alpha diversity and assign to variable alphad
   alphad <- diversity(subtaxaalpha, index = "shannon", MARGIN = 1, base = exp(1))
   alphad <- as.data.frame(alphad)
   alphad <- rownames_to_column(alphad, var = "X.SampleID")
   
   #merge melted df with map to get access to StudyDay variable
   merged_alpha<- merge(x=alphad, y=map, all.x=TRUE)
   
 render(input = "lib/MCTS_pdf_mcb.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         output_dir = "output/alpha_d_output") 
}



