library(rmarkdown)
library(knitr)
library(ggplot2)
library(reshape2)
library(tibble)

setwd("/Users/abby/Documents/Projects/zebra") # set to zebra directory

#### Read in mapping file ####
map <- read.table("raw/new_map_with_treatment.txt", sep = "\t", header = TRUE, comment = "") #load mapping file

#### Read in taxa table #### 
taxa <- read.delim("raw/youbet.txt", row =1)

#### Summarizing ####
split <- strsplit(rownames(taxa),";")            # Split and rejoin on lv7
taxaStrings <- sapply(split,function(x) paste(x[1:7],collapse=";")) # level 7 is species, 8 is strain
taxa <- rowsum(taxa,taxaStrings)                  # Collapse by taxonomy name

# #### Manicure the samplenames, grab latest (mct study only and only when re-running things from a new taxa table)
# colnames(taxa) = gsub(".S[0-9]+.R1.001","",colnames(taxa));      # Clean old plate IDs
# taxa = taxa[,order(colnames(taxa))];              # Sort nicely by sample ID
# taxa = taxa[,-(grep("L0",colnames(taxa))-1)];     # Keep new runs only
# colnames(taxa) = gsub(".S[0-9]+.L001.R1.001","",colnames(taxa)); # Clean new plate IDs

#### Massaging ####
taxa = taxa[,colSums(taxa)>=20000]
taxa = sweep(taxa,2,colSums(taxa),'/');           # Normalize
taxa = taxa[order(rowMeans(taxa),decreasing=T),]; # Sort by avg. abundance
taxa = taxa[rowMeans(taxa) >= 0.0001,];            # Drop rare taxa (abundance)
taxa = taxa[rowSums(taxa > 0) > 10,];             # Drop rare taxa (prevalence)


### Get summary infomration for the pdf doc specific to nutrition
for (id in unique(map$UserName)){
  submap <- map[map$UserName == id,]
  subtaxa <- taxa[(colnames(taxa) %in% submap[,"X.SampleID"])]
  render(input = "lib/MCTS_pdf.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         output_dir = "output") 
}



