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
taxa = sweep(taxa,2,colSums(taxa),'/');           # Normalize
taxa = taxa[order(rowMeans(taxa),decreasing=F),]; # Sort by avg. abundance
taxa = taxa[rowMeans(taxa) >= 0.0001,];            # Drop rare taxa (abundance)
taxa = taxa[rowSums(taxa > 0) > 10,];             # Drop rare taxa (prevalence)


### Get summary infomration for the pdf doc specific to nutrition
for (id in unique(map$UserName)){
  submap <- map[map$UserName == id,]
  subtaxa <- taxa[(colnames(taxa) %in% submap[,"X.SampleID"])]
  #render(input = "lib/MCTS_pdf_mcb.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         #output_dir = "output/5.31.test2") 
}

#species (sp) of subtaxa made the first column rather than being rownames
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



###Begin Plotting###

ggplot(mergedf, aes(x = StudyDayNo, y = value, fill = rn)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = FALSE) +
  theme_classic() +
  theme(strip.text.y = element_text(angle = 0, size = 8, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "grey")) +
  guides(fill = guide_legend(reverse = TRUE, 
                             keywidth = 1, 
                             keyheight = 1, 
                             ncol = 1)) +
  ylab("Relative Abundance\n") +
  ggtitle("Main species within your gut per day")

for (id in unique(map$UserName)){
  submap <- map[map$UserName == id,]
  subtaxa <- taxa[(colnames(taxa) %in% submap[,"X.SampleID"])]
  render(input = "lib/MCTS_pdf_mcb.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         output_dir = "output/5.31.test3") 
}

###alpha diversity###



