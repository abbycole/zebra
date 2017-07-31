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
library(extrafont)

####Pre Processing####

#set to working directory to project zebra directory
setwd("C:/Users/arzan/Desktop/Bioinformatics/MCT/GIT/zebra")


#Read in sampleID and nutrition files and assign them to pertinent variables
nutrition_table <- read.table("raw/Totals_to_use.txt", sep = "\t", header = TRUE, comment = "")

SampleID_map <- read.table("raw/SampleID_map.txt", sep = "\t", header = TRUE, comment = "")

#merge the two to create the mapping file
map <-  merge(x=nutrition_table, y=SampleID_map, by= "X.SampleID")

#read in the food table and assign it to pertinent variable
food_map <-  read.table("raw/mct.food.otu.txt", sep = "\t", header = TRUE, comment = "")

#set the row names of the food map file to the taxonomy strings
rownames(food_map) <- food_map$taxonomy

#remove the taxonomy string column and the food names column (redundant)
food_map$taxonomy <- NULL #removing taxa string column
food_map$X.FOODID <- NULL #removing food names column

#Bring food taxa strings from row names to first column using setDT. Food_map is now analagous to subtaxasp
#food_map <- setDT(food_map, keep.rownames = TRUE)

#melt food_map to create a data drame called melted_food to merge by SampleIDs
melted_food<- melt(food_map)

#merge melted_food and map by variable (SampleID -- called variable in melted food file - so we use by.x)
merged_food<- merge(x=melted_food, y=map, by.x = "variable", by.y= "X.SampleID", all.x=TRUE)

colnames(merged_food)[3] <- "UserName" #Change col name from UserName.x to UserName - compatibility purposes
colnames(merged_food)[4] <- "StudyDayNo" #Change col name from StudyDayNo.x to StudyDayNo - compatibility purposes

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

colnames(map)[2] <- "UserName" #Change col name from UserName.x to UserName - compatibility purposes
colnames(map)[3] <- "StudyDayNo" #Change col name from StudyDayNo.x to StudyDayNo - compatibility purposes

for (id in unique(map$UserName)){
  
  #microbiome iteration
  submap <- map[map$UserName == id,] #create submap variable we can use to iterate through 
  subtaxa <- staxa[(colnames(staxa) %in% submap[,"X.SampleID"])] #create subtaxa variable - looks at each subject individually
  subtaxasp <- setDT(subtaxa, keep.rownames=TRUE) #make rownames first column of frame
  
  #food iteration
  sub_foodmap <- merged_food[merged_food$UserName == id,] #sub food variable allows us to access one subject at a time (like submap)
  sub_food <- food_map[(colnames(food_map) %in% sub_foodmap[,"variable"])]#sub food variable looks at each subject individually
  sub_foodsp <- setDT(sub_food, keep.rownames = TRUE) #make rownames first column of frame
  
  #alpha diversity iteration
  subtaxaalpha <- taxa[(colnames(taxa) %in% submap[,"X.SampleID"])] #subtaxa variable used for alpha diversity
 
  #beta diversity iteration
  betataxa<- staxa[(colnames(staxa) %in% map[,"X.SampleID"])] #Create table with just taxa and subjects
  betataxa <- t(betataxa) #transpose 
  
  #rendering
  render(input = "lib/MCTS_pdf_mcb_7_3.Rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         output_dir = "output/7.20.output") 
}











