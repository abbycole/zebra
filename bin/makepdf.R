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
#library(extrafont)

####Pre Processing####

#set to working directory to project zebra directory
setwd("/Users/abby/Documents/Projects/zebra/")
#setwd("C:/Users/arzan/Desktop/Bioinformatics/MCT/GIT/zebra")

#Read in sampleID and nutrition files and assign them to pertinent variables
nutrition_table <- read.table("raw/Totals_to_use.txt", sep = "\t", header = TRUE, comment = "")

#Remove outlier nutr variables
outliers = c("MCT.f.0021", "MCT.f.0044", "MCT.f.0050", "MCT.f.0056", "MCT.f.0058", "MCT.f.0060", "MCT.f.0076", "MCT.f.0116", "MCT.f.0465", "MCT.f.0486", "MCT.f.0601")
nutrition_table <- nutrition_table[!(nutrition_table$X.SampleID %in% outliers),]

#Load the maps
SampleID_map <- read.table("raw/SampleID_map.txt", sep = "\t", header = TRUE, comment = "")

#Drop blanks from map
# drop the blank samples
SampleID_map <- SampleID_map[grep("Blank", SampleID_map$X.SampleID, invert = TRUE),]

UserName_map <- read.table("raw/UserName_map.txt", sep = "\t", header = TRUE, comment = "")

#merge the two to create the mapping file
#map <-  merge(x=nutrition_table, y=SampleID_map, by= "X.SampleID")

# create REE and TEE variables
map <- SampleID_map %>% mutate(REE = ifelse(Gender=="M", (10*Weight+6.25*Height-5*Age+5), 
                                   ifelse(Gender=="F", (10*Weight+6.25*Height-5*Age-161),
                                          NA)), 
                      TEE = REE*Activity.Factor) # calculate TEE with activity factor

#read in the food table and assign it to pertinent variable
food_map <-  read.table("raw/mct.food.otu.txt", sep = "\t", header = TRUE, comment = "")

# read in the food taxonomy mapping file 
food_taxa <- read.table("raw/mct.taxonomy.txt", sep = "\t", header = TRUE, comment = "")

# remove water from the food_map df
# find descriptions of water
water <- food_taxa %>% filter(FoodID != "92410110" &
                                FoodID != "92410210"& 
                                FoodID != "92410250"& 
                                FoodID != "94000100"&
                                FoodID != "94100100"&
                                FoodID != "94100200" &
                                FoodID != "94100300" &
                                FoodID != "94210100" &
                                FoodID != "94210200" &
                                FoodID != "94210300" &
                                FoodID != "94220200" &
                                FoodID != "94300100")

#remove water
food_taxa <- food_taxa %>% filter(taxonomy %in% water$taxonomy)

# subset to remove water from the food table
food_map <- food_map[food_map$X.FOODID %in% food_taxa$Main.food.description,]

# clean up environment
rm(water)

#set the row names of the food map file to the taxonomy strings
rownames(food_map) <- food_map$taxonomy

#remove the taxonomy string column and the food names column (redundant)
food_map$taxonomy <- NULL #removing taxa string column
food_map$X.FOODID <- NULL #removing food names column

#drop outlier nutrition samples from food_map
food_map <- food_map[,!(colnames(food_map) %in% outliers)]

# collapse by food level
# Summarizing at different levels - makes changes to everything downstream
split <- strsplit(rownames(food_map),";")             # Split and rejoin on lv7 to get species level
foodStrings <- sapply(split,function(x) paste(x[1:1],collapse=";"))
food_map <- rowsum(food_map,foodStrings) 

# make food relative abundance
food_map <- sweep(food_map, 2, colSums(food_map), "/")

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
  subnutr <- nutrition_table[nutrition_table$UserName == id,]
  subtaxa <- staxa[(colnames(staxa) %in% submap[,"X.SampleID"])] #create subtaxa variable - looks at each subject individually
  # sort the subtaxa df by abundance
  subtaxa <- subtaxa[order(rowMeans(subtaxa), decreasing = TRUE),]
  subtaxasp <- setDT(subtaxa, keep.rownames=TRUE) #make rownames first column of frame
  
  
  #food iteration
  
  macro_totals <-  round(mean(subnutr$PROT),2) +  round(mean(subnutr$TFAT),2) + round(mean(subnutr$CARB),2)
  sub_foodmap <- merged_food[merged_food$UserName == id,] #sub food variable allows us to access one subject at a time (like submap)
  sub_food <- food_map[(colnames(food_map) %in% sub_foodmap[,"variable"])]#sub food variable looks at each subject individually
  sub_foodsp <- setDT(sub_food, keep.rownames = TRUE) #make rownames first column of frame
  
  #alpha diversity iteration
  subtaxaalpha <- taxa[(colnames(taxa) %in% submap[,"X.SampleID"])] #subtaxa variable used for alpha diversity
 
  #beta diversity iteration
  betataxa<- staxa[(colnames(staxa) %in% map[,"X.SampleID"])] #Create table with just taxa and subjects
  betataxa <- t(betataxa) #transpose 
  
  #blood draw values
  subblood <- UserName_map[UserName_map$UserName == id,]
  
  #rendering
  render(input = "lib/mypdf.Rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         output_dir = "output/output") 
}











