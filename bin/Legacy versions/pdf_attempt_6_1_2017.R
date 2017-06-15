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


###Instantiate Plots###


for (id in unique(map$UserName)){
  submap <- map[map$UserName == id,]
  subtaxa <- staxa[(colnames(staxa) %in% submap[,"X.SampleID"])]
  subtaxasp <- setDT(subtaxa, keep.rownames=TRUE) 
#   
#   #melt subtaxasp to get our dataframe in the long format for future usage
#   meltdf<- melt(subtaxasp)
#   
#   #merge to get access to Day var
#   mergedf<- merge(x=meltdf, y=map, by.x = "variable", by.y= "X.SampleID", all.x=TRUE)
#   
#   #convert our dataframe species (rn) column to a character
#   mergedf$rn <- as.character(mergedf$rn) 
#   
#   #series of gsub commands meant to neaten and clarify legend content
#   mergedf$rn <- gsub(".*s__", "", mergedf$rn)
#   mergedf$rn <- gsub("\\[", "",mergedf$rn)
#   mergedf$rn <- gsub("\\]", "",mergedf$rn)
#   mergedf$rn <- gsub("_", " ",mergedf$rn)
#   
#   
#   #create <10% abundance category
#   mergedf$rn[mergedf$value < 0.1] <- "<10% abundance"
#   
 
  
 render(input = "lib/MCTS_pdf_mcb.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         output_dir = "output/A") 
}


#### Alpha Diversity ####

taxa = taxa[,colSums(taxa)>=20000] #use taxa var to create OTU table (counts)

#for loop to yield subject alpha diversity to use for graphing
for (id in unique(map$UserName)){
  submapalpha <- map[map$UserName == id,] 
  subtaxaalpha <- taxa[(colnames(taxa) %in% submapalpha[,"X.SampleID"])]
  
#   #transpose
#    subtaxaalpha <- t(subtaxaalpha)
#  
#   #calculate alpha diversity and assign to variable alphad
#    alphad <- diversity(subtaxaalpha, index = "shannon", MARGIN = 1, base = exp(1))
#    alphad <- as.data.frame(alphad)
#    alphad <- rownames_to_column(alphad, var = "X.SampleID")
#    
#    #merge melted df with map to get access to StudyDay variable
#    merged_alpha<- merge(x=alphad, y=map, all.x=TRUE)
   
 render(input = "lib/MCTS_pdf_mcb.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         output_dir = "output/B") 
}


#### Beta Diversity ####

taxa = taxa[,colSums(taxa)>=20000] #use taxa var to create OTU table (counts)
taxa = sweep(taxa,2,colSums(taxa),'/');   
taxa = taxa[order(rowMeans(taxa),decreasing=F),]; # Sort by avg. abundance
taxa = taxa[rowMeans(taxa) >= 0.0001,];            # Drop rare taxa (abundance)
taxa = taxa[rowSums(taxa > 0) > 10,];             # Drop rare taxa (prevalence)

betataxa<- taxa[(colnames(taxa) %in% map[,"X.SampleID"])] #Create table with just taxa and subjects

betataxa <- t(betataxa) #transpose the table to get to OTU format

#calculate beta diversity and assign to variable betad - method is "w" for now
#betad <-betadiver(betataxa, method="w", order = FALSE, help = FALSE)

betad<-vegdist(betataxa, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
        na.rm = FALSE) 

#betad as matrix 
betad <- as.matrix(dist(betad))

#column name is "IDs"
IDs <-colnames(betad)

# Run the pcoa() function on the beta diversity table, and store the vectors generated as a dataframe 
PCOA <- data.frame(pcoa(betad)$vectors)

#create vector with placeholders

new_names <- rep("",ncol(PCOA))

#create for loop which gives column names (currently named "axis") the name of PC1, PC2, PC3.. so forth
#Each column for the PCOA matrix is a new PC (number).

for(i in 1:ncol(PCOA)){
  new_names[i] <- paste("PC",i,sep="")  
}

#code above uses a for loop - from the range of 1 to the number of columns in the PCOA
#matrix, the for loop accesses the new name variable and names it PC + i (which is 1,2,3  depending on column
#as indicated by for loop - for loop increases in increments of 1 (i goes up by 1) and goes thru columns)

#new_names (PC1, PC2, etc..) replace the current colnames
names(PCOA) <- new_names

#creates column that is sampleID
PCOA$X.SampleID <- IDs
PCOA$UserName <- IDs

#merge map and PCOA by SampleID - select both lines and run together - won't work otherwise
PCOA <- merge(PCOA, map, by = "X.SampleID")
PCOA <- merge(map,PCOA, by = "UserName")

#### Plot PCOA ####


#instantiate shape and color palettes 

shape_pal <- c(1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8)
shape_pal2 <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
shape_pal3 <- c(15,19,17,18,15,19,17,18,15,19,17,18,15,19,17,18,15,19,17,18,15,19,17,18,15,19,17,18,15,19,17,18,15,19,17,18,15,19,17,18)
colpal<- colorRampPalette(brewer.pal(9,"Set3"))(37)


#use ggplot to plot
ggplot(PCOA) +
  geom_point(aes(x=PC1,y=PC2, col=UserName, shape=UserName, size=4),alpha=0.2) + #gives axes and item of focus
  labs(title="PCOA plot") +
  scale_shape_manual(values=shape_pal3) +
  scale_color_manual(values=colpal) +
  theme_classic()+ 
  guides(fill = guide_legend(override.aes= list(alpha = 1/10)), colour = guide_legend(override.aes = list(size=5))) 
  #+ guides(fill = guide_legend(override.aes= list(alpha = 0.4)))

render(input = "lib/MCTS_pdf_mcb.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
       output_dir = "output/C") 





