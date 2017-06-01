library(rmarkdown)
library(knitr)
library(ggplot2)
library(reshape2)
library(tibble)
library(data.table)
library(vegan)

setwd("C:/Users/arzan/Desktop/Bioinformatics/MCT/zebra-master/zebra-master") # set to zebra directory

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
taxa = taxa[order(rowMeans(taxa),decreasing=F),]; # Sort by avg. abundance
taxa = taxa[rowMeans(taxa) >= 0.0001,];            # Drop rare taxa (abundance)
taxa = taxa[rowSums(taxa > 0) > 10,];             # Drop rare taxa (prevalence)


### Get summary infomration for the pdf doc specific to nutrition
for (id in unique(map$UserName)){
  submap <- map[map$UserName == id,]
  subtaxa <- taxa[(colnames(taxa) %in% submap[,"X.SampleID"])]
  render(input = "lib/MCTS_pdf.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
         output_dir = "output") 
}

for (id in unique(map$UserName)){
  submap <- map[map$UserName == id,]
  subtaxa <- taxa[(colnames(taxa) %in% submap[,"X.SampleID"])]
  #render(input = "lib/MCTS_pdf.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
        # output_dir = "output") 
}


meltdf<- melt(subtaxa) #melt to get in long format

mergedf <- merge(x=meltdf, y=map, by.x = "variable", by.y= "X.SampleID", all.x=TRUE) #merge to get access to Day var

mergedf2$rn <- as.character(mergedf2$rn) #to char
mergedf2$rn[mergedf2$value < 0.1] <- "<10% abundance"

mergedf2$rn <- gsub(".*s__", "", mergedf2$rn)
mergedf2$rn <- gsub("\\[", "",mergedf2$rn)
mergedf2$rn <- gsub("\\]", "",mergedf2$rn)
mergedf2$rn <- gsub("_", " ",mergedf2$rn)


###Begin Plotting###

ggplot(mergedf2, aes(x = StudyDayNo, y = value, fill = rn)) +
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

#alpha diversity



###experimentation####

subtaxa2 <- setDT(subtaxa, keep.rownames=TRUE)
meltdf2<- melt(subtaxa2) #melt to get in long format

mergedf2 <- merge(x=meltdf2, y=map, by.x = "variable", by.y= "X.SampleID", all.x=TRUE) #merge to get access to Day var




abund_species <- ""
name <- ""


for (i in seq_len(ncol(subtaxa)))
{
  abund <- subtaxa[,i] #assign abundance values to variable
  abund_species[i] <- head(abund,10)
  name[i] <- rownames(subtaxa[i,])
  names(abund_species) = name
 #names(abund_species[i]) = rownames(subtaxa[i,]) #assigning names..
 
}

#for (i in seq_along(abund_species))
#{
 # names(abund_species[i]) = rownames(sub) 
#}


dfm <- melt(submap[,c('StudyDayNo','PROT','TFAT','FIBE','CAFF','FIBE','SUGR')],id.vars = 1)
ggplot(dfm,aes(x = StudyDayNo,y = value)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "stack") + scale_x_discrete(labels = abbreviate)


# make ggplot bar chart of top 10 most abundant species per day 

 colsub <-colnames(subtaxa)
 colone <- colsub[1]
 
 ggplot(subtaxa,aes(x=subtaxa[,1], y=colone))+
 geom_bar(stat="identity")
 
 geom_bar(mapping = NULL, data = NULL, stat = "count",
          position = "stack", ..., width = NULL, binwidth = NULL, na.rm = FALSE,
          show.legend = NA, inherit.aes = TRUE)


ggplot(subtaxa,geom_point(aes(x=subtaxa[,1], y=colone, color=red)))

test <- subtaxa[1,]








mcb graph a
y = relative abundance
x= days

mcb graph backsolve(
  y=alpha div
  x days
  
  pkg reshape2
)
