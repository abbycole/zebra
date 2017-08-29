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
library(cowplot)


####Pre Processing####

setwd("C:/Users/arzan/Desktop/Bioinformatics/MCT/GIT/zebra")

#### LOAD PERTINENT DATA FILES ####

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

# make sure studydayno is a factor for later plotting
map$StudyDayNo <- as.factor(map$StudyDayNo)

#read in the food table and assign it to pertinent variable
#food_map <-  read.table("raw/mct.food.otu.txt", sep = "\t", header = TRUE, comment = "")
# changed to use dehydrated food weights
food_map <-  read.table("raw/mct.dhydrt.otu.txt", sep = "\t", header = TRUE, comment = "")

# read in the food taxonomy mapping file 
food_taxa <- read.table("raw/mct.taxonomy.txt", sep = "\t", header = TRUE, comment = "")



#### INSTANTIATE BETA DIVERSITY MATRIX ####


taxa <- read.delim("raw/youbet.txt", row =1)

#### Summarizing ####
split <- strsplit(rownames(taxa),";")                               # Split and rejoin on lv7
taxaStrings <- sapply(split,function(x) paste(x[1:7],collapse=";")) # level 7 is species, 8 is strain
taxa <- rowsum(taxa,taxaStrings)                                    # Collapse by taxonomy name

#### Massaging ####
taxa = taxa[,colSums(taxa)>=20000]                  #use taxa var to create OTU table (counts)
staxa = taxa                                         #staxa is utilized for relative abundances, taxa for counts
staxa = sweep(staxa,2,colSums(staxa),'/');           # Normalize
staxa = staxa[order(rowMeans(staxa),decreasing=F),]; # Sort by avg. abundance
staxa = staxa[rowMeans(staxa) >= 0.0001,];            # Drop rare taxa (abundance)
staxa = staxa[rowSums(staxa > 0) > 10,];             # Drop rare taxa (prevalence)

betataxa<- staxa[(colnames(staxa) %in% map[,"X.SampleID"])] #Create table with just taxa and subjects
betataxa <- t(betataxa) #transpose 

betad<-vegdist(betataxa, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
               na.rm = FALSE) 

#betad as matrix 
betad <- as.matrix(dist(betad))

#column name is "IDs"
IDs <-colnames(betad)




#### ALTER FOOD MAP FILE, PREPARE FOR MERGE WITH MAP ####


#set the row names of the food map file to the taxonomy strings
rownames(food_map) <- food_map$taxonomy

#remove the taxonomy string column and the food names column (redundant)
food_map$taxonomy <- NULL #removing taxa string column
food_map$X.FOODID <- NULL #removing food names column

#drop outlier nutrition samples from food_map
food_map <- food_map[,!(colnames(food_map) %in% outliers)]

split <- strsplit(rownames(food_map),";")             # Split and rejoin on lv7 to get species level
foodStrings <- sapply(split,function(x) paste(x[1:2],collapse=";"))
food_map <- rowsum(food_map,foodStrings) 

food_map <- sweep(food_map, 2, colSums(food_map), "/")

food_map <- t(food_map)

#Change the first column name of map to SampleID instead of X.SampleID
colnames(map)[1] <- "SampleID" 

#SampleIDs are also put in the row names
rownames(map) <- map$SampleID

#Use SetDT function to bring subject names in food map to col 1 for merging
food_map <- as.data.frame(food_map)
food_map <- setDT(food_map, keep.rownames = TRUE)
colnames(food_map)[1] <- "SampleID" 


#### MERGE MAP AND FOOD FILES (Food file) ####

#merge food and map
merged_food<- merge(x=food_map, y=map,by= "SampleID", all.x=TRUE)

#Subject IDs to colnames for merged_food
rownames(merged_food) <- merged_food$SampleID

#Assign common sampleIDs (sampleIDs shared between beta and metadata tables) to  IDs_keep.
#This will only keep samples shared in both tables
IDs_keep <- intersect(rownames(betad),rownames(merged_food))

#keep only samples shared in both
merged_food <- merged_food[IDs_keep,] #error here, perhaps
betad <- betad[IDs_keep,]

#We can now turn beta into a resemblence matrix  using as.dist()
beta_dist <- as.dist(betad)

#Merged_foods to data frame

merged_food <- as.data.frame(merged_food)



####Nutrition Table  ####

#Nutrition Table formatting
colnames(nutrition_table)[1] <- "SampleID" 

rownames(nutrition_table) <- nutrition_table$SampleID
#SampleID of metadata (1,2,3) are now put in rows rather than in the column


#Change all NAs in nutrition table and betad to 0
nutrition_table[is.na(nutrition_table)] <- 0
betad[is.na(betad)] <- 0


#Assign common sampleIDs (sampleIDs shared between beta and metadata tables) to  IDs_keep.
#This will only keep samples shared in both tables
IDs_keep_ntr <- intersect(rownames(betad),rownames(nutrition_table))
nutrition_table <- nutrition_table[IDs_keep_ntr,]
betad_ntr <- betad[IDs_keep_ntr,]


#We can now turn beta into a resemblence matrix  using as.dist()
beta_dist <- as.dist(betad_ntr)


#### BEGIN NUTRITION ADONIS ANALYSIS ####

adonis_PF_POULT <- adonis(betad_ntr ~ nutrition_table[,"PF_POULT"], data=nutrition_table, permutations = 999)
adonis_PF_POULT

adonis_dairy <- adonis(betad_ntr ~ nutrition_table[,"D_TOTAL"], data=nutrition_table, permutations = 999)

adonis_milk <- adonis(betad_ntr ~ nutrition_table[,"D_MILK"], data=nutrition_table, permutations = 999)

adonis_yogurt <- adonis(betad_ntr ~ nutrition_table[,"D_YOGURT"], data=nutrition_table, permutations = 999)

adonis_cheese <- adonis(betad_ntr ~ nutrition_table[,"D_CHEESE"], data=nutrition_table, permutations = 999)


#loop for diff taxa

ntaxa <- staxa
pvals_ntr <- c()

for(i in 1:nrow(ntaxa)){
  y <- as.numeric(ntaxa[i,])
  #kw_test <- kruskal.test(y ~ nutrition_table[,"D_TOTAL"])
  #kw_test <- kruskal.test(y ~ nutrition_table[,"D_MILK"])
  #kw_test <- kruskal.test(y ~ nutrition_table[,"D_YOGURT"])
  kw_test <- kruskal.test(y ~ nutrition_table[,"D_CHEESE"])
  pvals_ntr[i] <- kw_test$p.value
}

#view how many significant p values there are
sum(pvals_ntr <0.05)

#add pvalues next to taxa names

taxa_names <- rownames(staxa)
taxa_p <- as.data.frame(taxa_names)   
taxa_p$pvalues <- pvals_ntr

taxap <- data.frame(taxa_p[,-1], row.names=taxa_p[,1])
colnames(taxap)[1] <- "Pvalue"

sorted_taxap<- taxap[order(taxap$Pvalue), , drop = FALSE]
rownames(sorted_taxap[1])






#### BEGIN ADONIS FOOD ANALYSIS ####

#Food Analysis (level 2)


Cheeses <- adonis(beta_dist ~ merged_food[,"L1_Milk_and_Milk_Products;L2_Cheeses"], data=merged_food, permutations = 999)
Cheeses

Creams <- adonis(beta_dist ~ merged_food[,"L1_Milk_and_Milk_Products;L2_Creams_and_cream_substitutes"], data=merged_food, permutations = 999)
Creams

Dessert_sauce <- adonis(beta_dist ~ merged_food[,"L1_Milk_and_Milk_Products;L2_Milk_desserts_sauces_gravies"], data=merged_food, permutations = 999)
Dessert_sauce 

Milk <- adonis(beta_dist ~ merged_food[,"L1_Milk_and_Milk_Products;L2_Milks_and_milk_drinks" ], data=merged_food, permutations = 999)
Milk




#### Taxa Testing ####

#staxa <- staxa[,IDs_keep]

#ignava <- as.numeric(staxa[2,])
#hist(ignava, br=10)
#y <- as.numeric(ignava)
#kw_test <- kruskal.test(y ~ merged_food[,"L1_Milk_and_Milk_Products;L2_Milks_and_milk_drinks"])
#kw_test


#test all taxa for Milk classification using a loop

pvals <- c()

for(i in 1:nrow(staxa)){
  y <- as.numeric(staxa[i,])
  kw_test <- kruskal.test(y ~ merged_food[,"L1_Milk_and_Milk_Products;L2_Milks_and_milk_drinks"])
  pvals[i] <- kw_test$p.value
}

#view how many significant p values there are
sum(pvals <0.05)

#apply false discovery rate correction
pvals.fdr = p.adjust(pvals,"fdr")
sum(pvals.fdr < 0.05)


#plot lowest 5
#which(pvals <0.05)
plot.new()
par(mfrow=c(1,5))
lowest_five <-sort(pvals, decreasing = F)

for(i in 1:length(lowest_five)){
  z = as.numeric(staxa[i,])
  name1 = gsub('; .__',' ', rownames(staxa)[i])
  name2 = gsub('k__','', name1)
  name3 = gsub('  ', '', name2)
  name = strsplit(name3, " ", fixed = T)[[1]]
  names_tail = tail(name, n=2)
  boxplot(z ~merged_food[,"L1_Milk_and_Milk_Products;L2_Milks_and_milk_drinks"],
          main= names_tail,
          ylab = "Relative Abundance of most sig. taxa species pertaining to L2 Milk",
          col = c("blue","red","gold","green","black"),
          xaxt = "n")
}


z <- as.numeric(staxa[9,]) # can use a loop here as shown in pdf (plug in i for 9)


boxplot(z ~merged_food[,"L1_Milk_and_Milk_Products;L2_Milks_and_milk_drinks"],width = .001)







#Bowel Movement Frequency Evaluation

fecalsum <- c()
missed_stool <- c()
submap$fecal.status[is.na(submap$fecal.status)]


for (id in unique(map$UserName)){
  
  submap <- map[map$UserName == id,] #create submap variable we can use to iterate through 
  submap$fecal.status[is.na(submap$fecal.status)] <- 0 #NAs to 0 - an assumption that will be addressed later
  fecalsum[id]<- sum(submap$fecal.status) #sums amount of "ones", sum of all reported bowel movements for subejct
  missed_stool[id] <- sum(submap$fecal.status %in% 0) #sums amount of missed stools - perhaps a better metric to explore motility 

  }

hist(missed_stool) #use to qualify low vs high motility, as designated by Hypothesis 1
abline(v=median(missed_stool),col="red")


hist(fecalsum)
abline(v=median(fecalsum),col="gold")

missed_stool <- as.data.frame(missed_stool)

