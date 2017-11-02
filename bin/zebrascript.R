#!/usr/bin/env Rscript

####Initiation of Script####


suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(cowplot))
suppressPackageStartupMessages(require(installr))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(vegan))
suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(knitr))
suppressPackageStartupMessages(require(rmarkdown))
suppressPackageStartupMessages(require(tibble))


usage = '\n Welcome to the zebrascript, more coming!' 

option_list = list(
  make_option(c('-i', '--input'),
              help='REQUIRED: Input is the mapping file, a long-form tab-delimited table
              prefiltered for missing values',
              default=NA, type = 'character'),
  
  make_option(c('-t', '--taxa_table'),
              help='Provide taxa table for analysis involving bacterial taxa.',
              default=NA, type = 'character'),
  
  make_option(c('-n', '--nutr_table'),
              help='REQUIRED: Input is the nutrition table, a long-form tab-delimited table
              prefiltered for missing values',
              default=NA, type = 'character'),
  
  make_option(c('-u', '--user_map'),
              help='REQUIRED: Input is the username mapping file, a long-form tab-delimited table
              prefiltered for missing values',
              default=NA, type = 'character'),
  
  make_option(c('-f', '--food_map'),
              help='REQUIRED: Input is the food map, a long-form tab-delimited table
              prefiltered for missing values',
              default=NA, type = 'character'))   
  
 # add_option(c('-am', '--alpha_metric')),
             #type = 'choice',
             #options = c('shannon','chao1'),
             #default = 'shannon',
             #help = 'Choose alpha diversity metric - shannon, chao1'))

opt <- parse_args(OptionParser(usage=usage, option_list=option_list))

if (is.na(opt$input)) {
  stop('Missing required parameters. See usage and options (--help)')
}

####parse command-line information to variables####
infile <- opt$input 
taxa_table <- opt$taxa_table
nutrition_table <- opt$nutr_table
user_map <- opt$user_map
food_map <- opt$food_map
#alpa_metric <- opt$alpha_metric

##### Read data provided in to pertinent tables ####
df<- read.table(file= infile, sep = "\t", header = TRUE, comment = "")
UserName_map <- read.table("../raw/UserName_map.txt", sep = "\t", header = TRUE, comment = "")
####Mapping File set up####

# designate as the SampleID map
SampleID_map <- df 

# drop the blank samples
SampleID_map <- SampleID_map[grep("Blank", SampleID_map$X.SampleID, invert = TRUE),]

# create REE and TEE variables
map <- SampleID_map %>% mutate(REE = ifelse(Gender=="M", (10*Weight+6.25*Height-5*Age+5), 
                                            ifelse(Gender=="F", (10*Weight+6.25*Height-5*Age-161),
                                                   NA)), 
                               TEE = REE*Activity.Factor) # calculate TEE with activity factor

# make sure studydayno is a factor for later plotting
map$StudyDayNo <- as.factor(map$StudyDayNo)



####Nutrition Table Set up####
if(is.na(nutrition_table)==FALSE){  #nutrition_table 
  
nutrition_table <- read.table(file=nutrition_table, sep = "\t", header = TRUE, comment = "")
#Remove outlier nutr variables
outliers = c("MCT.f.0021", "MCT.f.0044", "MCT.f.0050", "MCT.f.0056", "MCT.f.0058", "MCT.f.0060", "MCT.f.0076", "MCT.f.0116", "MCT.f.0465", "MCT.f.0486", "MCT.f.0601")
nutrition_table <- nutrition_table[!(nutrition_table$X.SampleID %in% outliers),]
} else {
  print("Nutrition Table not provided. Pertinent nutrition analyses will not be included in the report.")
}


#### Food Map Set Up ####

if(is.na(food_map)==FALSE){ #food table evaluator
  
food_map <-  read.table(file=food_map, sep = "\t", header = TRUE, comment = "")
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
} else {
  print("Food mapping file not provided - report will exclude pertinent food data")
}


#### Taxa Table Set up ####

if(is.na(taxa_table)==FALSE){ #evaluator
  
taxa <- read.delim(file=taxa_table, row =1)
  
    ### Summarizing ###
split <- strsplit(rownames(taxa),";")                               # Split and rejoin on lv7
taxaStrings <- sapply(split,function(x) paste(x[1:7],collapse=";")) # level 7 is species, 8 is strain
taxa <- rowsum(taxa,taxaStrings)                                    # Collapse by taxonomy name

### Massaging ###
taxa = taxa[,colSums(taxa)>=20000]                  #use taxa var to create OTU table (counts)

#species level taxa#

staxa = taxa                                         #staxa is utilized for relative abundances, taxa for counts
staxa = sweep(staxa,2,colSums(staxa),'/');           # Normalize
staxa = staxa[order(rowMeans(staxa),decreasing=F),]; # Sort by avg. abundance
staxa = staxa[rowMeans(staxa) >= 0.0001,];            # Drop rare taxa (abundance)
staxa = staxa[rowSums(staxa > 0) > 10,];             # Drop rare taxa (prevalence)

###genus level taxa ###

gtaxa = taxa
splitg <- strsplit(rownames(gtaxa),";")                               # Split and rejoin on lv7
gtaxaStrings <- sapply(splitg,function(x) paste(x[1:6],collapse=";")) # level 7 is species, 8 is strain
gtaxa <- rowsum(gtaxa,gtaxaStrings)                                    # Collapse by taxonomy name
gtaxa = sweep(gtaxa,2,colSums(gtaxa),'/')
gtaxa = gtaxa[order(rowMeans(gtaxa),decreasing=T),]

#selet top four
gtaxa <- gtaxa[1:4,]
#traspose to add to map for later use
gtaxa <- t(gtaxa)
#colnames(sptaxa) <- gsub("^.*\\.","", colnames(sptaxa) )
colnames(gtaxa) <- gsub(".*;g__?", "", colnames(gtaxa))
colnames(gtaxa) <- gsub("_", "", colnames(gtaxa))
colnames(gtaxa) <- gsub(";.*","",colnames(gtaxa))

#gtaxa and map append#

myrownames <- rownames(gtaxa)
gtaxa <- apply(gtaxa, 2, function(x) cut(x, breaks = seq(0, max(x), length.out = 15)))
rownames(gtaxa) <- myrownames
gtaxa <- as.data.frame(gtaxa)
gtaxa <- rownames_to_column(gtaxa, "X.SampleID")

###Instantiate Plots (with gtaxa)###

colnames(map)[2] <- "UserName" #Change col name from UserName.x to UserName - compatibility purposes
colnames(map)[3] <- "StudyDayNo" #Change col name from StudyDayNo.x to StudyDayNo - compatibility purposes


# add relative abundace of species taxa to the map for plotting
map <- left_join(map,gtaxa, by = "X.SampleID")
} else {
  print("Taxa table not provided. Pertinent taxa analyses will not be included in the report.")
}

for (id in unique(map$UserName)){
  
  #mapping
  submap <- map[map$UserName == id,] #create submap variable we can use to iterate through 
  
  #nutrition iteration
  if(is.na(nutrition_table)==FALSE){ #evaluator
  subnutr <- nutrition_table[nutrition_table$UserName == id,]
  macro_totals <-  round(mean(subnutr$PROT),2) +  round(mean(subnutr$TFAT),2) + round(mean(subnutr$CARB),2)
  } else (
    return(NULL)
  )
  
  #food iteration#
  if(is.na(food_map)==FALSE){ #evaluator
  sub_foodmap <- merged_food[merged_food$UserName == id,] #sub food variable allows us to access one subject at a time (like submap)
  sub_food <- food_map[(colnames(food_map) %in% sub_foodmap[,"variable"])]#sub food variable looks at each subject individually
  sub_foodsp <- setDT(sub_food, keep.rownames = TRUE) #make rownames first column of frame
  } else {
    print("Mission Failed")
  }
  
  #taxa iteration#
  if(is.na(taxa_table)==FALSE){ #evaluator
  subtaxa <- staxa[(colnames(staxa) %in% submap[,"X.SampleID"])] #create subtaxa variable - looks at each subject individually
  subtaxa <- subtaxa[order(rowMeans(subtaxa), decreasing = TRUE),]# sort the subtaxa df by abundance
  subtaxasp <- setDT(subtaxa, keep.rownames=TRUE) #make rownames first column of frame
  #alpha diversity iteration
  subtaxaalpha <- taxa[(colnames(taxa) %in% submap[,"X.SampleID"])] #subtaxa variable used for alpha diversity
  taxalpha <- taxa[(colnames(taxa) %in% map[,"X.SampleID"])] #matrix utilized to calculate study average alpha div.
  #beta diversity iteration
  betataxa<- staxa[(colnames(staxa) %in% map[,"X.SampleID"])] #Create table with just taxa and subjects
  betataxa <- t(betataxa) #transpose 
  } else {
    print("Mission Failed")
  }
 
  #blood iteration#
  subblood <- UserName_map[UserName_map$UserName == id,] #blood draw values
  
  #rendering
  render(input = "../lib/scripthtml.Rmd",output_file = paste0('report.', id, '.html'),"html_document",
         output_dir = "Script_output") 
}








#Write information to text files

#sink()

#write.table(map, file="map.txt", sep = "\t",row.names=FALSE)

#sink()


