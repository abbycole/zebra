#!/usr/bin/env Rscript


require(optparse)
require(dplyr)

if (!require("optparse")) {
  install.packages("optparse", repos="http://cran.rstudio.com/") 
  library("optparse")
}

if (!require("dplyr")) {
  install.packages("dplyr", repos="http://cran.rstudio.com/") 
  library("dplyr")
}


usage = '\n Welcome to the zebrascript, more coming!' 

option_list = list(
  make_option(c('-i', '--input'),
              help='REQUIRED: Input is the mapping file, a long-form tab-delimited table
              prefiltered for missing values',
              default=NA, type = 'character'))

opt <- parse_args(OptionParser(usage=usage, option_list=option_list))

if (is.na(opt$input)) {
  stop('Missing required parameters. See usage and options (--help)')
}

infile <- opt$input #parse command-line entry to infile var


# Read the data file
df<- read.table(file= infile, sep = "\t", header = TRUE, comment = "")

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


sink()

write.table(map, file="map.txt", sep = "\t",row.names=FALSE)

sink()


