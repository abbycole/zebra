
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
library(cowplot)

###Instantiation###
#set to working directory to project zebra directory
#setwd("/Users/abby/Documents/Projects/zebra/")
setwd("C:/Users/arzan/Desktop/Bioinformatics/MCT/GIT/zebra")


shinyApp(
  ui = fluidPage(
    sliderInput("slider", "Slider", 1, 100, 50),
    downloadButton("report", "Generate report"),
    
    sidebarLayout(
      
      sidebarPanel(
       
        #Import mapping file
        fileInput("file_map", "Upload mapping file",
                  buttonLabel = "Browse",
                  multiple = FALSE,
                  accept = ".txt",
                  placeholder = "Tab-delimited!"),
                  
        
        #Import nutrition table
        fileInput("file_nutr", "Upload nutrition table",
                  buttonLabel = "Browse",
                  multiple = FALSE,
                  accept = ".txt",
                  placeholder = "Tab-delimited!"),
       
        
        #Import user mapping file
        fileInput("file_user", "Upload username mapping file",
                  buttonLabel = "Browse",
                  multiple = FALSE,
                  accept = ".txt",
                  placeholder = "Tab-delimited!"),
        
        #Import taxa file
        fileInput("file_taxa", "Upload taxa table",
                  buttonLabel = "Browse",
                  multiple = FALSE,
                  accept = ".txt",
                  placeholder = "Tab-delimited!"),
      
        
        #Import food-mapping file
        fileInput("file_fmap", "Upload food map",
                  buttonLabel = "Browse",
                  multiple = FALSE,
                  accept = ".txt",
                  placeholder = "Tab-delimited!"),
       
        
        #Import food taxa file
        fileInput("file_ftaxa", "Upload food taxa",
                  buttonLabel = "Browse",
                  multiple = FALSE,
                  accept = ".txt",
                  placeholder = "Tab-delimited!")
        ),
    
      mainPanel(
        
        selectInput("alpha_metric","Choose metric for alpha diversity calculations:",
                    c("shannon", "simpson","chao1" ))
        
      )
     
    )
    ),
  
  #End UI initialization, start server initialization 
  server = function(input, output) {
  
    output$report <- downloadHandler(
      
      # For PDF output, change this to "test.pdf"
      filename = "test.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        # tempReport <- file.path(tempdir(), "test.Rmd")
        # file.copy("test.Rmd", tempReport, overwrite = TRUE)
        
        
        #Require inputs - switch on/off depending on what user can provide
        req(input$file_map) 
        req(input$file_user) 
        req(input$file_nutr) 
        req(input$file_taxa) 
        req(input$file_fmap) 
        req(input$file_ftaxa) 
        
        #Metric/User choice set up:
        
        #User input for alpha metric becomes value for alpha metric
        alpha_metric <- input$alpha_metric 
       
        
        # Set up parameters to pass to Rmd document
        params <- list(n = input$slider)
        
        #read in mapping file
        SampleID_map <- read.table(input$file_map$datapath, sep = "\t", header = TRUE, comment = "")
        
        #Read in sampleID and nutrition files and assign them to pertinent variables
        nutrition_table <- read.table(input$file_nutr$datapath, sep = "\t", header = TRUE, comment = "")
        
        #Remove outlier nutr variables
        outliers = c("MCT.f.0021", "MCT.f.0044", "MCT.f.0050", "MCT.f.0056", "MCT.f.0058", "MCT.f.0060", "MCT.f.0076", "MCT.f.0116", "MCT.f.0465", "MCT.f.0486", "MCT.f.0601")
        nutrition_table <- nutrition_table[!(nutrition_table$X.SampleID %in% outliers),]
        
        #Load the maps
        #SampleID_map <- read.table("raw/SampleID_map.txt", sep = "\t", header = TRUE, comment = "")
        
        #Drop blanks from map
        # drop the blank samples
        SampleID_map <- SampleID_map[grep("Blank", SampleID_map$X.SampleID, invert = TRUE),]
        
        UserName_map <- read.table(input$file_user$datapath, sep = "\t", header = TRUE, comment = "")
        
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
        food_map <-  read.table(input$file_fmap$datapath, sep = "\t", header = TRUE, comment = "")
        
        # read in the food taxonomy mapping file 
        food_taxa <- read.table(input$file_ftaxa$datapath, sep = "\t", header = TRUE, comment = "")
        
        # # remove water from the food_map df
        # # find descriptions of water
        # water <- food_taxa %>% filter(FoodID != "92410110" &
        #                                 FoodID != "92410210"& 
        #                                 FoodID != "92410250"& 
        #                                 FoodID != "94000100"&
        #                                 FoodID != "94100100"&
        #                                 FoodID != "94100200" &
        #                                 FoodID != "94100300" &
        #                                 FoodID != "94210100" &
        #                                 FoodID != "94210200" &
        #                                 FoodID != "94210300" &
        #                                 FoodID != "94220200" &
        #                                 FoodID != "94300100")
        # 
        # #remove water
        # food_taxa <- food_taxa %>% filter(taxonomy %in% water$taxonomy)
        # 
        # # subset to remove water from the food table
        # food_map <- food_map[food_map$X.FOODID %in% food_taxa$Main.food.description,]
        # 
        # # clean up environment
        # rm(water)
        
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
        taxa <- read.delim(input$file_taxa$datapath, row =1)
        
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
        
        ###phylum level taxa ###
        ptaxa = taxa
        split <- strsplit(rownames(ptaxa),";")                               # Split and rejoin on lv7
        ptaxaStrings <- sapply(split,function(x) paste(x[1:2],collapse=";")) # level 7 is species, 8 is strain
        ptaxa <- rowsum(ptaxa,ptaxaStrings)                                    # Collapse by taxonomy name
        ptaxa = sweep(ptaxa,2,colSums(ptaxa),'/')
        ptaxa = ptaxa[order(rowMeans(ptaxa),decreasing=T),]
        
        #selet top four
        ptaxa <- ptaxa[1:4,]
        #traspose to add to map for later use
        ptaxa <- t(ptaxa)
        colnames(ptaxa) <- gsub(".*p__?", "", colnames(ptaxa))
        ptaxa <- as.data.frame(ptaxa)
        ptaxa <- rownames_to_column(ptaxa, "X.SampleID")
        
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
        
        myrownames <- rownames(gtaxa)
        
        gtaxa <- apply(gtaxa, 2, function(x) cut(x, breaks = seq(0, max(x), length.out = 15)))
        
        rownames(gtaxa) <- myrownames
        
        gtaxa <- as.data.frame(gtaxa)
        gtaxa <- rownames_to_column(gtaxa, "X.SampleID")
        
        
        ###Instantiate Plots###
        
        colnames(map)[2] <- "UserName" #Change col name from UserName.x to UserName - compatibility purposes
        colnames(map)[3] <- "StudyDayNo" #Change col name from StudyDayNo.x to StudyDayNo - compatibility purposes
        
        # add realative abundace of phylum taxa to the map for plotting
        map <- left_join(map,ptaxa)
        
        # add realative abundace of species taxa to the map for plotting
        map <- left_join(map,gtaxa, by = "X.SampleID")
        ###End Instantiation###
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        
        # for(id in 1:3){
        # 
        # render(input = "test.Rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document",
        #        output_dir = "../output/testput",envir = new.env(), 
        #        params = list(
        #          var1= metdf[1,3],
        #          var2= metdf[2,3]
        #        )) }
          
          #access params in test.rmd using the header
          #next stpe is to import the actual zebra code and incorporate it. test.rmd will hold mypdf.rmd
          #while app.R will try and instantiate everything. 
        
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
          taxalpha <- taxa[(colnames(taxa) %in% map[,"X.SampleID"])] #matrix utilized to calculate study average alpha div.
          
          #beta diversity iteration
          betataxa<- staxa[(colnames(staxa) %in% map[,"X.SampleID"])] #Create table with just taxa and subjects
          betataxa <- t(betataxa) #transpose 
          
          #blood draw values
          subblood <- UserName_map[UserName_map$UserName == id,]
          
          #render
          render(input = "bin/set_script.Rmd",output_file = paste0('report.', id, '.html'),"html_document",
                 output_dir = "output/zebra_output", envir = new.env()) }
          
          
        
      }
    )
  }
)

