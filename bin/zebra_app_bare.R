 #Load packages required for the project
require(shiny)
require(rmarkdown)
require(knitr)
require(ggplot2)
require(reshape2)
require(tibble)
require(data.table)
require(vegan)
require(ape)
require(RColorBrewer)
require(dplyr)
require(cowplot)
require(kableExtra)

###Instantiation###
#set to working directory to project zebra directory
#setwd("/Users/abby/Documents/Projects/zebra/")
setwd("C:/Users/arzan/Documents/GitHub/zebra/")




shinyApp(
  ui = fluidPage(
    selectInput("beta_metric","Choose metric for beta diversity calculations:",
                c("bray", "euclidean","binomial","chao")),
    selectInput("alpha_metric","Choose metric for alpha diversity calculations:",
                c("shannon", "simpson","chao1")),
    textInput("idsample", "Enter Sample label", value = "X.SampleID", width = NULL, placeholder = "X.SampleID"),
    
    sidebarLayout(
      
      sidebarPanel(
       
        #Import mapping file
        fileInput("file_map", "Upload mapping file*",
                  buttonLabel = "Browse",
                  multiple = FALSE,
                  accept = ".txt",
                  placeholder = "Tab-delimited!"),
                  
        
        #Import nutrition table
        fileInput("file_nutr", "Upload nutrition table*",
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
        
        #Import username file table
        fileInput("file_umap", "Upload user mapping file",
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
                  placeholder = "Tab-delimited!"),
        
        #Import figure captions file
        fileInput("file_cap", "Upload captions for figures",
                  buttonLabel = "Browse",
                  multiple = FALSE,
                  accept = ".txt",
                  placeholder = "Place captions in an excel file")
        ),
    
      mainPanel(
        downloadButton("report", "Generate report")
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
        
        
        #Require inputs - switch on/off depending on what user can provideu
         req(input$file_map) 
         req(input$file_nutr) 
        # req(input$file_taxa) 
        # req(input$file_fmap) 
        # req(input$file_ftaxa) 
        
        #Input processing 
        
        map_sw   <- input$file_map$datapath
        user_sw  <- input$file_umap$datapath
        nutr_sw <-  input$file_nutr$datapath
        taxa_sw  <- input$file_taxa$datapath
        fmap_sw  <- input$file_fmap$datapath
        ftaxa_sw <- input$file_ftaxa$datapath
        idsample <- input$idsample
       
        
        #Metric/User choice set up:
        alpha_metric <- input$alpha_metric   #User input for alpha metric becomes value for alpha metric
        beta_metric  <- input$beta_metric    #User input for beta metric becomes value for beta metric
       
        #Caption set up
        
       # captions <- read.table(input$file_cap$datapath,sep = "\t", header = TRUE, comment = "")
        captions <- read.table(input$file_cap$datapath,sep = "\t", header = TRUE, comment = "")
       
        ##### Read data provided in to pertinent maps ####
   
        #read in mapping file
        SampleID_map <- read.table(input$file_map$datapath, sep = "\t", header = TRUE, comment = "")
       
        #Drop blanks from map
        SampleID_map <- SampleID_map[grep("Blank", SampleID_map[,idsample], invert = TRUE),]
        
        # create REE and TEE variables
        map <- SampleID_map %>% mutate(REE = ifelse(Gender=="M", (10*Weight+6.25*Height-5*Age+5), 
                                                    ifelse(Gender=="F", (10*Weight+6.25*Height-5*Age-161),
                                                           NA)), 
                                       TEE = REE*Activity.Factor) # calculate TEE with activity factor
        
        # make sure studydayno is a factor for later plotting
        map$RecordDayNo <- as.factor(map$RecordDayNo)
        
        
        ###UserName Map set up###
        if(is.null(user_sw)==FALSE){ #evaluator
        UserName_map <- read.table(input$file_umap$datapath, sep = "\t", header = TRUE, comment = "")
        } else {
          print("No Username file")
        }
        
        
        ####Nutrition Table Set up####
        if(is.null(nutr_sw)==FALSE){  #nutrition_table 
          
          nutrition_table <- read.table(input$file_nutr$datapath, sep = "\t", header = TRUE, comment = "")
          #Remove outlier nutr variables
          outliers = c("MCT.f.0021", "MCT.f.0044", "MCT.f.0050", "MCT.f.0056", "MCT.f.0058", "MCT.f.0060", "MCT.f.0076", "MCT.f.0116", "MCT.f.0465", "MCT.f.0486", "MCT.f.0601")
          nutrition_table <- nutrition_table[!(nutrition_table[,idsample] %in% outliers),]
        } else {
          print("Nutrition Table not provided. Pertinent nutrition analyses will not be included in the report.")
        }
        
        
        #### Food Map Set Up ####
        
        if(is.null(fmap_sw)==FALSE){ #food table evaluator
          
          food_map <-  read.table(file=input$file_fmap$datapath, sep = "\t", header = TRUE, comment = "")
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
          #boop
          merged_food<- merge(x=melted_food, y=map, by.x = "variable", by.y= idsample, all.x=TRUE)
          
          colnames(merged_food)[3] <- "UserName" #Change col name from UserName.x to UserName - compatibility purposes
          colnames(merged_food)[4] <- "RecordDayNo" #Change col name from StudyDayNo.x to StudyDayNo - compatibility purposes
        } else {
          print("Food mapping file not provided - report will exclude pertinent food data")
        }
        
        
        #### Taxa Table Set up ####
        
        if(is.null(taxa_sw)==FALSE){ #evaluator
          
          taxa <- read.delim(input$file_taxa$datapath, row =1)
          
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
          gtaxa <- rownames_to_column(gtaxa, idsample)
          
          ###Instantiate Plots (with gtaxa)###
          
          colnames(map)[2] <- "UserName" #Change col name from UserName.x to UserName - compatibility purposes
          colnames(map)[3] <- "RecordDayNo" #Change col name from StudyDayNo.x to StudyDayNo - compatibility purposes
          
          
          # add relative abundace of species taxa to the map for plotting
          map <- left_join(map,gtaxa, by = idsample)
        } else {
          print("Taxa table not provided. Pertinent taxa analyses will not be included in the report.")
        }
        
        for (id in unique(map$UserName)){
          
          #mapping
          submap <- map[map$UserName == id,] #create submap variable we can use to iterate through 
          
          #nutrition iteration
          if(is.null(nutr_sw)==FALSE){ #evaluator
            subnutr <- nutrition_table[nutrition_table$UserName == id,]
            macro_totals <-  round(mean(subnutr$PROT),2) +  round(mean(subnutr$TFAT),2) + round(mean(subnutr$CARB),2)
          } else (
            return(NULL)
          )
          
          #food iteration#
          if(is.null(fmap_sw)==FALSE){ #evaluator
            sub_foodmap <- merged_food[merged_food$UserName == id,] #sub food variable allows us to access one subject at a time (like submap)
            sub_food <- food_map[(colnames(food_map) %in% sub_foodmap[,"variable"])]#sub food variable looks at each subject individually
            sub_foodsp <- setDT(sub_food, keep.rownames = TRUE) #make rownames first column of frame
          } else {
            print("Mission Failed")
          }
          
          #taxa iteration#
          if(is.null(taxa_sw)==FALSE){ #evaluator
            subtaxa <- staxa[(colnames(staxa) %in% submap[,idsample])] #create subtaxa variable - looks at each subject individually
            subtaxa <- subtaxa[order(rowMeans(subtaxa), decreasing = TRUE),]# sort the subtaxa df by abundance
            subtaxasp <- setDT(subtaxa, keep.rownames=TRUE) #make rownames first column of frame
            #alpha diversity iteration
            subtaxaalpha <- taxa[(colnames(taxa) %in% submap[,idsample])] #subtaxa variable used for alpha diversity
            taxalpha <- taxa[(colnames(taxa) %in% map[,idsample])] #matrix utilized to calculate study average alpha div.
            #beta diversity iteration
            betataxa<- staxa[(colnames(staxa) %in% map[,idsample])] #Create table with just taxa and subjects
            betataxa <- t(betataxa) #transpose 
          } else {
            print("Mission Failed")
          }
          
          #blood iteration#
          
          if(is.null(user_sw)==FALSE){ #evaluator
          subblood <- UserName_map[UserName_map$UserName == id,] #blood draw values
          } else {
            print("No Username file")
          }
          
          #render
          render(input = "bin/zebra_markdown_bare.Rmd",output_file = paste0('report.', id, '.html'),"html_document",
                 output_dir = "output/demonstration_captions", envir = new.env()) }
          
          
        
      }
    )
  }
)

