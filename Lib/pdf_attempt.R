setwd("C:\\Users\\Surface Pro 3\\Documents\\GitHub\\zebra\\raw") # set to zebra directory

map <- read.table("new_map_with_treatment.txt", sep = "\t", header = TRUE, comment = "") #load mapping file

setwd("C:\\Users\\Surface Pro 3\\Documents\\GitHub\\zebra\\lib")
library("rmarkdown")
library("knitr")
for (id in unique(map$UserName)){
  subgroup <- map[map$UserName == id,]
  render("MCTS_pdf.rmd",output_file = paste0('report.', id, '.pdf'),"pdf_document")    
}
