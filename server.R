
library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(leafcutter)
library(reshape2)
library(gridExtra)
library(intervals) # needed for pretty strand arrow placement
library(foreach)
library(shinycssloaders)
library(grid)
library(gtable)
library(ggrepel)

if (!exists("introns")){
  load("sQTL_results.Rdata")
  defaultValue <- 1
}else{
  defaultValue <- NULL
}

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  output$all_clusters <- DT::renderDataTable({
    datatable( clusters[,c("gene","coord","N","FDR","annotation")],
               escape = FALSE,
               rownames = FALSE,
               colnames = c('Genomic location'='coord','Gene'='gene','N'='N','Annotation'='annotation','q'='FDR'),
               selection = 'single', 
               caption = "Click on a row to plot the corresponding visualization. N: number of introns within a cluster. q: Benjamini–Hochberg q-value.",
               fillContainer = FALSE,
               options = list(
                 pageLength = 15,
                 columnDefs = list(list(className = 'dt-center', targets = 0:4) )
               )
    ) 
  })
     

  
})
