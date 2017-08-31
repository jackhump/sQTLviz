
library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(leafcutter)
library(reshape2)
library(gridExtra)
# library(intervals) # needed for pretty strand arrow placement
# library(foreach)
library(shinycssloaders)
# library(grid)
# library(gtable)
# library(ggrepel)

if (!exists("introns")){
  load("../sQTL_results.Rdata")
  defaultValue <- 1
}else{
  defaultValue <- NULL
}

source("make_sQTL_cluster_plot.R")
make_sQTL_cluster_plot( "clu_34225",
                   main_title = "CAST",
                   vcf=vcf,
                   vcf_meta=vcf_meta,
                   exons_table = exons_table,
                   counts = clusters,
                   introns = annotatedClusters,
                   cluster_ids = annotatedClusters$clusterID,
                   snp_pos = "chr5:96076487",
                   snp = "rs7724759" )

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  # ALL CLUSTER X SNP TABLE
  output$all_clusters <- DT::renderDataTable({
    datatable( resultsToPlot,
                rownames = FALSE ) #,
               # #escape=FALSE,
               # #colnames = c('Genomic location'='coord','Gene'='gene','N'='N','Annotation'='annotation','q'='FDR'),
               # selection = 'single', 
               # caption = "Click on a row to plot the corresponding visualization.  q: Benjamini–Hochberg q-value.",
               # fillContainer = FALSE,
               # options = list(
               #   pageLength = 15,
               #   columnDefs = list(list(className = 'dt-center', targets = 0:5) )
               # )
#) 
  })
  
  # REACTIVITY
  
  values <- reactiveValues(default = defaultValue) # RBFOX1 in the Brain vs Heart dataset
  # REACTIVE VALUE IS UPDATED BY INPUT
  observeEvent(input$all_clusters_rows_selected,{
    #print("new row selected!")
    values$default <- input$all_clusters_rows_selected # if all_clusters_rows_selected changes then update value - this sets everything!
    #print(paste0("VALUE: ", values$default ))
  })
  
  # USE REACTIVE VALUE TO GENERATE ALL VARIABLES NEEDED
  
  mydata <- eventReactive(values$default,{
    sel <- values$default
    gene  <- resultsToPlot[ sel, ]$gene
    SNP <- resultsToPlot[ sel, ]$SNP
    SNP_pos <- resultsToPlot[ sel, ]$SNP_pos
    #gene <- gsub("<.*?>", "", gene) # strip out html italic tags
    #width <- getGeneLength(gene)
    clusterID <- row.names(resultsToPlot[ sel, ])
    cluster_pos <- clusters[ sel, ]$cluster_pos
    return(list(gene = gene, SNP=SNP, SNP_pos=SNP_pos, cluster_pos = cluster_pos, clusterID = clusterID) )
  })
  
  # PLOTTING
  
  output$select_cluster_plot <- renderPlot({
    plotTitle <- c(mydata()$gene, as.character(mydata()$cluster) )
    suppressWarnings( print(
      make_sQTL_cluster_plot( mydata()$cluster,
                         main_title = plotTitle,
                         vcf = vcf,
                         vcf_meta = vcf_meta,
                         exons_table = exons_table,
                         counts = clusters,
                         introns = annotatedClusters,
                         cluster_ids = annotatedClusters$clusterID,
                         snp_pos = mydata()$SNP_pos,
                         snp = mydata()$SNP )
    ))
  }, width = "auto", height = "auto",  res = 90
  )
  # make_sQTL_cluster_plot( "clu_34225",
  #                    main_title = "CAST",
  #                    vcf=vcf,
  #                    vcf_meta=vcf_meta,
  #                    exons_table = exons_table,
  #                    counts = clusters,
  #                    introns = annotatedClusters,
  #                    cluster_ids = annotatedClusters$clusterID,
  #                    snp_pos = "chr5:96076487",
  #                    snp = "rs7724759" )
  

  
})
