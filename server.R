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
library(ggbeeswarm)
library(stringr)

if (!exists("introns")){
  load("sQTL_results.Rdata")
  defaultValue <- 1
}else{
  defaultValue <- NULL
}

source("make_sQTL_cluster_plot.R")
source("make_sQTL_gene_plot.R")
#source("/Users/Jack/google_drive/Work/PhD_Year_3/leafcutter/leafcutter/R/make_gene_plot.R")
source("make_sQTL_box_plot.R")
# sel <- 5
# junction_to_plot <- sigJunctions[ sigJunctions$clu == row.names(resultsToPlot)[sel], ]
# sigJunction <- junction_to_plot[ which( junction_to_plot$bpval == min(junction_to_plot$bpval) ), ]$pid
# make_sQTL_cluster_plot( row.names(resultsToPlot)[sel],
#                    main_title = "test",
#                    vcf=vcf,
#                    vcf_meta=vcf_meta,
#                    exons_table = exons_table,
#                    counts = clusters,
#                    introns = annotatedClusters,
#                    cluster_ids = annotatedClusters$clusterID,
#                    snp_pos = resultsToPlot[sel,]$SNP_pos,
#                    snp = resultsToPlot[sel,]$SNP,
#                    sigJunction = sigJunction )

# current issue - SNP coord is 0 - how wide to make it?

#source("/Users/Jack/google_drive/Work/PhD_Year_3/leafcutter/leafcutter/R/make_gene_plot.R")
#source("make_sQTL_gene_plot.R")
# sel <- 2
# make_gene_plot(resultsToPlot[sel,]$gene,
#                counts = clusters,
#                introns = annotatedClusters,
#                exons_table = exons_table,
#                cluster_list = NULL,
#                clusterID = row.names(resultsToPlot)[sel],
#                introns_to_plot = introns_to_plot,
#                snp_pos = resultsToPlot[sel,]$SNP_pos,
#                snp = resultsToPlot[sel,]$SNP
# )

# 
# sel <- 19
# junction_to_plot <- sigJunctions[ sigJunctions$clu == row.names(resultsToPlot)[sel], ]
# junction_to_plot <- junction_to_plot[ which( junction_to_plot$bpval == min(junction_to_plot$bpval) ), ]$pid
# make_sQTL_box_plot(
#     cluster_to_plot =  row.names(resultsToPlot)[sel],
#     junction_to_plot = junction_to_plot,
#     all_junctions = all_junctions,
#     main_title = NA,
#     vcf = vcf,
#     vcf_meta = vcf_meta,
#     exons_table = exons_table,
#     counts = clusters,
#     introns = annotatedClusters,
#     cluster_ids = annotatedClusters$clusterID,
#     snp_pos = resultsToPlot[sel,]$SNP_pos,
#     snp = resultsToPlot[sel,]$SNP )
# 
# # Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Three different selections:
  # 1) PD GWAS SNPs
  # 2) Yang's TWAS SNPs
  # 3) All SNPs in the sQTL analysis
  
  #output$resultsTable <-
  
  #output$test_text <- renderText({paste0("You are viewing tab \"", input$navBarPage, "\"")})
  
  
  # ALL CLUSTER X SNP TABLE
  output$all_clusters <- DT::renderDataTable({
    #subsetChoice <- eval(parse(text=input$datasetChoice) )
    # clicked tab gives you the thing to subset
    subsetChoice <- eval(parse(text=input$navBarPage ) )
    df <- subset(resultsToPlot, row.names(resultsToPlot) %in% row.names(subsetChoice) )
    if( input$navBarPage == "GWASresults"){
      df <- df %>%
        mutate("GWAS P" = GWASresults$`GWAS P` ) %>%
        select("GWAS P", everything() )
    }
    datatable( df,
                rownames = FALSE,
               # #escape=FALSE,
               # #colnames = c('Genomic location'='coord','Gene'='gene','N'='N','Annotation'='annotation','q'='FDR'),
                selection = 'single',
               # caption = "Click on a row to plot the corresponding visualization.  q: Benjamini–Hochberg q-value.",
                fillContainer = FALSE ,
                options = list(
                  #TODO: why doesn't this work?
                  language = list(
                    searchPlaceholder = "for a SNP or gene..."
                  )
                  #   pageLength = 15,
               #   columnDefs = list(list(className = 'dt-center', targets = 0:5) )
                )
    )
#) 
  })
  
  # JUNCTION TABLE
  output$junctionTable <- DT::renderDataTable({
    jtable <- dplyr::filter(junctionTable, clu == mydata()$clusterID) %>%
      dplyr::select(-clu)
    datatable(jtable,
              escape = FALSE,
              rownames = FALSE,
              fillContainer = FALSE,
              options <- list( searching = FALSE, paging = FALSE, info = FALSE ))
  })
  
  # REACTIVITY
  
  values <- reactiveValues(default = defaultValue) 
  # REACTIVE VALUE IS UPDATED BY INPUT
  observeEvent(input$all_clusters_rows_selected,{
    print("new row selected!")
    values$default <- input$all_clusters_rows_selected # if all_clusters_rows_selected changes then update value - this sets everything!
    print(paste0("VALUE: ", values$default ))
  })
  
  # SECOND REACTIVE VALUE - HOW TO SUBSET TABLE
  
  
  # USE REACTIVE VALUE TO GENERATE ALL VARIABLES NEEDED
  
  mydata <- eventReactive(values$default,{
    subsetChoice <- eval(parse(text=input$navBarPage) )
    df <- subset(resultsToPlot, row.names(resultsToPlot) %in% row.names(subsetChoice) )
    sel <- values$default
    print(sel)
    gene  <- df[ sel, ]$gene
    SNP <- df[ sel, ]$SNP
    SNP_pos <- df[ sel, ]$SNP_pos
    #gene <- gsub("<.*?>", "", gene) # strip out html italic tags
    #width <- getGeneLength(gene)
    clusterID <- row.names(df)[sel]
    print(clusterID)
    cluster_pos <- df[ sel, ]$cluster_pos
    # get the most significant junction in the selected cluster
    junction <- sigJunctions[ sigJunctions$clu == row.names(df)[sel], ]
    junction <- junction[ which( junction$bpval == min(junction$bpval) ), ]$pid # causing problems - what is pid?
    return(list(gene = gene, SNP=SNP, SNP_pos=SNP_pos, cluster_pos = cluster_pos, clusterID = clusterID, width = "auto", junction=junction) )
  })
  
  # PLOTTING
  
  output$select_cluster_plot <- renderPlot({
    #plotTitle <- c(mydata()$gene, as.character(mydata()$clusterID) )
    suppressWarnings( print(
      make_sQTL_cluster_plot(
                         cluster_to_plot =  mydata()$clusterID,
                         main_title = NA,
                         vcf = vcf,
                         vcf_meta = vcf_meta,
                         exons_table = exons_table,
                         counts = clusters,
                         introns = annotatedClusters,
                         cluster_ids = annotatedClusters$clusterID,
                         snp_pos = mydata()$SNP_pos,
                         snp = mydata()$SNP,
                         sigJunction = mydata()$junction)
    ))
  }, width = "auto", height = "auto",  res = 60
  )
  
  # WHOLE GENE PLOTTING
  observeEvent(values$default,{ 
    output$select_gene_plot <- renderPlot({
      suppressWarnings( print( 
        make_gene_plot(mydata()$gene,
                       counts = clusters,
                       introns = annotatedClusters,
                       exons_table = exons_table,
                       cluster_list = NULL,
                       clusterID = mydata()$clusterID,
                       introns_to_plot = introns_to_plot,
                       snp_pos = mydata()$SNP_pos,
                       snp = mydata()$SNP )
      )
      )
    }, width = mydata()$width, height = "auto", res = 90 # try changing height param
    )
  })

  # BOX PLOTS OF GENOTYPE AGAINST NORMALISED JUNCTION COUNTS
  
  output$select_box_plot <- renderPlot({
    #plotTitle <- c(mydata()$gene, as.character(mydata()$clusterID) )
    suppressWarnings( print(
      make_sQTL_box_plot(
        cluster_to_plot =  mydata()$clusterID,
        junction_to_plot = mydata()$junction,
        all_junctions = all_junctions,
        main_title = NA,
        vcf = vcf,
        vcf_meta = vcf_meta,
        exons_table = exons_table,
        counts = clusters,
        introns = annotatedClusters,
        cluster_ids = annotatedClusters$clusterID,
        junctionTable = junctionTable,
        snp_pos = mydata()$SNP_pos,
        snp = mydata()$SNP )
    ))
  }, width = "auto", height = "auto",  res = 90
  )
  
  # VIEW CLUSTER IN UCSC
  
  output$view_cluster_UCSC <- renderUI({
    coord <- mydata()$cluster_pos
    print("coord:")
    print(coord)
    snp <- mydata()$SNP
    url <- paste0( "http://genome.ucsc.edu/cgi-bin/hgTracks?&org=human&db=hg19&position=", 
                   coord,"&hgFind.matches=", snp  )
  return(tags$a(href = url, "view on UCSC", target = "_blank", class = "btn btn_default", id = "UCSC" ) )
  })
  
})
