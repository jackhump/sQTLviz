
library(shiny)
library(DT)
library(shinycssloaders)
library(shinyjs)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  tags$style(type="text/css", "
#UCSC_button {
  margin-top: 10px;
  margin-bottom: 10px;
  display: block;
  text-align: center;
}
#UCSC {
  border-color: #ccc;
  color: #333;
}
             
#UCSC:hover {
  background-color: #E6E6E6;
}

#titlePanel {
  text-align: center;
}

#logo {
    margin: 0 auto;
    width: 4%;
    position: fixed;
    z-index: -1;
}

#title {
    size:huge;
}

"),
             
  
  # Application title
  div(id = "titlePanel",
    #HTML("<img id=logo src=squirtle.png>"),
    h1("sQTLviz", id = "title")
    ),
    
    selectInput(inputId = "datasetChoice",
                label = "Which dataset?", 
                choices = list("Yang SNPs" = "YangResults",
                            "PD GWAS SNPs" = "GWASresults",
                            "SQTLs" = "resultsToPlot") ),
    
    fluidRow(
      column(8,offset=2,
        div(
          DT::dataTableOutput("all_clusters")
        )
      )
    ),
    withSpinner(div(
      div(
        plotOutput("select_gene_plot",width="100%", height = "200px") 
      ),
      fluidRow(
        column(7, offset = 1,
          div(
           plotOutput("select_cluster_plot", width = "100%", height = 400) 
          )
        ),
        column(3,
          div(
            plotOutput("select_box_plot", width = "100%", height = 400)
          )
        )
      ),
      div( id = "UCSC_button",
           htmlOutput("view_cluster_UCSC", inline = TRUE)
           )
    ), type = 8
    )
  )
)
