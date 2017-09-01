
library(shiny)
library(DT)
library(shinycssloaders)
library(shinyjs)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("sQTLviz"),
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
      )
    ), type = 8
  )
))
