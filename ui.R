
library(shiny)
library(DT)
library(shinycssloaders)
library(shinyjs)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("sQTLviz"),
  
  # Sidebar with a slider input for number of bins 
  # sidebarLayout(
  #   sidebarPanel(
      div(
        withSpinner(DT::dataTableOutput("all_clusters"))
      ),
      div(
        withSpinner(plotOutput("select_cluster_plot", width = "100%") )
      )
    #),
    
    # Show a plot of the generated distribution
    #mainPanel(
    #   plotOutput("distPlot")
    #)
  #)
))
