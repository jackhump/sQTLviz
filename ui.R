
library(shiny)
library(DT)
library(shinycssloaders)
library(shinyjs)

ui <- tagList(
  useShinyjs(),
  navbarPage(
    #title  = actionLink("aboutLink", "LeafViz"),
    title = "sQTLviz",
    #title = a(id = "gitLink", href="https://github.com/davidaknowles/leafcutter/tree/master/leafviz","LeafViz", target = "_blank"), 
    id = "navBarPage",
    windowTitle = "sQTLviz",
    tabPanel("CommonMind Consortium sQTLs", 
             value = "resultsToPlot",
             # padding-top: 70px
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

")),
    tabPanel("Li et al TWAS SNPs", value = "YangResults"),
    tabPanel("Nalls et al PD GWAS SNPs", value = "GWASresults" ),
    tabPanel("About", value = "NA"),
    
             
  
  # Application title
  div(id = "titlePanel",
    #HTML("<img id=logo src=squirtle.png>"),
    h1("sQTLviz: explore splicing QTLs", id = "title")
    ),
    
    
    fluidRow(
      column(8,offset=2,
        div(
          DT::dataTableOutput("all_clusters")
        )
      )
    ),
    withSpinner(div(
      column(7, offset = 1,
        h3(id = "title", "Whole gene visualization")
      ),
      div(
        plotOutput("select_gene_plot",width="100%", height = "350px") 
      ),
      fluidRow(
        column(7, offset = 1,
          div(
            h3(id = "title", "Cluster visualization"),
            h5(id = "subtitle", "Most significant junction is bolded")
          ),
          div(
           plotOutput("select_cluster_plot", width = "100%", height = "500px") 
          )
        ),
        column(3,
          h3(id = "title", "Most significant junction"),
          div(
            plotOutput("select_box_plot", width = "100%", height = "500px")
          )
        )
      ),
      fluidRow(
        column(8,offset=2,
          div(
        DT::dataTableOutput("junctionTable")
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
