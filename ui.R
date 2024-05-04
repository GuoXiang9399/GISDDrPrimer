###############################################################################
###############################################################################
###############################################################################
  library(shiny)
  library(shinydashboard)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(DT)
  library(ggplot2)
  library(viridis)
  library(TmCalculator)
  library(Biostrings)
###############################################################################
# Define UI for application that draws a histogram
  dashboardPage(
    skin = "yellow",
    dashboardHeader(title = "GISDDprimer"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("About",tabName = "About", icon = icon("cube")),
        menuItem("Dashborad",tabName = "PrimerDashborad",icon = icon("dashboard")),
        menuItem("Primer Bank",tabName = "PrimerBank", icon = icon("th")),
        menuItem("Primer Evaluation",tabName = "PrimerEvaluation",icon = icon("virus")),
        menuItem("Links", tabName = "Links", icon = icon("cube"))
      )),
    dashboardBody(
      tabItems(
        tabItem(tabName = "PrimerDashborad",
                h2("Dashborad for DENV primers"),
                fluidRow(
                  infoBoxOutput("NumberBox"),
                  infoBoxOutput("PaperBox")
                ),
                fluidRow(
                  box(title = "TOP Affiliations",  
                      background = "yellow",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      imageOutput("Plot1", height = 420)),
                  box(title = "TOP Affiliations",  
                      background = "yellow",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      imageOutput("Plot2", height = 420))
                )),
        tabItem(tabName = "PrimerBank",
                box(selectInput("Method", "Choose a method:",
                              c("All","Detection","Sanger sequence","High through sequence")),
                    selectInput("Serotype", "Choose a serotype:",
                              c("All","Universal", "D1", "D2","D3","D4")),
                    selectInput("Target", "Choose a target:",
                              c("All","C-prM","E","whole genome","Others")),
                    dataTableOutput("table"),
                    downloadButton("downloadData", "Download"),
                    width = 12, style = "height:800px; overflow-y: scroll;overflow-x: scroll;"
                    ),
                ),
        tabItem(tabName = "PrimerEvaluation",
                h2("Primer Evaluation"),
                fluidRow(
                  box(textInput("upprimer","Your Forward Primer (5'->3'):"),
                      h5("example primer from Lanciotti_1992: ",
                         HTML("<span style='color: blue; font-weight: bold;'>TCAATATGCTGAAACGCGCGAGAAACCG</span>")) ),
                  box(textInput("downprimer","Your Reverse Primer (5'->3'):"),
                      h5("example primer from Lanciotti_1992: ",
                         HTML("<span style='color: blue; font-weight: bold;'>TTGCACCAACAGTCAATGTCTTCAGGTTC</span>")))),
                fluidRow(
                  box(selectInput("Serotype", "Choose a serotype:",
                              c("All","Universal", "D1", "D2","D3","D4")))),
                box(
                  dataTableOutput("PrimerOutTable"),width = 12  )
                ),
        tabItem(tabName = "Links",
                h2(HTML("www.baidu.com")),
                fluidRow(
                  box(solidHeader = TRUE,
                      collapsible = TRUE,
                      h2("primalscheme"),
                      h5("https://primalscheme.com/")),
                  box(solidHeader = TRUE,
                      collapsible = TRUE,
                      h2("primalscheme"))
                ),
                
                h5("Our study has established a reproduceable and comparable global genotyping framework of DENV with contextualizing spatiotemporal epidemiological information before. The defned framework was discriminated with three hierarchical layers of genotype, subgenotype and clade with respective mean pairwise distances 2–6%, 0.8–2%, and ≤0.8%. This framework reveals that the persisting traditional endemic sourcing, the emerging epidemic difusing, and the probably hidden epidemics are the crucial drivers of the rapid global spread of dengue. ")
        ),
        tabItem(tabName = "About",
                h2("GISDDrPrimer"),
                h4("GISDDrPrimer is a Shiny R based application that allows searching for published/novel primer pairs targeting the DENV fragment and visualizing their alignment to the reference genome."),
                h4("GISDDrPrimer shows the amplicon and any variation found in the amplicon as well as in the primer binding regions.  GISDDrPrimer also provides a list of in-house designed conserved and degenerate primer pairs across the viral genome and presents information on occurrence and lineage of mutations with the degenerate codes, alignment to the genome.")
                )
        )
      )
    )
       
