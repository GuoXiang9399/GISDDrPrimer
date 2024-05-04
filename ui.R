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
        menuItem("Primer Evaluation",tabName = "PrimerEvaluation",icon = icon("virus"))
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
                  box(textInput("upprimer","Your Forward Primer (5'->3'):")),
                  box(textInput("downprimer","Your Reverse Primer (5'->3'):"))),
                box(
                  selectInput("Serotype", "Choose a serotype:",
                              c("All","Universal", "D1", "D2","D3","D4")),
                  dataTableOutput("PrimerOutTable"),width = 12
                  )
                ),
        tabItem(tabName = "About",
                h2("GISDDrPrimer"),
                h4("GISDDrPrimer is a Shiny R based application that allows searching for published/novel primer pairs targeting the DENV fragment and visualizing their alignment to the reference genome."),
                h4("GISDDrPrimer shows the amplicon and any variation found in the amplicon as well as in the primer binding regions.  GISDDrPrimer also provides a list of in-house designed conserved and degenerate primer pairs across the viral genome and presents information on occurrence and lineage of mutations with the degenerate codes, alignment to the genome.")
                )
        )
      )
    )
       
