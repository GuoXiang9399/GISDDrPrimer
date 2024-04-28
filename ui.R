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
###############################################################################
# Define UI for application that draws a histogram
  dashboardPage(
    skin = "yellow",
    dashboardHeader(title = "GISDDprimer"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Dashborad",tabName = "PrimerDashborad",icon = icon("dashboard")),
        menuItem("Primer Bank",tabName = "PrimerBank", icon = icon("th")),
        menuItem("Primer Evaluation",tabName = "PrimerEvaluation",icon = icon("virus")),
        menuItem("About",tabName = "About", icon = icon("cube"))
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
                    dataTableOutput("table"),width = 12,
                    style = "height:800px; overflow-y: scroll;overflow-x: scroll;"
                    )
                ),
        tabItem(tabName = "PrimerEvaluation",
                h2("Primer Evaluation"),
                fluidRow(
                  box(textInput("upprimer","Forward primer:")),
                  box(textInput("downprimer","Reverse primer:"))),
                box(
                  dataTableOutput("PrimerOutTable"),width = 12
                )
                ),
        tabItem(tabName = "About",
                h2("GISDDrPrimer")
                )
        )
      )
    )
       
