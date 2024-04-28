###############################################################################
###############################################################################
###############################################################################
  library(shiny)
  library(shinydashboard)
###############################################################################
# Define UI for application that draws a histogram
  dashboardPage(
    dashboardHeader(
      title = "GISDDprimer"
    ),
    dashboardSidebar(
      selectInput("Method", "Choose a method:",
                  c("All",
                    "Detection", 
                    "Sanger sequence", 
                    "High through sequence")),
      selectInput("Serotype", "Choose a serotype:",
                  c("All",
                    "Universal", 
                    "D1", 
                    "D2",
                    "D3",
                    "D4")),
      selectInput("Target", "Choose a target:",
                  c("All",
                    "C-prM", 
                    "E", 
                    "whole genome", 
                    "Others")),
      selectInput("Set", "Choose a target:",
                  c("All")),
      menuItem("Download", tabName = "downloadData", icon = icon("download"))  
    ),
    dashboardBody(
      fluidRow(
        column(width = 12,
          box(dataTableOutput("table"), width = NULL)),
        column(width = 12,
          box(tabItem(tabName = "downloadData",  
              fluidRow(  
                column(12,  
                       actionButton("downloadBtn", "Download Filtered Data")  
                )  
              ) 
      ))     
               ),
        
        column(width = 12,
          box(plotOutput("Plot"), width = NULL)),
        ),
      
      )
    )
    
