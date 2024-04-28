
  library(ggplot2)
  
  data <- read.csv("Data/GISDDprimer.csv")  
  
  
  ggplot(GISDDprimer)+
    geom_point(aes(Position_start, Primer))
  
  
  
  
  
  
  
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
  
  
  
  
  ,
  column(width = 12,
         box(tabItem(tabName = "downloadData",  
                     fluidRow(  
                       column(12,  
                              actionButton("downloadBtn", "Download Filtered Data")  
                       )  
                     ) 