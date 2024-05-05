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
    dashboardHeader(title = "GISDDprimer",
                    dropdownMenu(type = "notifications",
                                 notificationItem( text = "Last update of data: 2022-09-14",icon = icon("calendar"))),
                    dropdownMenu(type = "messages",icon = icon("at"),
                                 messageItem(from = "Contact information",message = tags$div("guoxiang19939@163.com" , tags$br(),"601765739@qq.com",style = "display: inline-block; vertical-align: middle;")))),
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
                  infoBoxOutput("SetBox"),
                  infoBoxOutput("PaperBox")
                ),
                fluidRow(
                  box(title = "TOP Affiliations",  
                      background = "yellow",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      imageOutput("Plot1", height = 350)),
                  box(title = "TOP Affiliations",  
                      background = "yellow",
                      solidHeader = TRUE,
                      collapsible = TRUE,
                      imageOutput("Plot2", height = 350))
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
                              c("All","Universal", "D1", "D2","D3","D4"))),
                  box(actionButton("do3", h4(strong("Search")), icon("search-plus")))),
                box(
                  dataTableOutput("PrimerOutTable"),width = 12  )
                ),
        tabItem(tabName = "Links",
                h2("Primer tools"),
                fluidRow(
                  box(solidHeader = TRUE,
                      collapsible = TRUE,
                      h2("Primer3web"),
                      h4("https://primer3.ut.ee/"),
                      h5("an integrating masking of template sequence with primer design software.")),
                  box(solidHeader = TRUE,
                      collapsible = TRUE,
                      h2("MFEprimer"),
                      h4("https://www.mfeprimer.com/"),
                      h5("The MFEprimer Web Server version 4.0 serves as a comprehensive platform for PCR primer quality control and design tools."))),
                fluidRow(
                  box(solidHeader = TRUE,
                      collapsible = TRUE,
                      h2("primalscheme"),
                      h4("https://primalscheme.com/"),
                      h5("primalscheme is a tool for designing primer panels for multiplex PCR. It uses a greedy algorithm to find primers for tiling amplicon generation for multiple reference genomes. It works best on viral isolates of known lineages e.g. outbreak strains.")),
                  box(solidHeader = TRUE,
                      collapsible = TRUE,
                      h2("Oligo"),
                      h4("https://www.oligo.net/"),
                      h5("OLIGO Primer Analysis Software is the essential tool for designing and analyzing sequencing and PCR primers, synthetic genes, and various kinds of probes including siRNA and molecular beacons."))),
                fluidRow(
                  box(solidHeader = TRUE,
                      collapsible = TRUE,
                      h2("qPrimerDB"),
                      h4("https://qprimerdb.biodb.org"),
                      h5("The qPrimerDB is the most comprehensive qPCR primer database available to date, with a web front-end providing gene-specific and pre-computed primer pairs over 1000 important organisms, including human, mouse, zebrafish, rice, maize, fungi, and microsporidia.")),
                  box(solidHeader = TRUE,
                      collapsible = TRUE,
                      h2("Primer Biosoft"),
                      h4("http://www.premierbiosoft.com/"),
                      h5("a group of computer scientists and biologists dedicated to accelerating research in life sciences.")
                      ))
                ),
        tabItem(tabName = "About",
                h2("GISDDrPrimer"),
                h4("GISDDrPrimer is a Shiny R based application that allows searching for published/novel primer pairs targeting the DENV fragment and visualizing their alignment to the reference genome."),
                #imageOutput("D1Genome"),
                h4("GISDDrPrimer shows the amplicon and any variation found in the amplicon as well as in the primer binding regions.  GISDDrPrimer also provides a list of in-house designed conserved and degenerate primer pairs across the viral genome and presents information on occurrence and lineage of mutations with the degenerate codes, alignment to the genome.")
                )
        )
      )
    )
       
