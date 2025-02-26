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
  function(input, output, session) {  
    ###########################################################################
    #filter data
    filteredData <- reactive({  
      data <- datasetInput()  
      if (input$Method != "All") {  
        data <- data[data$Method == input$Method, ]  
      }  
      if (input$Serotype != "All") {  
        data <- data[data$Serotype == input$Serotype, ]  
      }  
      if (input$Target != "All") {  
        data <- data[data$Target == input$Target, ]  
      } 
      return(data)  
    })  
    #primer input
    db_input <- reactive({
      data <- readDNAStringSet("Data/GISDD_All_20231204_Ge.fas","fasta")
      if (input$Serotype != "All") {  
        data <- readDNAStringSet("Data/GISDD_All_20231204_Ge.fas","fasta") 
      }  
      if (input$Serotype != "D1") {  
        data <- readDNAStringSet("Data/GISDD_D1_20231204_Ge.fas","fasta") 
      }  
      if (input$Serotype != "D2") {  
        data <- readDNAStringSet("Data/GISDD_D2_20231204_Ge.fas","fasta") 
      }  
      if (input$Serotype != "D3") {  
        data <- readDNAStringSet("Data/GISDD_D3_20231201_Ge.fas","fasta") 
      }  
      if (input$Serotype != "D4") {  
        data <- readDNAStringSet("Data/GISDD_D4_20231204_Ge.fas","fasta") 
      }  
	  return(data)  
    })
    # Output the filtered data as a table  
    output$table <- renderDataTable({  
      filteredData()  
    })  
    ###########################################################################
    #function
    datasetInput <- function() {  
      read_excel("Data/GISDDprimer.xlsx")  
    }  
    find_Gc <- function(input) {
      primer <- strsplit(as.character(input), "")
      rv_copy2 <- primer
      rv_copy <- primer
      for (m in 1:length(primer[[1]])) {
        # for max primer GC
        if (primer[[1]][m] %in% c("g", "R", "S", "K", "B", "D")) {
          rv_copy2[[1]][m] = "g"
        } else if (primer[[1]][m] %in% c("c", "Y", "M", "V", "H")) {
          rv_copy2[[1]][m] = "g"
        } else if (primer[[1]][m] %in% c("a", "W")) {
          rv_copy2[[1]][m] = "a"
        } else {
          rv_copy2[[1]][m] = "t"
        }
        # for min primer GC
        if (primer[[1]][m] %in% c("g", "S")) {
          rv_copy[[1]][m] = "g"
        } else if (primer[[1]][m] %in% c("t", "Y", "K", "B", "D", "H")) {
          rv_copy[[1]][m] = "a"
        } else if (primer[[1]][m] %in% c("a", "R", "M", "V")) {
          rv_copy[[1]][m] = "a"
        } else {
          rv_copy[[1]][m] = "c"
        }
      }
      max_Gc <-
        specify_decimal(calculate_GC(as.character(rv_copy2)), 2)
      min_Gc <-
        specify_decimal(calculate_GC(as.character(rv_copy)), 2)
      return(c(max_Gc, min_Gc))
    }
    specify_decimal <- function(x, k) {
      trimws(format(round(x, k), nsmall = k)) # make desired decimal in text part
    }
    calculate_GC <- function(input) {
      inp = strsplit(input, "")
      w = length(which(inp[[1]] == "a"))
      x = length(which(inp[[1]] == "t"))
      y = length(which(inp[[1]] == "g"))
      z = length(which(inp[[1]] == "c"))
      GC = (y + z) / (y + z + w + x) * 100
      return(GC)
    }
    matchrate <- function(primer){
      primer.str <- DNAString(primer)
      sequences <- readDNAStringSet("Data/db_demo.fas","fasta")
      max.mismatch = 3
      min.mismatch = 0
      matches <- vmatchPattern(primer,
                               sequences,
                               max.mismatch = max.mismatch, 
                               min.mismatch = min.mismatch,fixed = T)
      nmatch_per_seq <- elementNROWS(matches)  
      #nmatch_per_seq <- nmatch_per_seq[nmatch_per_seq > 0]  
      if(length(nmatch_per_seq) > 0){  
        return(sum(nmatch_per_seq) / length(nmatch_per_seq))  
      } else {  
        return(1)  
      }  
    }
    location <- function(primer){
      primer.str <- DNAString(primer)
      sequences <- db_input()
      max.mismatch = 3
      min.mismatch = 0
      matches <- vmatchPattern(primer,
                               sequences,
                               max.mismatch = max.mismatch, 
                               min.mismatch = min.mismatch,fixed = T)
      end <- paste(unique(end(matches)))
      as.numeric(min(end))-unique(matches@width0)+1

    }
    OfftargetSubg <- function(primer){
      primer.str <- DNAString(primer)
      sequences <- db_input()
      max.mismatch = 3
      min.mismatch = 0
      matches <- vmatchPattern(primer,
                               sequences,
                               max.mismatch = max.mismatch, 
                               min.mismatch = min.mismatch,fixed = T)
      if (length(matches) == 0) {  
        return("/")  
      }
      end <- paste(unique(end(matches)))
      Offtarget <- data.frame(isolate=matches@NAMES,offtarget=end)
      Offtarget <- subset(Offtarget,offtarget=="NA")
      Offtarget <- separate(Offtarget,isolate,into=c("Accession","label"),
                            sep=".1_")
      GISDD <- read_excel("Data/GISDD.v1.2.3.simple.xlsx")
      Offtarget <- left_join(Offtarget,GISDD)
      Offtarget <- unite(Offtarget,Virus_Type,Clade,col="Label",sep="_",remove=F)
      Label <-unique(Offtarget$Label) 
      Label <- toString(Label)
      print(Label) 
      if (length(print(as.character(Offtarget$isolate))) == 0) {  
        return("NA")  
      }
    }
    ###########################################################################
    #box
    output$NumberBox <- renderInfoBox({
      data <- datasetInput()
      infoBox(
        "Primer Number", 
        paste0(length(data$Sequence)),
        icon = icon("list"),
        color = "green"
        )
    })
    #box
    output$SetBox <- renderInfoBox({
      data <- datasetInput()
      infoBox(
        "Set Number", 
        paste0(length(unique(data$Set))),
        icon = icon("list"),
        color = "green"
      )
    })
    #box
    output$PaperBox <- renderInfoBox({
      data <- datasetInput()
      infoBox(
        "Involved Paper Number", 
        paste0(length(unique(data$Referrence))),
        icon = icon("list"),
        color = "green"
      )
    })
    ###########################################################################
    # Output the plot
    ###########################################################################
    # Output the plot  
    output$Plot1 <- renderPlot({  
      data <- filteredData()  
      ggplot(data) +  
        theme_classic() +  
        xlab("Position")+ylab("Primer number")+
        geom_histogram(aes(x = Position_Start),
                       color="black", fill="#51A261", linewidth=0.50) +  
        scale_y_continuous(expand = c(0,0),breaks = c(seq(0,100,by=5)))+
        scale_x_continuous(limits = c(0, 12000), expand = c(0, 0),  
                         breaks = seq(0, 12000, by = 1000))  +
        theme(legend.position = "none",
              axis.text = element_text(size=12),
              axis.title = element_text(size=15))
    })  
    # Output the plot  
    output$Plot2 <- renderPlot({  
      data <- filteredData()   
      data <- group_by(data, Target) %>% summarise(Number=n())
      data <- subset(data, Target!="NA")
      ggplot(data) +  
        theme_classic() +
        xlab("Target")+ylab("Primer number")+
        geom_col(aes(x = reorder(Target,Number)  , y = Number),
                 color="black",linewith=0.50,width=0.60,fill="#51A261") +  
        scale_y_continuous(expand = c(0,0),breaks = c(seq(0,100,by=20)))+
        scale_x_discrete(limits=c("5UTR","C/PrM","E","NS1","NS2A","NS2B","NS3","NS4A","NS4B","NS5","3UTR"))+
        theme(legend.position = "none",
              axis.text = element_text(size=12),
              axis.title = element_text(size=15))
    })  
    # Output the plot  
    output$Plot3 <- renderPlot({  
      data <- filteredData()   
      data <- group_by(data, Year) %>% summarise(Number=n())
      data <- subset(data, Year!="NA")
      ggplot(data) +  
        theme_classic() +
        xlab("Year")+ylab("Primer number")+
        geom_col(aes(x = Year  , y = Number),
                 color="black",linewith=0.50,width=0.60,fill="#51A261") +  
        scale_x_continuous(breaks = c(seq(0,3000,by=4)))+
        scale_y_continuous(expand = c(0,0))+
        theme(legend.position = "none",
              axis.text = element_text(size=12),
              axis.title = element_text(size=15))
    })  
    # Output the plot  
    output$Plot4 <- renderPlot({  
      data <- filteredData()   
      data <- group_by(data, Journal) %>% summarise(Number=n())
      data <- subset(data, Journal!="NA")
      ggplot(data) +  
        theme_classic() +
        xlab("Journal")+ylab("Primer number")+
        geom_col(aes(x = Journal  , y = Number),
                 color="black",linewith=0.50,width=0.60,fill="#51A261") +  
        theme(legend.position = "none",
              axis.text = element_text(size=12),
              axis.title = element_text(size=15))
    })  
    ###########################################################################
    # primer analysis
    output$PrimerOutTable <- renderDataTable({
      shiny::validate(shiny::need(input$upprimer, "Enter a valid sequence for your forward primer."))
      shiny::validate(shiny::need(input$downprimer, "Enter a valid sequence for your reverse primer."))
      shiny::validate(shiny::need(input$do3, "Click the 'Analyze' button, and wait about 20-30s."))
      namelist <- c(
        "Primers",
        "Primer Length",
        "Primer start",
        "Tm (GC content)",
        "Tm (nearest neighbor thermodynamics)",
        "Off-target ratio",
        "Off-target subgenotype"
      )
      upseq <- paste(input$upprimer)
      downseq <- paste(input$downprimer)
      downseq <- Biostrings::DNAString(downseq)
      downseq <- Biostrings::reverseComplement(downseq)
      downseq <- toString(downseq)
      Forward_primer <- c(
        upseq,
        nchar(upseq),
        location(upseq),
        Tm_GC(upseq,ambiguous=TRUE,variant="Primer3Plus",Na=50,mismatch=TRUE)[["Tm"]],
        Tm_NN(upseq,Na=50)[["Tm"]],
        paste(1-matchrate(upseq)) ,
        OfftargetSubg(upseq)#"5R"
      )
      Reverse_primer <- c(
        downseq,
        nchar(downseq),
        location(downseq),
        Tm_GC(downseq,ambiguous=TRUE,variant="Primer3Plus",Na=50,mismatch=TRUE)[["Tm"]],
        Tm_NN(downseq,Na=50)[["Tm"]],
        paste(1-matchrate(downseq)),
        OfftargetSubg(downseq)#"5R"
      )
      OutTable <- data.frame(Forward_primer = Forward_primer,
                 Reverse_primer = Reverse_primer)
      rownames(OutTable) <- namelist
      OutTable
    },filter="none")
    ###########################################################################
    # Download handler for the filtered data  
    output$downloadData <- downloadHandler(  
      filename = function() {  
        paste("filtered_data", ".csv", sep = "")  
      },  
      content = function(file) {  
        write.csv(filteredData(), file, row.names = FALSE)  
    })
    ###########################################################################
  }
