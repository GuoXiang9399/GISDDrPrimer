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
    matchrate_D1 <- function(primer){
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
    location_D1 <- function(primer){
      primer.str <- DNAString(primer)
      sequences <- readDNAStringSet("Data/db_demo.fas","fasta")
      max.mismatch = 3
      min.mismatch = 0
      matches <- vmatchPattern(primer,
                               sequences,
                               max.mismatch = max.mismatch, 
                               min.mismatch = min.mismatch,fixed = T)
      end <- paste(unique(end(matches)))
      as.numeric(min(end))-unique(matches@width0)+1
    }
    range_st_D1 <- function(input) {
      D1 <- "ACCACAAAGAAAGACTTAGGGATTGGCCATGTGGCTGTTGAAAACCACCACCATGCCGCAATGCTGGACGTGGACTTACATCCAGCTTCAGCCTGGACCCTCTATGCAGTGGCCACAACAATTATCACTCCCATGATGAGGCACACAATCGAAAACACAACGGCAAACATTTCCCTGACAGCCATTGCAAACCAGGCAGCTATATTGATGGGACTTGACAAAGGATGGCCAATATCAAAGATGGACATAGGAGTTCCACTTCTCGCCTTGGGGTGCTATTCCCAGGTGAATCCACTGACGCTGACAGCGGCGGTATTGATGCTAGTGGCTCATTACGCTATAATTGGACCTGGACTGCAAGCAAAAGCTACTAGAGAAGCTCAAAAAAGGACAGCGGCCGGAATAATGAAAAATCCAACCGTTGATGGAATCGTTGCAATAGATTTGGACCCTGTGGTTTATGATGCGAAATTTGAGAAACAACTAGGCCAAATAATGTTGCTGATACTATGCACATCACAGATCCTCTTGATGCGGACTACATGGGCCTTGTGCGAATCCATCACGTTGGCCACTGGACCTCTGACCACGCTCTGGGAGGGATCTCCAGGAAAATTTTGGAACACCACGATAGCGGTTTCCATGGCAAACATTTTCAGAGGAAGTTATCTAGCAGGAGCAGGTCTGGCCTTCTCATTAATGAAATCTCTAGGAGGAGGTAGGAGAGGCACGGGAGCCCAAGGGGAAACACTGGGAGAGAAATGGAAAAGACGACTGAACCAACTGAGCAAGTCAGAATTTAACACCTATAAAAGGAGTGGGATTATGGAAGTGGACAGATCCGAAGCCAAAGAGGGACTGAAAAGAGGAGAAACAACCAAACATGCAGTGTCGAGAGGAACCGCTAAATTGAGGTGGTTTGTGGAGAGGAACCTTGTGAAACCAGAAGGGAAAGTCATAGACCTCGGTTGTGGAAGAGGTGGCTGGTCATATTATTGCGCTGGGCTGAAGAAAGTCACAGAAGTGAAGGGATATACAAAAGGAGGACCTGGACATGAAGAACCAATCCCAATGGCGACCTATGGATGGAACCTAGTAAAGCTGCATTCCGGGAAAGACGTATTCTTTATACCACCTGAGAAATGTGACACCCTTTTGTGTGATATTGGTGAGTCCTCTCCAAACCCAACTATAGAGGAAGGAAGAACGCTACGCGTCCTAAAGATGGTGGAACCATGGCTCAGAGGAAACCAATTTTGCATAAAAATTCTGAATCCCTACATGCCAAGTGTGGTGGAAACTCTGGAGCAAATGCAAAGAAAACATGGAGGGATGCTAGTGCGAAATCCACTTTCAAGAAATTCCACTCATGAAATGTATTGGGTTTCATGTGGAACAGGAAACATTGTGTCAGCAGTAAACATGACATCCAGAATGTTGCTAAATCGATTCACAATGGCTCACAGGAAACCAACATATGAAAGAGACGTGGACCTAGGCGCCGGAACAAGACACGTGGCAGTGGAACCAGAGGTAGCCAACCTAGATATCATTGGCCAGAGGATAGAGAACATAAAACATGAACACAAGTCAACATGGCATTATGATGAGGACAATCCATACAAAACATGGGCCTATCATGGATCATATGAGGTCAAGCCATCAGGATCAGCCTCATCCATGGTCAATGGCGTGGTGAAACTGCTCACCAAACCATGGGATGTCATCCCCATGGTCACACAA"  
      seq_up <- paste(input)  
      alignment <- pairwiseAlignment(D1, seq_up, type="global") 
      alignment@pattern@range@start
    }
    ###########################################################################
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
    # Output the filtered data as a table  
    output$table <- renderDataTable({  
      filteredData()  
    })  
    ###########################################################################
    #box
    output$NumberBox <- renderInfoBox({
      data <- datasetInput()
      infoBox(
        "Primer Number", 
        paste0(length(data$Sequence)),
        icon = icon("list"),
        color = "yellow"
        )
    })
    #box
    output$SetBox <- renderInfoBox({
      data <- datasetInput()
      infoBox(
        "Set Number", 
        paste0(length(unique(data$Set))),
        icon = icon("list"),
        color = "yellow"
      )
    })
    #box
    output$PaperBox <- renderInfoBox({
      data <- datasetInput()
      infoBox(
        "Involved Paper Number", 
        paste0(length(unique(data$Referrence))),
        icon = icon("list"),
        color = "yellow"
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
                       color="black", fill="#9BD5E7", linewidth=0.50) +  
        scale_y_continuous(expand = c(0,0),breaks = c(seq(0,100,by=1)))+
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
        geom_col(aes(x = reorder(Target,-Number)  , y = Number,fill=Target),
                 color="black",linewith=0.5,width=0.65) +  
        scale_y_continuous(expand = c(0,0),breaks = c(seq(0,100,by=1)))+
        theme(legend.position = "none",
              axis.text = element_text(size=12),
              axis.title = element_text(size=15))
    })  
    ###########################################################################
    # primer analysis
    output$PrimerOutTable <- renderDataTable({
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
      #downseq <- Biostrings::DNAString(downseq)
      #downseq <- Biostrings::reverseComplement(downseq)
      Forward_primer <- c(
        upseq,
        nchar(upseq),
        location_D1(upseq),
        Tm_GC(upseq,ambiguous=TRUE,variant="Primer3Plus",Na=50,mismatch=TRUE)[["Tm"]],
        Tm_NN(upseq,Na=50)[["Tm"]],
        paste(1-matchrate_D1(upseq)) ,
        "5R"
      )
      Reverse_primer <- c(
        downseq,
        nchar(downseq),
        location_D1(downseq),
        Tm_GC(downseq,ambiguous=TRUE,variant="Primer3Plus",Na=50,mismatch=TRUE)[["Tm"]],
        Tm_NN(downseq,Na=50)[["Tm"]],
        paste(1-matchrate_D1(downseq)),
        "5R"
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
