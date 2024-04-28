###############################################################################
###############################################################################
###############################################################################
  library(shiny)
  library(shinydashboard)
  library(ggplot2)
###############################################################################
  function(input, output, session) {  
    ###########################################################################
    datasetInput <- function() {  
      read.csv("Data/GISDDprimer.csv")  
    }  
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
    output$Plot1 <- renderPlot({  
      data <- filteredData()  
      ggplot(data) +  
        theme_bw() +  
        geom_histogram(aes(x = Position_Start)) +  
        scale_x_continuous(limits = c(0, 12000), expand = c(0, 0),  
                         breaks = seq(0, 12000, by = 1000))  
    })  
    # Output the plot  
    output$Plot2 <- renderPlot({  
      data <- filteredData()  
      data <- group_by(data, Target) %>% summarise(Number=n())
      data <- subset(data, Target!="NA")
      ggplot(data) +  
        theme_bw() +  
        geom_col(aes(x =Target, y = Number)) +  
        #scale_x_continuous(limits = c(0, 12000), expand = c(0, 0),  
         #                  breaks = seq(0, 12000, by = 1000))  
        theme()
    })  
    ###########################################################################
    # primer analysis
    output$PrimerOutTable <- renderDataTable({
      data.frame(up= paste(input$upprimer) )
    })
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
