###############################################################################
###############################################################################
###############################################################################
  library(shiny)
  library(shinydashboard)
  library(ggplot2)
###############################################################################
  function(input, output, session) {  
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
      if (input$Set != "All") {  
        data <- data[data$Target == input$Set, ]  
      }  
      return(data)  
    })  
    # Output the filtered data as a table  
    output$table <- renderDataTable({  
      filteredData()  
    })  
  # Output the plot  
  output$Plot <- renderPlot({  
    data <- filteredData()  
    
    ggplot(data) +  
      theme_bw() +  
      geom_point(aes(x = Position_start, y = Primer)) +  
      geom_point(aes(x = Position_end, y = Primer)) +  
      scale_x_continuous(limits = c(0, 12000), expand = c(0, 0),  
                         breaks = seq(0, 12000, by = 1000))  
  })  
  
  # Download handler for the filtered data  
  output$downloadData <- downloadHandler(  
    filename = function() {  
      paste("filtered_data", ".csv", sep = "")  
    },  
    content = function(file) {  
      write.csv(filteredData(), file, row.names = FALSE)  
    }  
  )  
}
