
  library(ggplot2)
  
  data <- read.csv("Data/GISDDprimer.csv")  
  
  
  ggplot(GISDDprimer)+
    geom_point(aes(Position_start, Primer))
  
  
  
  
  