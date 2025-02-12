# library(ggplot2)
# 
# InputTab <- readRDS(file = "segundoData/segundoData.rds")[[1]]
# selVar <- "OS_DAYS"
# infoTab <-  read.table(file = "abbreviationList.csv",header = F,sep = ";")
# 
# asd <- basicPlotter(selVar = "Age",InputTab = readRDS(file = "cuartoData/CuartoData.rds")[[1]])
# 
# asd$sumPlot
# asd$sumTab


basicPlotter <- function(selVar,InputTab, infoTab){
  redTab <- InputTab[,colnames(InputTab) == selVar]
  
  if (class(redTab) %in% c("character","factor")) {
    ploTab <- as.data.frame(table(redTab))
    colnames(ploTab) <- c("Var", "N")
    sumPlot <- ggplot(ploTab, aes(x=Var, y=N)) +
      geom_segment( aes(x=Var, xend=Var, y=0, yend=N), color="grey") +
      geom_point( color="orange", size=4) +
      theme_light() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      xlab("") +
      ylab("Number of patients")
    
    sumTab <- t(data.frame(unclass(table(redTab)), check.names = FALSE))
    
  } else if (class(redTab) == "numeric") {
    
    redTab <- InputTab[,colnames(InputTab) == selVar]
    redTab <- data.frame("Var" = redTab)
    
    ploTab <- as.data.frame(table(redTab))
    colnames(ploTab) <- c("Var", "N")
    ploTab$Var <- as.numeric(ploTab$Var)
    
    sumPlot <-  ggplot(redTab, aes(x=Var)) +
      geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
      ggtitle(paste0("Distribution of patients according to ",selVar))+
      theme_light() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank()) +
      xlab(paste0(selVar)) +
      ylab("Number of patients")
    
    sumTab <- t(data.frame(unclass(summary(redTab$Var)), check.names = FALSE))
    
  }
  
  if (selVar %in% infoTab$V1) {
    dataToDisplay <- paste0("Explanation of ",'"',as.character(infoTab[infoTab$V1 == selVar,1]),'"',": ",as.character(infoTab[infoTab$V1 == selVar,2]))
  } else {dataToDisplay <- NULL}
  
  res <- list(sumPlot = sumPlot, sumTab = sumTab, dataToDisplay = dataToDisplay)
  return(res)
  
}






 