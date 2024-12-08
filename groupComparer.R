groupComparer <- function(proteinName, xAxis,  jitterVariable,histVars, histParams, expTable, clinTable,proteins){
  if (length(proteinName) == 1) {
    selectedVariables <-  which(colnames(clinTable) %in% c(xAxis, jitterVariable))
    selectedHistVariables <- which(colnames(clinTable) %in% histVars)
    reducedClin <- clinTable[,c(which(colnames(clinTable) == "Sample_ID"),selectedVariables, selectedHistVariables)]
    reducedClin2 <- reducedClin[reducedClin[,which(colnames(reducedClin) == histVars)] >= as.numeric(min(histParams)),]
    reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == histVars)] <= as.numeric(max(histParams)),]
    reducedClin2 <- na.omit(reducedClin2)
    
    protExp <- expTable[expTable$UNIPROT_ID == proteinName,]
    protExp <- as.data.frame(t(protExp))
    protExp$Sample_ID <- rownames(protExp)
    protExp <- protExp[-1,]
    colnames(protExp) <- c("ProteinExpression", "Sample_ID")
    
    reducedClin3 <- merge(reducedClin2, protExp)
    
    tableForPlot <- data.frame("proteinIntensity" = reducedClin3$ProteinExpression,
                               "xAxisVariable" = reducedClin3[,colnames(reducedClin3) == xAxis],
                               "jitterVariable" = reducedClin3[,colnames(reducedClin3) == jitterVariable])
    tableForPlot$proteinIntensity <- as.numeric(tableForPlot$proteinIntensity)
    
    boxPlot <- ggplot(data = tableForPlot, aes(y = proteinIntensity, x = xAxisVariable)) +
      geom_violin(alpha = 0.6, scale = "width")+
      geom_jitter(aes(colour = jitterVariable),fill = jitterVariable) +
      {if (jitterVariable == xAxis) {
        scale_color_manual(values = as.character(rep.int(x = "black",times = length(unique(tableForPlot$jitterVariable)))))}}+
      labs(colour=jitterVariable, fill = xAxis, x = xAxis, 
           y = paste0("Log2 intensity of \n", names(proteins[proteins %in% proteinName])))+
      theme_minimal() +
      theme(legend.position="bottom",
            axis.text.x = element_text(size = 13, colour = "black"),
            axis.text.y = element_text(size = 13, colour = "black"),
            axis.title.y = element_text(size = 13, face = "bold"),
            axis.title.x = element_text(size = 13, face = "bold"))
    downloadTable <- tableForPlot
    colnames(downloadTable)[2] <- xAxis
    colnames(downloadTable)[3] <- jitterVariable
    res <- list(boxPlot = boxPlot, downloadTable = downloadTable)
    
    
    
  } else if (length(proteinName) > 1) {
    selectedVariables <-  which(colnames(clinTable) %in% c(xAxis, jitterVariable))
    selectedHistVariables <- which(colnames(clinTable) %in% histVars)
    reducedClin <- clinTable[,c(which(colnames(clinTable) == "Sample_ID"),selectedVariables, selectedHistVariables)]
    
    reducedClin2 <- reducedClin[reducedClin[,which(colnames(reducedClin) == histVars)] >= as.numeric(min(histParams)),]
    reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == histVars)] <= as.numeric(max(histParams)),]
    reducedClin2 <- na.omit(reducedClin2)
    
    protExp <- expTable[expTable$UNIPROT_ID %in% proteinName,]
    protExp <- as.data.frame(t(protExp))
    colnames(protExp) <- protExp[1,]
    protExp$Sample_ID <- rownames(protExp)
    protExp <- protExp[-1,]
    
    reducedClin3 <- merge(reducedClin2, protExp)
    
    tableForPlot <- data.frame()
    for (protein in 1:(ncol(protExp)-1)) {
      theProt <- colnames(protExp)[protein]
      table <- data.frame("proteinIntensity" = reducedClin3[,colnames(reducedClin3) == theProt],
                          "xAxisVariable" = reducedClin3[,colnames(reducedClin3) == xAxis],
                          "jitterVariable" = reducedClin3[,colnames(reducedClin3) == jitterVariable],
                          "protein" = theProt)
     table$proteinIntensity <- as.numeric(table$proteinIntensity)
      tableForPlot <- rbind(tableForPlot, table)
      
      
    }
    
    boxPlot <- ggplot(data = tableForPlot, aes(x = protein, y = proteinIntensity, fill = xAxisVariable)) +
      geom_vline(xintercept=seq(1.5, length(proteinName)-0.5, 1), colour="black")+
      geom_violin(alpha = 1, scale = "count")+
      geom_jitter(aes(colour = jitterVariable),fill = jitterVariable) +
      {if (jitterVariable == xAxis) {
        scale_color_manual(values = as.character(rep.int(x = "black",times = length(unique(tableForPlot$jitterVariable)))))}} +
      labs(colour=jitterVariable, fill = xAxis, y= NULL, x = "Protein intensities (log2)")+
      guides(fill=guide_legend(title= xAxis)) +
      coord_flip()+
      theme(legend.position="bottom",
            axis.text.x = element_text(size = 13, colour = "black"),
            axis.text.y = element_text(size = 13, colour = "black", face = "bold"),
            axis.title.y = element_text(size = 13, face = "bold"),
            axis.title.x = element_text(size = 13, face = "bold"))
    downloadTable <- tableForPlot
    colnames(downloadTable)[2] <- xAxis
    colnames(downloadTable)[3] <- jitterVariable
    
    res <- list(boxPlot = boxPlot, downloadTable = downloadTable)
  }
  
  return(res)
  
  
}