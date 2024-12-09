# library(ggplot2)
# library(plotly)
# ##Cero data
# clinicalTableCero <- readRDS(file = "ceroData/ceroData.rds")[[1]]
# expTableCero <- readRDS(file = "ceroData/ceroData.rds")[[2]]
# proteinNamesCero <- readRDS(file = "ceroData/ceroData.rds")[[3]]
# #Primero data
# clinicalTablePrimero <- readRDS(file = "primeroData/primeroData.rds")[[1]]
# expTablePrimero <- readRDS(file = "primeroData/primeroData.rds")[[2]]
# proteinNamesPrimero <- readRDS(file = "primeroData/primeroData.rds")[[3]]
# ##Segundo data
# clinicalTableSegundo <- readRDS(file = "segundoData/segundoData.rds")[[1]]
# expTableSegundo <- readRDS(file = "segundoData/segundoData.rds")[[2]]
# proteinNamesSegundo <- readRDS(file = "segundoData/segundoData.rds")[[3]]
# #Tercero data
# clinTableTercero <- readRDS(file = "terceroData/terceroData.rds")[[1]]
# expTableTercero <- readRDS(file = "terceroData/terceroData.rds")[[2]]
# proteinsNamesTercero <- readRDS(file = "terceroData/terceroData.rds")[[3]]
# #Cuarto Data
# clinTableCuarto <- readRDS(file = "cuartoData/cuartoData.rds")[[1]]
# expTableCuarto <- readRDS(file = "cuartoData/cuartoData.rds")[[2]]
# proteinNamesCuarto <- readRDS(file = "cuartoData/cuartoData.rds")[[3]]
# #
# #
# proteinName <-  proteinNamesCuarto[8594:8596]
# xAxis <- "SampleType"
# jitterVariable <- "Age"
# ## histParams <- list("Tumor_T" = c(0,99),
# ##                    "Necrosis_T" = c(0,37.5),
# ##                    "AdjTissue_T" = c(0,95))
# histVars <- "Tumor_content"
# histParams <- c(0,99)
# expTable <- expTableCuarto
# clinTable <- clinTableCuarto
# proteins <- proteinNamesCuarto


groupComparer <- function(proteinName, xAxis,  jitterVariable,histVars, histParams, expTable, clinTable,proteins){
  if (length(proteinName) == 1) {
    selectedVariables <-  which(colnames(clinTable) %in% c(xAxis, jitterVariable))
    selectedHistVariables <- which(colnames(clinTable) %in% histVars)
    reducedClin <- clinTable[,c(which(colnames(clinTable) == "Sample_ID"),selectedVariables, selectedHistVariables)]
    #reducedClin <- na.omit(reducedClin)
    
    #first hist param filter
    reducedClin2 <- reducedClin[reducedClin[,which(colnames(reducedClin) == histVars)] >= as.numeric(min(histParams)),]
    reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == histVars)] <= as.numeric(max(histParams)),]
    # #second hist param filter
    # reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == names(histParams)[2])] >= as.numeric(min(histParams[[2]])),]
    # reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == names(histParams)[2])] <= as.numeric(max(histParams[[2]])),]
    # #third hist param filter
    # reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == names(histParams)[3])] >= as.numeric(min(histParams[[3]])),]
    # reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == names(histParams)[3])] <= as.numeric(max(histParams[[3]])),]
    
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
    
    # if (identical(expTable,expTableCuarto)) {
    #   if (xAxis == "SampleType") {
    #     tableForPlot$xAxisVariable <- factor(x = tableForPlot$xAxisVariable,levels = c("NT","TM", "PT", "CM", "LN", "DM", "LR"))
    #     
    #   }
    #   if (jitterVariable == "SampleType") {
    #     tableForPlot$jitterVariable <- factor(x = tableForPlot$jitterVariable,levels = c("NT","TM", "PT", "CM", "LN", "DM", "LR"))
    #     
    #   }
    #   
    # }
    
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
    
    #boxPlot
    
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
      #"Gene" = protNames[protNames$protein == theProt, 2])
      
      table$proteinIntensity <- as.numeric(table$proteinIntensity)
      tableForPlot <- rbind(tableForPlot, table)
      
      
    }
    
    # if (xAxis %in% c("Patient",	"BRAF",	"NRAS",	"cKIT",	"Origin",	"Gender",	"Primary.loc",	"Clark",	"Exulceration",
    #                  "Histological.subtype",	"Treatments",	"Treatments.grouped")) {
    #   tableForPlot$xAxisVariable <- as.factor(tableForPlot$xAxisVariable)
    # }
    
    
    boxPlot <- ggplot(data = tableForPlot, aes(x = protein, y = proteinIntensity, fill = xAxisVariable)) +
      geom_violin(alpha =1)+
      geom_vline(xintercept=seq(1.5, length(proteinName)-0.5, 1), colour="black")+
      #geom_violin(alpha = 1)+
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
    
    #ggplotly(boxPlot)
    downloadTable <- tableForPlot
    colnames(downloadTable)[2] <- xAxis
    colnames(downloadTable)[3] <- jitterVariable
    
    res <- list(boxPlot = boxPlot, downloadTable = downloadTable)
    
    
    
  }
  
  return(res)
  
  
}

# 
# test <- groupComparer(proteinName = proteinNamesCuarto[8596],
#                         xAxis = "SampleType",
#                         jitterVariable = "SampleType",
#                         histParams = c(0,100),
#                         expTable = expTableCuarto,
#                         clinTable = clinTableCuarto,
#                         proteins = proteinNamesCuarto,
#                         histVars = "Tumor_content")
# test$downloadTable
# test$boxPlot
