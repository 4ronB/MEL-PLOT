allPlotter <- function(selProt, allClins,proteinNamesCero, proteinNamesPrimero,
                       proteinNamesSegundo, proteinsNamesTercero, proteinNamesCuarto,
                       expTableCero, expTablePrimero, expTableSegundo, expTableTercero, expTableCuarto) {
  #Cero
  if (selProt %in% expTableCero$UNIPROT_ID == TRUE) {
    clinTabCero <- allClins[allClins$Study == "Cero",]
    clinTabCero <- clinTabCero[order(clinTabCero$Sample_ID),]
    selExps <- as.data.frame(t(expTableCero[expTableCero$UNIPROT_ID == selProt,-1]))
    selExps$IDs <- rownames(selExps)
    colnames(selExps) <- c("ProtExp", "IDs")
    selExps <- selExps[selExps$IDs %in% clinTabCero$Sample_ID,]
    selExps <- selExps[order(selExps$IDs),]
    identical(selExps$IDs, clinTabCero$Sample_ID)
    clinTabCero$ProtExp <- selExps$ProtExp
    
  } else if (selProt %in% expTableCero$UNIPROT_ID == FALSE) {
    clinTabCero <- NA
  }
  #Segundo
  if (selProt %in% expTableSegundo$UNIPROT_ID == TRUE) {
    clinicalTableSegundo <- allClins[allClins$Study == "Segundo",]
    selExps <- as.data.frame(t(expTableSegundo[expTableSegundo$UNIPROT_ID == selProt,-1]))
    selExps$IDs <- rownames(selExps)
    colnames(selExps) <- c("ProtExp", "IDs")
    selExps <- selExps[selExps$IDs %in% clinicalTableSegundo$Sample_ID,]
    clinicalTableSegundo <- clinicalTableSegundo[clinicalTableSegundo$Sample_ID %in% selExps$IDs,]
    clinicalTableSegundo <- clinicalTableSegundo[order(clinicalTableSegundo$Sample_ID),]
    selExps <- selExps[order(selExps$IDs),]
    identical(selExps$IDs, clinicalTableSegundo$Sample_ID)
    clinicalTableSegundo$ProtExp <- selExps$ProtExp
    
  } else if (selProt %in% expTableSegundo$UNIPROT_ID == FALSE) {
    clinicalTableSegundo <- NA
  }
  #Cuarto
  if (selProt %in% expTableCuarto$UNIPROT_ID == TRUE) {
    clinTabCuarto <- allClins[allClins$Study == "Cuarto",]
    clinTabCuarto <- clinTabCuarto[order(clinTabCuarto$Sample_ID),]
    selExps <- as.data.frame(t(expTableCuarto[expTableCuarto$UNIPROT_ID == selProt,-1]))
    selExps$IDs <- rownames(selExps)
    colnames(selExps) <- c("ProtExp", "IDs")
    selExps <- selExps[order(selExps$IDs),]
    identical(selExps$IDs, clinTabCuarto$Sample_ID)
    clinTabCuarto$ProtExp <- selExps$ProtExp
    
  } else if (selProt %in% expTableCuarto$UNIPROT_ID == FALSE) {
    clinTabCuarto <- NA
  }
  plotTab <- rbind(clinTabCero, 
                   clinicalTableSegundo, 
                   clinTabCuarto)
  
  plotTab$SampleType <- factor(x = plotTab$SampleType, 
                            levels = c("LN", "PT", "CM", "DM","LR","NT","TM"))
  
  plotTab$Study <- factor(x = plotTab$Study, levels = c("Cero", "Segundo", "Cuarto"))
  plotTab <- plotTab[,colnames(plotTab) %in% c("Sample_ID","SampleType","Study","ProtExp")]
  plotTab <- plotTab[!is.na(plotTab$SampleType),]
  plotTab <- plotTab[!is.na(plotTab$ProtExp),]

  allSumPlot <- ggplot(data = plotTab, aes(x = SampleType,y = ProtExp, fill = SampleType)) +
    geom_boxplot() +
    geom_jitter(width = 0.25) +
    labs(x = NULL, y = paste0("Log2 intensity of\n",selProt)) +
    facet_wrap(~Study, scales = "free")+
    theme(legend.position = "bottom",
      panel.background = element_rect(fill = NA, colour = "grey50"),
      strip.text.x = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"),
      axis.text.x = element_text(size = 13, colour = "black"),
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title.y = element_text(size = 18,  vjust = 2, face = "bold"),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18))

  tabToDownload <- plotTab[,-1]

  return(list(allSumPlot = allSumPlot, tabToDownload = tabToDownload))
  
}




