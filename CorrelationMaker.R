multiCorr <- function(proteinName, origin, cutOff, clinicalTable, expTable, proteins){
  if (origin == "All") {
    origNames <- clinicalTable[,colnames(clinicalTable) == "Sample_ID"]
  } else if (origin != "All") {
    origNames <- clinicalTable[clinicalTable$SampleType == origin,colnames(clinicalTable) == "Sample_ID"]
  }
  
  #create a matrix of the raw data with filtered NAs, including only those patients who have expression data
  youNeedThesePatients <- expTable[expTable$UNIPROT_ID %in% proteinName,]
  youNeedThesePatients <- t(youNeedThesePatients[,colnames(youNeedThesePatients) %in% origNames])
  youNeedThesePatients <- na.omit(youNeedThesePatients)
  
  rownames(expTable) <- expTable$UNIPROT_ID 
  expTable <- expTable[,!colnames(expTable) == "UNIPROT_ID"]
  
  if (length(youNeedThesePatients) >= 7) {
    tableForAnalysis <- as.matrix(na.omit(expTable[,colnames(expTable) %in% rownames(youNeedThesePatients)]))
    
    #make the selected gene to a numeric vector and perform correlation test
    selectedGene <- as.vector(unlist(youNeedThesePatients))
    corResult <- apply(tableForAnalysis, 1, function(x) cor.test(x, selectedGene, method = "spearman", exact = F)[c(4,3)])
    
    #store the results as a dataframe and convert the list to a dataframe
    result = as.data.frame(matrix(nrow = nrow(tableForAnalysis), ncol = 3))
    colnames(result) <- c("Protein Name" , "Correlation coefficient", "p-value")
    
    result[, 1] = rownames(tableForAnalysis[rownames(tableForAnalysis) %in% names(corResult),])
    result[, 2] = round(x = sapply(corResult, "[[", 1), digits = 2)
    result[, 3] = format(x = sapply(corResult, "[[", 2), digits = 3, scientific = T)
    
    #perform p adjustment, and add gene names and description
    result$pAdj <- format(x = p.adjust(p = result$`p-value`, method = "fdr"), digits = 3, scientific = T)
    result <- merge(x = result, y = data.frame("UNIPROT" = proteins, "Genes" = names(proteins)),
                         by.x = "Protein Name", by.y = "UNIPROT")
    result <- result[result$`Correlation coefficient` >= min(cutOff) & result$`Correlation coefficient` <=max(cutOff),]
    result <- result[order(result$`Correlation coefficient` ,decreasing = T),]
    
  } else if (length(youNeedThesePatients) <= 7) {
    stop("Available patient number with this protein is too low to perform statistical analysis")
    
  }
  
  return(result)
}
singleCorr <- function(proteinA, proteinB, origin, 
                                clinicalTable, expTable, proteins){
  if (origin == "All") {
    origNames <- clinicalTable[,colnames(clinicalTable) == "Sample_ID"]
  } else if (origin != "All") {
    origNames <- clinicalTable[clinicalTable$SampleType == origin,colnames(clinicalTable) == "Sample_ID"]
  }
  
  rownames(expTable) <- expTable$UNIPROT_ID 
  expTable <- expTable[,!colnames(expTable) == "UNIPROT_ID"]
  
  proteinAsamples <- as.data.frame(t(expTable[rownames(expTable) %in% proteinA,]))
  proteinBsamples <- as.data.frame(t(expTable[rownames(expTable) %in% proteinB,]))
  proteinAsamples$ID <- row.names(proteinAsamples)
  proteinBsamples$ID <- row.names(proteinBsamples)
  
  proteinAsamples <- proteinAsamples[proteinAsamples$ID %in% origNames,]
  proteinBsamples <- proteinBsamples[proteinBsamples$ID %in% origNames,]
  
  proteinAsamples <- proteinAsamples[!is.na(proteinAsamples[,1]),]
  proteinBsamples <- proteinBsamples[!is.na(proteinBsamples[,1]),]
  
  #find common patients
  commons <- intersect(rownames(proteinAsamples), rownames(proteinBsamples))
  
  #analysis do not run if sample N is too low
  if (length(commons) <= 7) {
    stop("Available patient number with this protein is too low to perform statistical analysis")
  } else if (length(commons) >= 7) {
    #make expression vectors from patients with common values 
    proteinAExp <- proteinAsamples[proteinAsamples$ID %in% commons, 1]
    proteinBExp <- proteinBsamples[proteinBsamples$ID %in% commons, 1]
    corrTestPearson <- cor.test(x = proteinAExp, y = proteinBExp, method = "pearson", exact = F)
    corrTestSpearman <- cor.test(x = proteinAExp, y = proteinBExp, method = "spearman", exact = F)
    summRes <- data.frame("protANum" = length(proteinAExp),
                          "protBNum" = length(proteinBExp),
                          "R" = c(corrTestPearson$estimate, corrTestSpearman$estimate),
                          "pValue" = c(corrTestPearson$p.value, corrTestSpearman$p.value),
                          "Correlation type" = c("Pearson", "Spearman"))
    
    colnames(summRes) <- c(paste0(proteinA, " N"), paste0(proteinB, " N"), "R", "p", "Analysis type")
    
    
    #create df for visualization
    dat <- data.frame("proteinAExp" = proteinAExp, "proteinBExp" = proteinBExp)
    
    scatPlot <- ggplot(dat, aes(x = proteinAExp, y = proteinBExp)) +
      geom_smooth(method = "lm",linewidth = 2) +
      geom_point(size = 4) +
      xlab(paste0("Log2 Intensity of ", proteinA)) + 
      ylab(paste0("Log2 Intensity of \n", proteinB)) +
      theme(legend.position = "none", 
            panel.background = element_rect(fill = NA, colour = "grey50"),  
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 18, colour = "black"),
            axis.text.y = element_text(size = 18, colour = "black"),
            axis.title.y = element_text(size = 18, face = "bold",  vjust = 2),
            axis.title.x = element_text(size = 18, face = "bold",  vjust = 2))
    res <- list(scatPlot = scatPlot, SumStat = summRes)
    }
  return(res)

}

