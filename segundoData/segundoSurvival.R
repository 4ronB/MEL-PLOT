survPlotterSegundo <- function(variable,histParams, expTable, clinicalTableSegundo) {
  library(survival)
  library(survminer)
  redCclinicalTable <- clinicalTableSegundo[clinicalTableSegundo$SampleType == "LN",]
  reducedClin2 <- redCclinicalTable[redCclinicalTable[,which(colnames(redCclinicalTable) == "Tumor_content")] >= as.numeric(min(histParams)),]
  reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == "Tumor_content")] <= as.numeric(max(histParams)),]
  reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == "necrosis")] >= 0,]
  reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == "necrosis")] <= 99,]
  reducedClin2 <- reducedClin2[!is.na(reducedClin2$Type),]
  
  
  usedTable <- data.frame("ID" = reducedClin2$Sample_ID,
                          "Time" = reducedClin2$OS_DAYS, 
                          "Event" = reducedClin2$dss.events,
                          "Remarks" = reducedClin2$`Remarks original`,
                          stringsAsFactors = T)
  usedTable <- usedTable[order(usedTable$ID),]
  
  usedExp <- as.data.frame(t(expTable[expTable$UNIPROT_ID == variable,])[-1,])
  usedExp$IDs <- rownames(usedExp)
  colnames(usedExp) <- c("protExp", "IDs")
  usedExp <- usedExp[order(usedExp$IDs),]
  usedExp <- usedExp[usedExp$IDs %in% intersect(as.character(unlist(usedTable$ID)), usedExp$IDs),]
  usedTable <- usedTable[usedTable$ID %in% intersect(as.character(unlist(usedTable$ID)), usedExp$IDs),]
  if (identical(as.character(unlist(usedTable$ID)), usedExp$IDs) == F) {
    stop("An error has occured, clincal and expression samples are different")
    
  }  
  usedTable <- cbind(usedTable, as.numeric(usedExp$protExp))
  colnames(usedTable) <- c("ID", "Time", "Event", "Remarks", "Value")
  
  usedTable <- usedTable[!is.na(usedTable$Value),]
  usedTable <- usedTable[!is.na(usedTable$Event),]
  usedTable$Event <- ifelse(usedTable$Event == 1, yes = 2, no = 1)
  usedTable <- usedTable[!is.na(usedTable$Time),]
  
  if (length(usedTable$Value) >= 10) {
    cutoff <- function(datavector, cutpoint, clintable) {
      term <- cut(x = datavector, breaks = c(min(datavector), cutpoint, max(datavector)), labels = F, include.lowest = T)
      cox <- summary(coxph(Surv(Time, Event) ~ term, data = usedTable))
      c(cutpoint, cox$sctest[3])
    }
    bestcutoff <- function(datavector, clintable) {
      breaks <- quantile(datavector, probs = seq(0.25, 0.75, by = 0.01))
      cutoff.table <- t(sapply(breaks, function(z) cutoff(datavector = datavector, cutpoint = z, clintable = usedTable)))
      cutoff.table <- cbind(cutoff.table, p.adjust(cutoff.table[, 2], method = "fdr"))
      colnames(cutoff.table) <- c("cutoff", "pvalue", "adj_pvalue")
      cutoff.table[order(cutoff.table[, 3]), "cutoff"][1]
    }
    cutoff.point <-  as.numeric(bestcutoff(datavector = usedTable$Value, clintable = usedTable))
    exp_category <-  c()
    for (c in 1:length(usedTable$Value)){
      if (usedTable$Value[[c]] >= cutoff.point[1]){
        exp_category[[c]] <- "high"
      } else {
        exp_category[[c]] <- "low"
      }
    }
    exp_category <- factor(exp_category, levels = c("low", "high"))
    usedTable$expCat <- exp_category
    survFit <- survfit(Surv(Time, Event) ~ expCat, data = usedTable)
    cox_result <- coxph(Surv(Time, Event) ~ expCat, data = usedTable)
    
    
    results_table <- data.frame("Protein" = variable,
                                "p-values" = format(as.numeric(summary(cox_result)$sctest['pvalue']), digits = 3, scientific = T),
                                "HR" = as.numeric(round(summary(cox_result)$conf.int[1], digits = 2)),
                                "Cutoff-Value" = cutoff.point,
                                "MedSurvLowExp" = median(usedTable[usedTable$Value < cutoff.point,2]),
                                "MedSurvHighExp" = median(usedTable[usedTable$Value >= cutoff.point,2]))
    
    plot <- ggsurvplot(survFit, data = usedTable, pval = TRUE,
                       xlab = "Time (days)",
                       risk.table = T,
                       legend.title = paste0(variable, "\n expression"),
                       legend.labs = c("Low", "High"))
    
    tabToDownload <- usedTable[,colnames(usedTable) %in% c("Time","Event","Value","expCat")]
    colnames(tabToDownload) <- c("Time","Event","ProteinExp","expCat")
    res <- list(survPlot = plot, resTab = results_table, tabToDownload = tabToDownload)
    
  } else if (length(usedTable$Value) < 10) {
    stop("An error has occured, available patient number is too low to perform anaylsis with this protein, please select another protein")
    }
}