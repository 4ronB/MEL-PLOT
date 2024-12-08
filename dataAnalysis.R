library(readxl)
library(ggplot2)
library(dunn.test)
library(enrichplot)
library(GOplot)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(Biobase)
library(msigdbr)
library(singscore)
library(GSEABase)
library(rcartocolor)
library(cowplot)
library(rstatix)
allClinTab <- read_xlsx(path = "clinData.xlsx", sheet = 1)
allClinTab <- as.data.frame(allClinTab)

###Cero data
clinicalTableCero <- as.data.frame(readRDS(file = "CeroData.rds")[[1]])
expTableCero <- readRDS(file = "CeroData.rds")[[2]]
expTableCero <- expTableCero[,c(1,which(colnames(expTableCero) %in% clinicalTableCero$Samples.ID))]
clinicalTableCero <- clinicalTableCero[clinicalTableCero$Samples.ID %in% colnames(expTableCero),]
proteinNamesCero <- read_xlsx("protNames.xlsx", sheet = 3, col_names = F)
proteinNamesCero <- setNames(object = proteinNamesCero$...1, nm = proteinNamesCero$...2)
#Primero data
expTablePrimero <- readRDS(file = "PrimeroData.rds")[[1]]
clinTablePrimero <- readRDS(file = "PrimeroData.rds")[[2]]
proteinNamesPrimero <- read_xlsx("protNames.xlsx", sheet = 2, col_names = F)
proteinNamesPrimero <- setNames(object = proteinNamesPrimero$...1, nm = proteinNamesPrimero$...2)
##Segundo data
clinicalTableSegundo <- readRDS(file = "SegundoData.rds")[[1]]
expTableSegundo <- readRDS(file = "SegundoData.rds")[[2]]
proteinNamesSegundo <- read_xlsx("protNames.xlsx", sheet = 4, col_names = F)
proteinNamesSegundo <- setNames(object = proteinNamesSegundo$...1,nm = proteinNamesSegundo$...2)
#Tercero data
clinTableTercero <- readRDS(file = "TerceroData.rds")[[1]]
expTableTercero <- readRDS(file = "../terceroData/TerceroData.rds")[[2]]
proteinsTercero <- read_xlsx("../protNames.xlsx", sheet = 5, col_names = F)
proteinsTercero <- setNames(object = proteinsTercero$...1,nm = proteinsTercero$...2)
#Cuarto Data
clinTableCuarto <- readRDS(file = "CuartoData.rds")[[1]]
expTableCuarto <- readRDS(file = "CuartoData.rds")[[2]]
proteinNamesCuarto <- read_xlsx("protNames.xlsx", sheet = 6, col_names = F)
proteinNamesCuarto <- setNames(object = proteinNamesCuarto$...1,nm = proteinNamesCuarto$...2)

allNames <- read_xlsx("protNames.xlsx", sheet = 7, col_names = F)
allNames <- setNames(object = allNames$...1,nm = allNames$...2)
##########
##########diff exp analysis between  tissue types
##########
#Cero
table(allClinTab[allClinTab$Study == "Cero", 2])

bigTab <- data.frame(matrix(nrow = nrow(expTableCero),ncol = 28))

for (protein in 1:nrow(expTableCero)) {
  selProt <- expTableCero[protein,1]

  primNames <- clinicalTableCero[clinicalTableCero$SampleType == "Prim_Cut", 2]
  metLNnames <- clinicalTableCero[clinicalTableCero$SampleType == "Met_LN", 2]
  metCutnames <- clinicalTableCero[clinicalTableCero$SampleType == "Met_Cut", 2]

  primExp <- expTableCero[expTableCero$Samples.ID == selProt,colnames(expTableCero) %in% primNames]
  primExp <- primExp[!is.na(primExp)]
  primExp <- setNames(object = primExp,nm = rep(x = "Prim_Cut",times = length(primExp)))

  metLnExp <- expTableCero[expTableCero$Samples.ID == selProt,colnames(expTableCero) %in% metLNnames]
  metLnExp <- metLnExp[!is.na(metLnExp)]
  metLnExp <- setNames(object = metLnExp,nm = rep(x = "Met_LN",times = length(metLnExp)))

  metCutExp <- expTableCero[expTableCero$Samples.ID == selProt,colnames(expTableCero) %in% metCutnames]
  metCutExp <- metCutExp[!is.na(metCutExp)]
  metCutExp <- setNames(object = metCutExp,nm = rep(x = "Met_Cut",times = length(metCutExp)))

  if (length(primExp) > 5 & length(metLnExp) > 5 & length(metCutExp) > 5) {

    effDat <- data.frame("Exp" = c(primExp, metLnExp, metCutExp),
                         "Names" = c(names(primExp),names(metLnExp),names(metCutExp)))

    kw <- kruskal.test(Exp~Names,data = effDat)
    dunnTest <- dunn.test(x = effDat$Exp,g = effDat$Names,label = T, kw = T)
    effectSize <- kruskal_effsize(data = effDat, formula = Exp~Names)
    sumTab <- data.frame("ProteinName" = selProt,
                         "kwP" = kw$p.value,
                         "Met_LN-Met_Cut_pAdj" = dunnTest$P.adjusted[1],
                         "Prim-Met_Cut_pAdj" = dunnTest$P.adjusted[2],
                         "Prim-Met_LN_pAdj" = dunnTest$P.adjusted[3],
                         "effSize" = effectSize$effsize,
                         "effSizeMagnitude" = as.character(effectSize$magnitude),

                         "Prim-Met_LN_FC" = mean(metLnExp) - mean(primExp),
                         "Prim-Met_Cut_FC" = mean(metCutExp) - mean(primExp),
                         "Met_LN-Met_Cut_FC" = mean(metCutExp) - mean(metLnExp),

                         "Mean_Prim" = mean(as.numeric(primExp)),
                         "Mean_Met_LN" = mean(as.numeric(metLnExp)),
                         "Mean_Met_Cut" = mean(as.numeric(metCutExp)),

                         "Min_Prim" = fivenum(as.numeric(primExp))[1],
                         "Q1_Prim" = fivenum(as.numeric(primExp))[2],
                         "Med_Prim" = fivenum(as.numeric(primExp))[3],
                         "Q3_Prim" = fivenum(as.numeric(primExp))[4],
                         "Max_Prim" = fivenum(as.numeric(primExp))[5],

                         "Min_Met_LN" = fivenum(as.numeric(metLnExp))[1],
                         "Q1_Met_LN" = fivenum(as.numeric(metLnExp))[2],
                         "Med_Met_LN" = fivenum(as.numeric(metLnExp))[3],
                         "Q3_Met_LN" = fivenum(as.numeric(metLnExp))[4],
                         "Max_Met_LN" = fivenum(as.numeric(metLnExp))[5],

                         "Min_Met_Cut" = fivenum(as.numeric(metCutExp))[1],
                         "Q1_Met_Cut" = fivenum(as.numeric(metCutExp))[2],
                         "Med_Met_Cut" = fivenum(as.numeric(metCutExp))[3],
                         "Q3_Met_Cut" = fivenum(as.numeric(metCutExp))[4],
                         "Max_Met_Cut" = fivenum(as.numeric(metCutExp))[5]
    )

  } else {

  sumTab <- data.frame("ProteinName" = selProt,
                       "kwP" = NA,
                       "Met_LN-Met_Cut_pAdj" = NA,
                       "Prim-Met_Cut_pAdj" = NA,
                       "Prim-Met_LN_pAdj" = NA,
                       "effSize" = NA,
                       "effSizeMagnitude" = NA,

                       "Prim-Met_LN_FC" = mean(metLnExp) - mean(primExp),
                       "Prim-Met_Cut_FC" = mean(metCutExp) - mean(primExp),
                       "Met_LN-Met_Cut_FC" = mean(metCutExp) - mean(metLnExp),

                       "Mean_Prim" = mean(as.numeric(primExp)),
                       "Mean_Met_LN" = mean(as.numeric(metLnExp)),
                       "Mean_Met_Cut" = mean(as.numeric(metCutExp)),

                       "Min_Prim" = fivenum(as.numeric(primExp))[1],
                       "Q1_Prim" = fivenum(as.numeric(primExp))[2],
                       "Med_Prim" = fivenum(as.numeric(primExp))[3],
                       "Q3_Prim" = fivenum(as.numeric(primExp))[4],
                       "Max_Prim" = fivenum(as.numeric(primExp))[5],

                       "Min_Met_LN" = fivenum(as.numeric(metLnExp))[1],
                       "Q1_Met_LN" = fivenum(as.numeric(metLnExp))[2],
                       "Med_Met_LN" = fivenum(as.numeric(metLnExp))[3],
                       "Q3_Met_LN" = fivenum(as.numeric(metLnExp))[4],
                       "Max_Met_LN" = fivenum(as.numeric(metLnExp))[5],

                       "Min_Met_Cut" = fivenum(as.numeric(metCutExp))[1],
                       "Q1_Met_Cut" = fivenum(as.numeric(metCutExp))[2],
                       "Med_Met_Cut" = fivenum(as.numeric(metCutExp))[3],
                       "Q3_Met_Cut" = fivenum(as.numeric(metCutExp))[4],
                       "Max_Met_Cut" = fivenum(as.numeric(metCutExp))[5]
  )}

  bigTab[protein,] <- sumTab
}

colnames(bigTab) <- colnames(sumTab)
write.table(x = bigTab, file = "CeroTypeDiffExpSum_effectSize.csv", sep = ";", row.names = F)

### cuarto diff exp
table(allClinTab[allClinTab$Study == "Cuarto", 2])

bigTab <- data.frame(matrix(nrow = nrow(expTableCuarto),ncol = 28))

for (protein in 1:nrow(expTableCuarto)) {
  selProt <- expTableCuarto[protein,1]
  types <- names(table(allClinTab[allClinTab$Study == "Cuarto", 2]))
  types <- c("Non_Tum","Prim_Cut", "Met_LN")

  normNames <- allClinTab[allClinTab$Study == "Cuarto" & allClinTab$SampleType == "NT", 1]
  primNames <- allClinTab[allClinTab$Study == "Cuarto" & allClinTab$SampleType == "PT", 1]
  metLNnames <- allClinTab[allClinTab$Study == "Cuarto" & allClinTab$SampleType == "LN", 1]

  normExp <- expTableCuarto[expTableCuarto$PG.ProteinAccessions == selProt,colnames(expTableCuarto) %in% normNames]
  normExp <- normExp[!is.na(normExp)]
  normExp <- setNames(object = normExp,nm = rep(x = "NT",times = length(normExp)))

  primExp <- expTableCuarto[expTableCuarto$PG.ProteinAccessions == selProt,colnames(expTableCuarto) %in% primNames]
  primExp <- primExp[!is.na(primExp)]
  primExp <- setNames(object = primExp,nm = rep(x = "PT",times = length(primExp)))

  metLnExp <- expTableCuarto[expTableCuarto$PG.ProteinAccessions == selProt,colnames(expTableCuarto) %in% metLNnames]
  metLnExp <- metLnExp[!is.na(metLnExp)]
  metLnExp <- setNames(object = metLnExp,nm = rep(x = "LN",times = length(metLnExp)))

  if (length(primExp) > 5 & length(metLnExp) > 5 & length(normExp) > 5) {
    effDat <- data.frame("Exp" = c(normExp, primExp, metLnExp),
                         "Names" = c(names(normExp),names(primExp),names(metLnExp)))

    kw <- kruskal.test(Exp~Names,data = effDat)
    dunnTest <- dunn.test(x = effDat$Exp,g = effDat$Names,label = T, kw = T)
    effectSize <- kruskal_effsize(data = effDat, formula = Exp~Names)
    sumTab <- data.frame("ProteinName" = selProt,
                         "kwP" = kw$p.value,
                         "Non_Tum-Met_LN_pAdj" = dunnTest$P.adjusted[1],
                         "Prim-Met_LN_pAdj" = dunnTest$P.adjusted[2],
                         "Non_Tum-Prim_pAdj" = dunnTest$P.adjusted[3],
                         "effSize" = effectSize$effsize,
                         "effSizeMagnitude" = as.character(effectSize$magnitude),

                         "Non_Tum-Prim_FC" = mean(primExp) - mean(normExp),
                         "Prim-Met_LN_FC" = mean(metLnExp) - mean(primExp),
                         "Non_Tum-Met_LN_FC" = mean(metLnExp) - mean(normExp),

                         "Mean_Non_Tum" = mean(as.numeric(normExp)),
                         "Mean_Prim" = mean(as.numeric(primExp)),
                         "Mean_Met_LN" = mean(as.numeric(metLnExp)),

                         "Min_Non_Tum" = fivenum(as.numeric(normExp))[1],
                         "Q1_Non_Tum" = fivenum(as.numeric(normExp))[2],
                         "Med_Non_Tum" = fivenum(as.numeric(normExp))[3],
                         "Q3_Non_Tum" = fivenum(as.numeric(normExp))[4],
                         "Max_Non_Tum" = fivenum(as.numeric(normExp))[5],

                         "Min_Prim" = fivenum(as.numeric(primExp))[1],
                         "Q1_Prim" = fivenum(as.numeric(primExp))[2],
                         "Med_Prim" = fivenum(as.numeric(primExp))[3],
                         "Q3_Prim" = fivenum(as.numeric(primExp))[4],
                         "Max_Prim" = fivenum(as.numeric(primExp))[5],

                         "Min_Met_LN" = fivenum(as.numeric(metLnExp))[1],
                         "Q1_Met_LN" = fivenum(as.numeric(metLnExp))[2],
                         "Med_Met_LN" = fivenum(as.numeric(metLnExp))[3],
                         "Q3_Met_LN" = fivenum(as.numeric(metLnExp))[4],
                         "Max_Met_LN" = fivenum(as.numeric(metLnExp))[5]
                         )

  } else {

    sumTab <- data.frame("ProteinName" = selProt,
                         "kwP" = NA,
                         "Non_Tum-Met_LN_pAdj" = NA,
                         "Prim-Met_LN_pAdj" = NA,
                         "Non_Tum-Prim_pAdj" = NA,
                         "effSize" = NA,
                         "effSizeMagnitude" = NA,

                         "Non_Tum-Prim_FC" = mean(primExp) - mean(normExp),
                         "Prim-Met_LN_FC" = mean(metLnExp) - mean(primExp),
                         "Non_Tum-Met_LN_FC" = mean(metLnExp) - mean(normExp),

                         "Mean_Non_Tum" = mean(as.numeric(normExp)),
                         "Mean_Prim" = mean(as.numeric(primExp)),
                         "Mean_Met_LN" = mean(as.numeric(metLnExp)),

                         "Min_Non_Tum" = fivenum(as.numeric(normExp))[1],
                         "Q1_Non_Tum" = fivenum(as.numeric(normExp))[2],
                         "Med_Non_Tum" = fivenum(as.numeric(normExp))[3],
                         "Q3_Non_Tum" = fivenum(as.numeric(normExp))[4],
                         "Max_Non_Tum" = fivenum(as.numeric(normExp))[5],

                         "Min_Prim" = fivenum(as.numeric(primExp))[1],
                         "Q1_Prim" = fivenum(as.numeric(primExp))[2],
                         "Med_Prim" = fivenum(as.numeric(primExp))[3],
                         "Q3_Prim" = fivenum(as.numeric(primExp))[4],
                         "Max_Prim" = fivenum(as.numeric(primExp))[5],

                         "Min_Met_LN" = fivenum(as.numeric(metLnExp))[1],
                         "Q1_Met_LN" = fivenum(as.numeric(metLnExp))[2],
                         "Med_Met_LN" = fivenum(as.numeric(metLnExp))[3],
                         "Q3_Met_LN" = fivenum(as.numeric(metLnExp))[4],
                         "Max_Met_LN" = fivenum(as.numeric(metLnExp))[5]
    )}

  bigTab[protein,] <- sumTab

}

colnames(bigTab) <- colnames(sumTab)
write.table(x = bigTab, file = "cuartoTypeDiffExpSum_effectSize.csv", sep = ";", row.names = F)

#########################
######################### survival analysis of the common proteins using Segundo Met_LN data
#########################

survPlotterSegundo <- function(variable ,histParams, expTableSegundo, clinicalTableSegundo) {
  library(survival)
  library(survminer)
  redCclinicalTable <- clinicalTableSegundo[clinicalTableSegundo$Type == "Met_LN",]
  protSYMBOL <- filteredProteins_forSurvAn[filteredProteins_forSurvAn$ProteinID == variable,2]
  #Tumor_content filter
  reducedClin2 <- redCclinicalTable[redCclinicalTable[,which(colnames(redCclinicalTable) == "Tumor_content")] >= 0,]
  reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == "Tumor_content")] <= 99,]
  #necrosis filter
  reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == "necrosis")] >= 0,]
  reducedClin2 <- reducedClin2[reducedClin2[,which(colnames(reducedClin2) == "necrosis")] <= 99,]
  reducedClin2 <- reducedClin2[!is.na(reducedClin2$Type),]

  usedTable <- data.frame("ID" = reducedClin2$sample,
                          "Time" = reducedClin2$`OS.DAYS     days from primary diagnosis to death or censoring`,
                          "Event" = reducedClin2$dss.events,
                          "Remarks" = reducedClin2$`Remarks original`,
                          stringsAsFactors = T)
  usedTable <- usedTable[order(usedTable$ID),]

  usedExp <- as.data.frame(t(expTableSegundo[expTableSegundo$Accession == variable,])[-1,])
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

  if (length(usedTable$Value) >= 0) {
    
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
    results_table <- data.frame("ProteinID" = variable,
                                "SYMBOL" = protSYMBOL,
                                "p-values" = format(as.numeric(summary(cox_result)$sctest['pvalue']), digits = 3, scientific = T),
                                "HR" = as.numeric(round(summary(cox_result)$conf.int[1], digits = 2)),
                                "Cutoff-Value" = cutoff.point,
                                "MedSurvLowExp" = median(usedTable[usedTable$Value < cutoff.point,2]),
                                "MedSurvHighExp" = median(usedTable[usedTable$Value >= cutoff.point,2]))

    plot <- ggsurvplot(survFit, data = usedTable, pval = FALSE,
                       xlab = "Time (days)",
                       risk.table = T,
                       legend.title = paste0(protSYMBOL, "\n expression"),
                       legend.labs = c("Low", "High"))

    #plot
    res <- list(survPlot = plot, resTab = results_table)

  } else if (length(usedTable$Value) < 10) {
    print("available patient number is too low")
  }


}

bigTab <- data.frame()
filteredProteins_forSurvAn <- filteredSurv 

for (protein in 1:nrow(filteredProteins_forSurvAn)) {
  selProtein <- filteredProteins_forSurvAn$ProteinID[44]

  result <- survPlotterSegundo(variable = selProtein,
                               histParams = list("Tumor_content" = c(0,99),
                                                 "necrosis" = c(0,99)),
                               expTableSegundo = expTableSegundo,
                               clinicalTableSegundo = clinicalTableSegundo)
  result$resTab
  bigTab <- rbind(bigTab, result$resTab)
}

bigTab$pVal <- as.numeric(bigTab$p.values)
write.table(x = bigTab, file = "segundoSurvSum.csv", sep = ";", row.names = F)


# working with the tables

CeroDiffExpSum <- read.table(file = "CeroTypeDiffExpSum_effectSize.csv", header = T, sep = ";")
CeroDiffExpSum$pAdj <- p.adjust(CeroDiffExpSum$kwP, method = "fdr")
CeroDiffExpSumRed <- CeroDiffExpSum[CeroDiffExpSum$pAdj <= 0.05,]
CeroDiffExpSumRed <- CeroDiffExpSumRed[!is.na(CeroDiffExpSumRed$pAdj),]

cuartoDiffExpSum <- read.table(file = "cuartoTypeDiffExpSum_effectSize.csv", header = T, sep = ";")
cuartoDiffExpSum$pAdj <- p.adjust(p = cuartoDiffExpSum$kwP, method = "fdr")
cuartoDiffExpSumRed <- cuartoDiffExpSum[cuartoDiffExpSum$pAdj <= 0.05,]
cuartoDiffExpSumRed <- cuartoDiffExpSumRed[!is.na(cuartoDiffExpSumRed$pAdj),]

expTableSegundo <- readRDS(file = "SegundoData.rds")[[2]]
segundoSruv <- read.table(file = "segundoSurvSum.csv", header = T, sep = ";")
segundoSruv <- segundoSruv[!duplicated(segundoSruv$ProteinID),]

filteredProteins <- intersect(cuartoDiffExpSumRed$ProteinName, intersect(expTableSegundo$Accession, CeroDiffExpSumRed$ProteinName))
cuartoRedDiffExp <- cuartoDiffExpSum[cuartoDiffExpSum$ProteinName %in% filteredProteins,]
CeroRedDiffExp <- CeroDiffExpSum[CeroDiffExpSum$ProteinName %in% filteredProteins,]
#calculate effect size ranges
range(CeroRedDiffExp$effSize)
range(cuartoRedDiffExp$effSize)

filteredSurv <- segundoSruv[segundoSruv$ProteinID %in% filteredProteins,]
filteredSurv$p.adj <- p.adjust(filteredSurv$p.values,method = "fdr")
filteredSurvSignif <- filteredSurv[filteredSurv$p.adj < 0.05,]
cuartoRed <- cuartoDiffExpSum[cuartoDiffExpSum$ProteinName %in% filteredSurvSignif$ProteinID,]
CeroRed <- CeroDiffExpSum[CeroDiffExpSum$ProteinName %in% filteredSurvSignif$ProteinID,]

#proteins with HR < 0.5
smallHRs <- filteredSurvSignif[filteredSurvSignif$HR < 0.5,]
#proteins with HR > 2
bigHRs <- filteredSurvSignif[filteredSurvSignif$HR > 2,]

#Gene Ontology
ego <- enrichGO(gene = filteredSurv$ProteinID,
                keyType = "UNIPROT",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.2,
                minGSSize = 10,
                maxGSSize = 500,
                readable = TRUE)

dotplot(ego)

ego2 <- pairwise_termsim(ego2)
treePlot <- treeplot(ego)

egoSurv <- enrichGO(gene = smallHRs$ProteinID,
                    keyType = "UNIPROT",
                    OrgDb = org.Hs.eg.db,
                    ont = "CC",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2,
                    minGSSize = 10,
                    maxGSSize = 500,
                    readable = TRUE)

dotPlot_Surv <- dotplot(egoSurv)
dotPlot_Surv

ego2Surv <- pairwise_termsim(egoSurv)
treePlotSurv <- treeplot(ego2Surv)
treePlotSurv

######generate heatmap for cuarto and Cero
#################Cero
patientAnn <- allClinTab[allClinTab$Study == "Cero",]

geneAnnTab <- filteredSurvSignif[,c(1,2)]
CeroRedExp <- expTableCero[expTableCero$Samples.ID %in% filteredSurvSignif$ProteinID,]
CeroRedExp <- merge(x = geneAnnTab,y = CeroRedExp, by.x = "ProteinID", by.y = "Samples.ID")
rownames(CeroRedExp) <- CeroRedExp$SYMBOL
#ordering
CeroRedExp <- CeroRedExp[,colnames(CeroRedExp) %in% patientAnn$Sample_ID]
CeroRedExp <- CeroRedExp[,order(colnames(CeroRedExp))]
patientAnn <- patientAnn[order(patientAnn$Sample_ID),]

identical(colnames(CeroRedExp), patientAnn$Sample_ID)

CeroRedExp <- as.matrix(CeroRedExp)
CeroRedExp <- t(scale(t(CeroRedExp)))
col_fun = colorRamp2(c(min(CeroRedExp,na.rm = T), median(CeroRedExp,na.rm = T), max(CeroRedExp,na.rm = T)), c("royalblue", "white", "red2"))

ta <- HeatmapAnnotation(Type = patientAnn$SampleType,
                        BRAF = patientAnn$BRAFstatus,
                        Gender = patientAnn$Sex,
                        Age = patientAnn$Age,
                        TC = patientAnn$Tumor_content,
                        
                        col = list("Type" = setNames(object = c(carto_pal(12,"Vivid")[1:length(unique(patientAnn$SampleType))]),nm = as.character(unique(patientAnn$SampleType)))
                                   )
                        )

hm1<-Heatmap(CeroRedExp,
             name="log2(Intensity)",
             top_annotation = ta,
             show_column_dend = TRUE,
             show_row_dend = T,
             row_names_gp = gpar(fontsize(20), fontface = "bold"),
             show_column_names = FALSE,
             rect_gp = gpar(col = "grey30", lwd = 0.5),
             width = unit(8,"in"),
             height = unit(8,"in"),
             row_title = NULL
)

#################Cuarto
patientAnn <- allClinTab[allClinTab$Study == "Cuarto",]
patientAnn <- patientAnn[patientAnn$SampleType %in% c("NT","PT","LN"),]
geneAnnTab <- filteredSurvSignif[,c(1,2)]
cuartoRedExp <- expTableCuarto[expTableCuarto$PG.ProteinAccessions %in% filteredSurvSignif$ProteinID,]
cuartoRedExp <- merge(x = geneAnnTab,y = cuartoRedExp, by.x = "ProteinID", by.y = "PG.ProteinAccessions")
rownames(cuartoRedExp) <- cuartoRedExp$SYMBOL
#ordering
cuartoRedExp <- cuartoRedExp[,colnames(cuartoRedExp) %in% patientAnn$Sample_ID]
cuartoRedExp <- cuartoRedExp[,order(colnames(cuartoRedExp))]
patientAnn <- patientAnn[order(patientAnn$Sample_ID),]
identical(colnames(cuartoRedExp), patientAnn$Sample_ID)

cuartoRedExp <- as.matrix(cuartoRedExp)
cuartoRedExp <- t(scale(t(cuartoRedExp)))

col_fun = colorRamp2(c(min(cuartoRedExp,na.rm = T), median(cuartoRedExp,na.rm = T), max(cuartoRedExp,na.rm = T)), c("royalblue", "white", "red2"))
ta <- HeatmapAnnotation(Type = patientAnn$SampleType,
                        BRAF = patientAnn$BRAFstatus,
                        Gender = patientAnn$Sex,
                        Age = patientAnn$Age,
                        TC = patientAnn$Tumor_content)

hm1<-Heatmap(cuartoRedExp,
             col=col_fun,
             name="log2(Intensity)",
             top_annotation = ta,
             cluster_row_slices = F,
             show_row_dend = T,
             row_names_gp = gpar(fontsize(20), fontface = "bold"),
             show_column_dend = T,
             show_column_names = FALSE,
             rect_gp = gpar(col = "grey30", lwd = 0.5),
             width = unit(8,"in"),
             height = unit(8,"in"),
             row_title = NULL
)
draw(hm1, heatmap_legend_side = "right",annotation_legend_side = "right", merge_legend=T) 