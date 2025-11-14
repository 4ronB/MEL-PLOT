library(shiny)
library(readxl)
library(ggplot2)
library(shinycssloaders)
library(shinymanager)
library(shinyjs)
library(DT)
library(plotly)
library(scales)

#All Data
allClins <- readRDS(file = "AllData.rds")[[1]]
allNames <- readRDS(file = "AllData.rds")[[2]]
infoTab <-  read.table(file = "www/abbreviationList.csv",header = F,sep = ";")
###Cero data
clinicalTableCero <- readRDS(file = "ceroData/ceroData.rds")[[1]] 
expTableCero <- readRDS(file = "ceroData/ceroData.rds")[[2]]
proteinNamesCero <- readRDS(file = "ceroData/ceroData.rds")[[3]]
#Primero data
clinicalTablePrimero <- readRDS(file = "primeroData/primeroData.rds")[[1]] 
expTablePrimero <- readRDS(file = "primeroData/primeroData.rds")[[2]]
proteinNamesPrimero <- readRDS(file = "primeroData/primeroData.rds")[[3]]
##Segundo data
clinicalTableSegundo <- readRDS(file = "segundoData/segundoData.rds")[[1]] 
expTableSegundo <- readRDS(file = "segundoData/segundoData.rds")[[2]]
proteinNamesSegundo <- readRDS(file = "segundoData/segundoData.rds")[[3]]
#Tercero data
clinicalTableTercero <- readRDS(file = "terceroData/terceroData.rds")[[1]] 
expTableTercero <- readRDS(file = "terceroData/terceroData.rds")[[2]]
proteinNamesTercero <- readRDS(file = "terceroData/terceroData.rds")[[3]]
#Cuarto Data
clinicalTableCuarto <- readRDS(file = "cuartoData/cuartoData.rds")[[1]] 
expTableCuarto <- readRDS(file = "cuartoData/cuartoData.rds")[[2]]
proteinNamesCuarto <- readRDS(file = "cuartoData/cuartoData.rds")[[3]]

source("allDataComparer.R")
source("basicPlotter.R")
source("groupComparer.R")
source("CorrelationMaker.R")
source("ceroData/ceroSurvival.R")
source("segundoData/segundoSurvival.R")

devInfo <- tagList(
  br(),
  br(),
  "Developed by: Department of Bioinformatics and Department of Dermatology, Venereology and Dermatooncology at Semmelweis University in collaboration with the European Cancer Moonshot Lund Center and Lund University, Sweden.",
  br(),
  br(),
  img(src = "1200px-Logo_univsemmelweis.svg.png", height = 110, width = 110),
  img(src = "borklinika.png", height = 110, width = 110),
  img(src = "european-cancer-moonshot-logo.png", height = 110, width = 110),
  img(src = "LundUniversity_Logo.png", height = 110, width = 110),
  br(),
  br(),
  "All rights reserved, copyright  © 2021-2025. How to cite: ",
  a("Melanoma Proteomics Unveiled: Harmonizing Diverse Data Sets for Biomarker Discovery and Clinical Insights via MEL-PLOT \n
    Journal of Proteome Research 24 (6), 3117-3128",
    href = "https://pubs.acs.org/doi/full/10.1021/acs.jproteome.4c00749", target = "_blank"),
  br(),
  br()
)

Features_UI <- function(id, clinicalTab) {
  ns <- NS(id)
  tagList(
    fluidPage(fluidRow(selectizeInput(inputId = ns("featureSelect"), 
                                      label = "Choose a clinical feature",
                                      choices = colnames(clinicalTab),
                                      selected = "SampleType"
    )),
    fluidRow(column(plotOutput(ns("featurePlot")),width = 8)),
    fluidRow(tableOutput(ns("featureTable"))),
    fluidRow(textOutput(ns("featureDetail"))),
    br()
    ),
  )
}
Features_Server <- function(id, clinicalTab) {
  moduleServer(
    id,
    function(input, output, session) {
      doFeature <- eventReactive(input$featureSelect,{basicPlotter(selVar = input$featureSelect,InputTab = clinicalTab, infoTab = infoTab)})
      output$featurePlot <- renderPlot(doFeature()$sumPlot)
      output$featureTable <- renderTable(doFeature()$sumTab,rownames = FALSE)
      output$featureDetail <- renderText(doFeature()$dataToDisplay)
      
    }
  )
}

groupComparerUI <- function(id,clinTable, featureNames1,featureNames2, sliderName1) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(title = "Group Comparison ",
                   br(),
                   "By this analysis page the user can compare the expresison profile of the selected protein based on the selected clinical feature.",
                   br(),
                   br(),
                   selectizeInput(inputId = ns("groupProtein"),
                                  label = "Choose a protein",
                                  choices = NULL,
                                  options = list(`live-search` = TRUE,`delimiter` = " ",
                                                 create = T,
                                                 `maxItems` = 10), multiple = T),
                   
                   selectizeInput(inputId = ns("featureSelector"),
                                  label = "Choose a clinical feature for violin plot",
                                  choices = featureNames1,
                                  selected = "Sex"),
                   
                   selectizeInput(inputId = ns("furtherFeatureSelector"),
                                  label = "Choose a clinical feature for jitter",
                                  choices = featureNames2, 
                                  selected = "Age"),
                   
                   sliderInput(ns("AdjustSlider1"), label = paste0( "Adjust to Tumor content"),
                               min = round(min(clinTable[,colnames(clinTable) == sliderName1]), digits = 1),
                               max = round(max(clinTable[,colnames(clinTable) == sliderName1]), digits = 1),
                               value = c(round(min(clinTable[,colnames(clinTable) == sliderName1]), digits = 1),
                                         round(max(clinTable[,colnames(clinTable) == sliderName1]), digits = 1))
                   ),
                   actionButton(inputId = ns("GroupAction"), label = "Start analysis"),
                   devInfo,
                   
                   uiOutput(ns("downloadExpTabUiButton"))
      ),
      mainPanel(
        plotlyOutput(outputId = ns("groupCompPlot")),
        #plotOutput(outputId = ns("groupCompPlot")),
                br()
      )
    )
  )
}
groupComparerServer <- function(id, proteinNames,exptable,clinTable,histVars,downloadTabSource,fileNameSource) {
  moduleServer(
    id,
    function(input, output, session) {
      updateSelectizeInput(session, "groupProtein", choices = proteinNames, server = TRUE, options = list(placeholder = 'Type protein name or Uniprot ID'))
      doProt <- eventReactive(input$GroupAction,{input$groupProtein})
      doGroupComp <- eventReactive(input$GroupAction,{
        groupComparer(proteinName = doProt(),
                      xAxis = input$featureSelector,  
                      jitterVariable = input$furtherFeatureSelector, 
                      histParams = input$AdjustSlider1,
                      histVars = histVars,
                      expTable = exptable, 
                      clinTable = clinTable,
                      proteins = proteinNames
        )})
      
      output$groupCompPlot <- renderPlotly(ggplotly(doGroupComp()$boxPlot,
                                                      #width = 1500,
                                                      height = 650*length(doProt()),
                                                      #width = cdata$output_groupCompPlot_width, #
                                                      #height = cdata$output_groupCompPlot_height, #
                                                      tooltip = c("y", "x", "jitterVariable" )))
      
      #output$groupCompPlot <- renderPlot(doGroupComp()$boxPlot)
        
      
      output$downloadExpTabUiButton <- renderUI({
        ns <- NS(id)
        tagList(
          downloadButton(outputId = ns("groupCompDownload"), label = "Download expression data"))
      })
      output$groupCompDownload <- downloadHandler(filename = function() {paste(input$GroupAction,"_underlyingDataMM500", ".txt", sep = "\t")},
                                                  content = function(file) {write.csv(doGroupComp()$downloadTable, file, row.names = F)})
      
    }
  )
}

Correlation_UI <- function(id, Tissues) {
  ns <- NS(id)
  tagList(sidebarLayout(
    sidebarPanel(
      tabsetPanel(id = ns("Corr"), selected = "Multiple Correlation analysis",
                  tabPanel(title = "Multiple Correlation analysis",
                           br(),
                           "Using the multiple correlation panel, one can get insight to the corerlation profile of the selected protein.",
                           br(),
                           br(),
                           #select protein
                           selectizeInput(inputId = ns("corrProtein"),
                                          label = "Choose a protein",
                                          choices = NULL,
                                          options = list(`live-search` = TRUE)),
                           selectizeInput(inputId = ns("corrTissue"),
                                          label = "Select tissue origin",
                                          choices = Tissues),
                           sliderInput(inputId = ns("corrNumForMultipCorr"),
                                       label = "Set correlation coefficient interval",
                                       value = c(0,0.75), min = -1, max = 1, step = 0.05),
                           
                           actionButton(inputId = ns("CorrAction"), label = "Start analysis"),
                           devInfo
                  ),
                  tabPanel(title = "Single Correlation analysis",
                           br(),
                           "Using the single correlation panel, one can perform correlation analysis using two selected proteins.",
                           br(),
                           br(),
                           #select protein A
                           selectizeInput(inputId = ns("singCorrProteinA"),
                                          label = "Choose protein A",
                                          choices = NULL,
                                          options = list(`live-search` = TRUE)),
                           #select protein B
                           selectizeInput(inputId = ns("singCorrProteinB"),
                                          label = "Choose protein B",
                                          choices = NULL,
                                          options = list(`live-search` = TRUE)),
                           selectizeInput(inputId = ns("singCorrTissue"),
                                          label = "Select tissue origin",
                                          choices = Tissues),
                           actionButton(inputId = ns("singleCorrAction"), label = "Start analysis"),
                           devInfo
                  )
      )
    ),
    
    mainPanel(
      conditionalPanel(condition = paste0("input['", ns("Corr"), "'] == 'Multiple Correlation analysis'"),  
                       dataTableOutput(outputId = ns("corrTab"), width = "60%")),
      conditionalPanel(condition = paste0("input['", ns("Corr"), "'] == 'Single Correlation analysis'"),
                       withSpinner(plotOutput(outputId = ns("singleCorrPlot"), width = "60%")),
                       tableOutput(outputId = ns("singleCorrTab"))
      )
    )
    
  )
  
  )
}
Correlation_Server <- function(id, proteins, clinicalTable, expTable) {
  moduleServer(
    id,
    function(input, output, session) {
      
      
      updateSelectizeInput(session, 'corrProtein', choices = proteins, server = TRUE, options = list(placeholder = 'Type protein name or Uniprot ID'))
      
      doCorr <- eventReactive(input$CorrAction,{multiCorr(proteinName = input$corrProtein,
                                                          origin = input$corrTissue,
                                                          cutOff = input$corrNumForMultipCorr,
                                                          clinicalTable = clinicalTable,
                                                          expTable = expTable,
                                                          proteins = proteins)})
      #tableOutput
      output$corrTab <- DT::renderDataTable({doCorr()}, filter = "none",rownames = FALSE, options= list(bAutoWidth = FALSE))
      
      #single correlation
      updateSelectizeInput(session, 'singCorrProteinA', choices = proteins, server = TRUE, options = list(placeholder = 'Type protein name or Uniprot ID'))
      updateSelectizeInput(session, 'singCorrProteinB', choices = proteins, server = TRUE, options = list(placeholder = 'Type protein name or Uniprot ID'))
      
      doSingCorr <- eventReactive(input$singleCorrAction,{singleCorr(proteinA = input$singCorrProteinA,
                                                                     proteinB = input$singCorrProteinB,
                                                                     origin = input$singCorrTissue,
                                                                     clinicalTable = clinicalTable,
                                                                     expTable = expTable,
                                                                     proteins = proteins)
      })
      output$singleCorrPlot <- renderPlot(doSingCorr()$scatPlot)
      output$singleCorrTab <- renderTable(doSingCorr()$SumStat)
      
    }
  )
}



# User interface ----
ui <-
  fluidPage( 
  titlePanel(title = "Visualization of melanoma proteome"),
  navbarPage(title = "",
             tabPanel(title = "Exploratory analysis",
                        sidebarLayout(
                          sidebarPanel(title = "",width = 3,
                                       br(),
                                       "Using this option one can have a general overview of the expression profile of a selected protein in multiple studies.",
                                       br(),
                                       br(),
                                       selectizeInput(inputId = "groupProteinAll",
                                                      label = "Choose a protein",
                                                      choices = NULL,
                                                      options = list(`live-search` = TRUE)),
                                       br(),
                                       br(),
                                       actionButton(inputId = "groupActionAll", label = "Explore!"),
                                       devInfo
                                       ),
                          mainPanel(plotlyOutput(outputId = "groupPlotAll", width = "100%"),
                                    uiOutput("downloadExpTabUiButton_ExpAn"))
                          )
                      ),
             #Cero
             navbarMenu(title = "Cero study analysis",
                        tabPanel(strong("Features"),
                                 Features_UI(id = "ceroFeature", clinicalTab = clinicalTableCero[,colnames(clinicalTableCero)%in% c("SampleType","Tumor_content","Age","Sex","PFS(months)",
                                                                                                                                    "OS(months)","AJCC8_stage","MelanomaSubType","Clark_level","Breslow_thickness","BRAFstatus","Breslow_stage")]),
                                 includeHTML(path = "www/pilotText.txt")
                        ),
                        tabPanel(strong("Group comparison"),
                                 groupComparerUI(id="ceroGroupComparer",
                                                 clinTable = clinicalTableCero,
                                                 featureNames1 = c("SampleType","Sex","AJCC8_stage","MelanomaSubType", "Clark_level","BRAFstatus", "Breslow_stage"),
                                                 featureNames2 = c("SampleType","Tumor_content","Age","Sex","PFS(months)","OS(months)","AJCC8_stage","MelanomaSubType","Clark_level","Breslow_thickness","BRAFstatus","Breslow_stage"), 
                                                 sliderName1 = "Tumor_content")
                                 ),
                        tabPanel(strong("Survival analysis"),
                                 sidebarLayout(
                                   sidebarPanel(title = "Survival analysis ",
                                                br(),
                                                "Using the survival analysis option one can estimate the effect of the selected protein's expression on the overall survival.",
                                                br(),
                                                br(),
                                                selectizeInput(inputId = "survProteinPilot",
                                                               label = "Choose a protein",
                                                               choices = NULL,
                                                               options = list(`live-search` = TRUE)),
                                                selectizeInput(inputId = "groupSurvSelectorPilot",
                                                               label = "Choose a tissue type",
                                                               choices = c("Lymph Node", "Cutaneous Metastasis", "Primary"), 
                                                               selected = "Primary"),
                                                sliderInput("tcAdjustSliderPilotSurv", label = "Adjust to Tumor content",
                                                            min = 0,
                                                            max = round(max(na.omit(clinicalTableCero$Tumor_content)), digits = 1),
                                                            value = c(0, round(max(na.omit(clinicalTableCero$Tumor_content))))
                                                ),
                                                actionButton(inputId = "survActionPilot", label = "Start analysis"),
                                                devInfo
                                                ),
                                   mainPanel(column(width = 6, plotOutput(outputId = "survPlotPilot", width = "700px", height = "500px"),
                                                    br(),
                                                    br(),
                                                    tableOutput(outputId = "survTabPilot"),
                                                    uiOutput(outputId = "survPilotDataDownload")),
                                             br(),
                                             br(),
                                             br(),
                                             br())
                                 )),
                        tabPanel(strong("Correlation analysis"),
                                 Correlation_UI(id = "ceroCorr", Tissues = c("All", unique(clinicalTableCero$SampleType)))
                                 )
                        ),
             #Primero
             navbarMenu(title = "Primero study analysis",
                        tabPanel(strong("Features"),
                                 Features_UI(id = "primeroFeature", clinicalTab = clinicalTablePrimero[,colnames(clinicalTablePrimero) %in% c("Sex","AJCC8_stage","stage","SampleType","Age","NRAS","BRAFstatus")]),
                                 includeHTML(path = "www/primeroText.txt")
                        ),
                        tabPanel(strong("Group comparsion"),
                                 groupComparerUI(id="primeroGroupComparer",
                                                 clinTable = clinicalTablePrimero,
                                                 featureNames1 = c("Sex","AJCC8_stage","stage","SampleType","NRAS","BRAFstatus"),
                                                 featureNames2 = c("Age","AJCC8_stage","stage","SampleType","Sex","NRAS","BRAFstatus"),
                                                 sliderName1 = "Average_of_Tumor(%)")
                                 ),
                        tabPanel(strong("Correlation analysis"),
                                 Correlation_UI(id = "primeroCorr", Tissues = c("All", unique(clinicalTablePrimero$SampleType)))
                                 )),
             #SEGUNDO
             navbarMenu(title = "Segundo study analysis",
                        tabPanel(strong("Features"),
                                 Features_UI(id = "segundoFeature", clinicalTab = clinicalTableSegundo[,colnames(clinicalTableSegundo) %in% c("Sex","AJCC8_stage","SampleType","Age","Breslow_thickness","Breslow_stage","Clark_level","MelanomaSubtype","OS_DAYS","DistMetLocation","NRAS","BRAFstatus")]),
                                 includeHTML(path = "www/segundoText.txt")
                        ),
                        tabPanel(strong("Group comparsion"),
                                 groupComparerUI(id="segundoGroupComparer",
                                                 clinTable = clinicalTableSegundo,
                                                 featureNames1 =  c("Sex","AJCC8_stage","SampleType","Breslow_stage","MelanomaSubtype","NRAS","BRAFstatus"),
                                                 featureNames2 =  c("Sex","AJCC8_stage","SampleType","Age","Breslow_thickness","Breslow_stage","Clark_level","MelanomaSubtype","OS_DAYS","DistMetLocation","NRAS","BRAFstatus"), 
                                                 sliderName1 = "Tumor_content")
                                 ),
                        tabPanel(strong("Survival analysis"),
                                 sidebarLayout(
                                   sidebarPanel(title = "Survival analysis ",
                                                br(),
                                                "Using the survival analysis option one can estimate the effect of the selected protein's expression on the overall survival.",
                                                br(),
                                                br(),
                                                selectizeInput(inputId = "survProteinSegundo",
                                                               label = "Choose a protein",
                                                               choices = NULL,
                                                               options = list(`live-search` = TRUE)),
                                                sliderInput("tcAdjustSliderSegundoSurv", label = "Adjust to Tumor content",
                                                            min = 0,
                                                            max = round(max(clinicalTableSegundo$Tumor_content), digits = 1),
                                                            value = c(0, round(max(clinicalTableSegundo$Tumor_content)))
                                                ),
                                                actionButton(inputId = "survActionSegundo", label = "Start analysis"),
                                                devInfo
                                                ),
                                   mainPanel(column(width = 6, 
                                                    plotOutput(outputId = "survPlotSegundo", width = "100%"),
                                                    br(),
                                                    br(),
                                                    tableOutput(outputId = "survTabSegundo"),
                                                    uiOutput(outputId = "survSegundoDataDownload")),
                                             br(),
                                             br())
                                 )),
                        tabPanel(strong("Correlation analysis"),
                                 Correlation_UI(id = "segundoCorr", Tissues = c("All", unique(clinicalTableSegundo$SampleType)))
                                 )),
             #Tercero
             navbarMenu(title = "Tercero study analysis",
                        tabPanel(strong("Features"),
                                 Features_UI(id = "terceroFeature", clinicalTab = clinicalTableTercero[,colnames(clinicalTableTercero) %in% c("SampleType","BRAFstatus","NRAS","cKIT","Origin","Sex","Age","OS.months","Primary_localization","Clark_level","Breslow_thickness","MelanomaSubtype")]),
                                 includeHTML(path = "www/terceroText.txt")
                                 ),
                        tabPanel(strong("Group comparsion"),
                                 groupComparerUI(id="terceroGroupComparer",
                                                 clinTable = clinicalTableTercero,
                                                 featureNames1 = c("SampleType","BRAFstatus","NRAS","cKIT","Sex","Primary_localization","Clark_level","MelanomaSubtype","Origin"),
                                                 featureNames2 = c("SampleType","BRAFstatus","NRAS","cKIT","Origin","Sex","Age","OS.months","Primary_localization","Clark_level","Breslow_thickness","MelanomaSubtype"),
                                                 sliderName1 = "Tumor_T")
                                 ),
                        tabPanel(strong("Correlation analysis"),
                                 Correlation_UI(id = "terceroCorr", Tissues = c("All", unique(clinicalTableTercero$SampleType)))
                                 )),
             
             #CUARTO
             navbarMenu(title = "Cuarto study analysis",
                        tabPanel(strong("Features"),
                                 Features_UI(id = "cuartoFeature", clinicalTab = clinicalTableCuarto[,colnames(clinicalTableCuarto) %in% c("Sex","AJCC8_stage","MelanomaSubtype","Clark_level","Breslow_thickness","SampleType","BRAFstatus")]),
                                 includeHTML(path = "www/cuartoText.txt")
                                 
                                 ),
                        tabPanel(strong("Group comparsion"),
                                 groupComparerUI(id="cuartoGroupComparer",
                                                 clinTable = clinicalTableCuarto,
                                                 featureNames1 = c("Sex","AJCC8_stage","MelanomaSubtype","Clark_level","SampleType","BRAFstatus"),
                                                 featureNames2 = c("Age","Sex","AJCC8_stage","MelanomaSubtype","Clark_level","Breslow_thickness","SampleType","BRAFstatus","Tumor_content"),
                                                 sliderName1 = "Tumor_content")
                                 ),
                        tabPanel(strong("Correlation analysis"),
                                 Correlation_UI(id = "cuartoCorr", Tissues = c("All", unique(clinicalTableCuarto$SampleType)))
                                 )
                        ),
             
             tabPanel(title = "List of abbreviations",
                      includeHTML(path = "www/clinDataInfo.txt")
                        )#
             )
  )
# Server logic ----
server <- function(input, output,session){
  #########################################
  ##server functions for ALL data
  updateSelectizeInput(session, 'groupProteinAll', choices = allNames, server = TRUE, options = list(placeholder = 'Type protein name or Uniprot ID'))
  
  
  doGroupAll <- eventReactive(input$groupActionAll,{allPlotter(selProt = input$groupProteinAll,
                                                               allClins = allClins,
                                                               proteinNamesCero = proteinNamesCero,
                                                               proteinNamesPrimero = proteinNamesPrimero,
                                                               proteinNamesSegundo = proteinNamesSegundo, 
                                                               proteinsNamesTercero = proteinNamesTercero, 
                                                               proteinNamesCuarto = proteinNamesCuarto,
                                                               expTableCero = expTableCero, 
                                                               expTablePrimero = expTablePrimero, 
                                                               expTableSegundo = expTableSegundo, 
                                                               expTableTercero = expTableTercero, 
                                                               expTableCuarto = expTableCuarto)})
  output$groupPlotAll <- renderPlotly(ggplotly(doGroupAll()$allSumPlot))
  output$downloadExpTabUiButton_ExpAn <- renderUI({
    downloadButton(outputId = "expAn_dataDownload", label = "Download expression data")
  })
  output$expAn_dataDownload <- downloadHandler(filename = function() {paste("dataMM500", ".txt", sep = "\t")},
                                              content = function(file) {write.csv(doGroupAll()$tabToDownload, file, row.names = F)})
  
  ##server functions for CERO data 
  #features
  Features_Server(id = "ceroFeature", clinicalTab = clinicalTableCero)
  ##GroupComparison 
  
  groupComparerServer(id = "ceroGroupComparer", proteinNames = proteinNamesCero, exptable = expTableCero, 
                      clinTable = clinicalTableCero,histVars = "Tumor_content")
  ##Survival analysis CERO
  updateSelectizeInput(session, 'survProteinPilot', choices = proteinNamesCero, server = TRUE, options = list(placeholder = 'Type protein name or Uniprot ID'))
  
  doSurvAnPilot <- eventReactive(input$survActionPilot,{
    survPlotterCero(variable = input$survProteinPilot,
                     histParams = input$tcAdjustSliderPilotSurv,
                     clinicalTable = clinicalTableCero,
                     expTable = expTableCero,
                     Type = input$groupSurvSelectorPilot
                     )
  })
  
  output$survPilotDataDownload <- renderUI({
    req(input$survActionPilot)
    downloadButton("pilotSegundoDownload", label = "Download expression data") })
  
  output$pilotSegundoDownload <- downloadHandler(filename = function() {paste("ceroSurvDataMM500", ".txt", sep = "\t")},
                                                content = function(file) {write.csv(doSurvAnPilot()$tabToDownload, file, row.names = F)})
  
  
  
  #plot OUtput
  output$survPlotPilot <- renderPlot(doSurvAnPilot()$survPlot)
  output$survTabPilot <- renderTable(doSurvAnPilot()$resTab)
  
  ##correlation CERO
  Correlation_Server(id = "ceroCorr",proteins = proteinNamesCero ,clinicalTable = clinicalTableCero,expTable = expTableCero)
  
  ##########################################
  ##server functions for PRIMERO data
  #features
  Features_Server(id = "primeroFeature", clinicalTab = clinicalTablePrimero)
  ##GroupComparison 
  groupComparerServer(id = "primeroGroupComparer", proteinNames = proteinNamesPrimero, exptable = expTablePrimero, 
                      clinTable = clinicalTablePrimero, histVars = "Average_of_Tumor(%)")
  #correlation
  Correlation_Server(id = "primeroCorr",proteins = proteinNamesPrimero ,clinicalTable = clinicalTablePrimero,expTable = expTablePrimero)
  
  ##########################################
  ##server functions for SEGUNDO data 
  #features
  Features_Server(id = "segundoFeature", clinicalTab = clinicalTableSegundo)
  ##GroupComparison 
  groupComparerServer(id = "segundoGroupComparer", proteinNames = proteinNamesSegundo, exptable = expTableSegundo, 
                      clinTable = clinicalTableSegundo,histVars = "Tumor_content")
  ##Survival SEGUNDO 
  updateSelectizeInput(session, 'survProteinSegundo', choices = proteinNamesSegundo, server = TRUE, 
                       options = list(placeholder = 'Type protein name or Uniprot ID'))
  
  doSurvAnSegundo <- eventReactive(input$survActionSegundo,{
    survPlotterSegundo(variable = input$survProteinSegundo,
                       histParams = input$tcAdjustSliderSegundoSurv,
                       clinicalTableSegundo = clinicalTableSegundo,
                       expTable = expTableSegundo)
  })
  output$survPlotSegundo <- renderPlot(doSurvAnSegundo()$survPlot)
  output$survTabSegundo <- renderTable(doSurvAnSegundo()$resTab)
  
  output$survSegundoDataDownload <- renderUI({
    req(input$survActionSegundo)
    downloadButton("survSegundoDownload", label = "Download expression data") })
  
  output$survSegundoDownload <- downloadHandler(filename = function() {paste("segundoSurvDataMM500", ".txt", sep = "\t")},
                                                content = function(file) {write.csv(doSurvAnSegundo()$tabToDownload, file, row.names = F)})
  
  
  
  #correlation
  Correlation_Server(id = "segundoCorr",proteins = proteinNamesSegundo ,clinicalTable = clinicalTableSegundo,expTable = expTableSegundo)
  ##########################################
  ##server functions for TERCERO data
  #features
  Features_Server(id = "terceroFeature", clinicalTab = clinicalTableTercero)
  # #group comparer
  groupComparerServer(id = "terceroGroupComparer", proteinNames = proteinNamesTercero, exptable = expTableTercero, 
                      clinTable = clinicalTableTercero,histVars = "Tumor_T")
  #correlation
  Correlation_Server(id = "terceroCorr",proteins = proteinNamesTercero ,clinicalTable = clinicalTableTercero,expTable = expTableTercero)
  ##########################################
  ##server functions for CUARTO data 
  #features
  Features_Server(id = "cuartoFeature", clinicalTab = clinicalTableCuarto)
  #group comparer
  groupComparerServer(id = "cuartoGroupComparer", proteinNames = proteinNamesCuarto, exptable = expTableCuarto, 
                      clinTable = clinicalTableCuarto,histVars = "Tumor_content")
  #correlation
  Correlation_Server(id = "cuartoCorr",proteins = proteinNamesCuarto ,clinicalTable = clinicalTableCuarto,expTable = expTableCuarto)
}

# Run app ----
shinyApp(ui = ui, server = server)

