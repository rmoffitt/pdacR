require(pdacR)
library(bioDist)
library(ConsensusClusterPlus)
library(foreign)
library(ggplot2)
library(ggpubr)
library(gplots)
library(limma)
library(nlme)
library(preprocessCore)
library(RColorBrewer)
library(reshape2)
library(Rtsne)
library(shiny)
library(shinyjs)
library(shinythemes)
library(stringr)
library(survival)
library(survminer)
library(scales)
#library(umapr)
library(DESeq2)


print(sessionInfo())

ui <- fluidPage(shinyjs::useShinyjs(),
                titlePanel("PDAC gene expression visualization tool"),
                #----------------------------------------
                fluidRow(tabsetPanel(type = "tabs",selected = "Heatmap",id = "bigTab",
                                     tabPanel(title = "Heatmap",
                                              fluidRow(column(2, offset = 9,
                                                              fluidRow(uiOutput("ReactiveBarplotcolor")))
                                              ),
                                              fluidRow(column(8,plotOutput(outputId="heatmap",
                                                                           height = "700px")),
                                                       column(4,plotOutput(outputId = "barplot",
                                                                           height = "700px")))
                                     ),
                                     tabPanel(title = "Cartesian",
                                              fluidRow(column(2, offset = 9,
                                                              fluidRow(uiOutput("ReactiveBarplotcolor2")))
                                              ),
                                              fluidRow(fluidRow(column(9, plotOutput(outputId="cartesian",
                                                                                     height = "700px")),
                                                                column(3, plotOutput(outputId = "barplot2",
                                                                                     height = "700px"))),
                                                       fluidRow(column(3),
                                                                column(2,radioButtons(inputId = "projection",
                                                                                      label = "X Y Projection",
                                                                                      #removed UMAP
                                                                                      choices = c("tSNE","PCA"),
                                                                                      selected = "tSNE") ),
                                                                column(4,radioButtons(inputId = "painting",
                                                                                      label = "How points should be colored",
                                                                                      choices = c("Signatures (R+G+B)","Categorical"),
                                                                                      selected = "Categorical") )
                                                       )
                                              )
                                     ),
                                     tabPanel(title = "Survival",
                                              fluidRow(column(10,plotOutput(outputId = "survival",
                                                                            height = "700px")),
                                                       column(2,
                                                              fluidRow(uiOutput("ReactiveContinuousSurv.1")),
                                                              fluidRow(uiOutput("ReactiveContinuousSurv.2")))),
                                              fluidRow(column(2, radioButtons(inputId = "survivalType",
                                                                              label = "Method",
                                                                              choices = c("Kaplan Meier", "Cox regression"),
                                                                              selected = "Kaplan Meier")),
                                                       column(2, uiOutput("ReactiveSurvival")),
                                                       column(4, htmlOutput(outputId = "ReactiveSurvivaltwo"))
                                              )
                                     ),
                                     tabPanel(title = "Diff Expr",
                                              fluidRow(column(9,plotOutput(outputId = "volcano",
                                                                           height = "700px")),
                                                       column(1,
                                                              fluidRow(uiOutput("DESeqcontrastA"))),
                                                       column(1,
                                                              fluidRow(uiOutput("DESeqcontrastB"))),
                                                       column(1,
                                                              fluidRow(uiOutput("Experiment_Type"))),
                                                       column(2,
                                                              fluidRow(actionButton(inputId = "go",
                                                                                    label = "Run Diff Expr"))))
                                     )
                )

                ),
                #----------------------------------------
                # fluidRow( column(2,downloadLink(label = "save image as pdf", outputId = 'pdflink') ),
                #           column(2,downloadLink(label = "save data as table", outputId = 'tablelink') )),
                #----------------------------------------
                fluidRow(
                  column(3,
                         radioButtons(inputId = "Species",
                                      label = "Select Species",
                                      choices = c("Human", "Mouse"),
                                      selected = "Human"),
                         uiOutput("ReactiveDatasets"),
                         htmlOutput(outputId = "ReactiveCitation"),
                         textInput(inputId = "addPackageText",
                                   label = "Add private data sets",
                                   value = "",
                                   placeholder="Name of data"),

                         textOutput("text"),

                         actionButton(inputId = "addPackageButton",
                                      label="load"),
                         uiOutput("ReactiveSampleFilterSpecifics")
                  ),
                  column(2,
                         radioButtons(inputId = "Defaults",
                                      label = "Use suggested filters",
                                      choices = c("Yes", "No"),
                                      selected = "Yes"),
                         uiOutput("ReactiveSampleFilters")
                  ),
                  column(2,
                         fluidRow(uiOutput("SelectionButton"),
                                  uiOutput("ClearButton")),
                         fluidRow(htmlOutput(outputId = "GeneTitle")),
                         fluidRow(uiOutput("ReactiveGenesets"))
                  ),
                  column(2,
                         fluidRow(
                           column(2,radioButtons(inputId = "scaling",
                                                 label = "Scale Colors",
                                                 choices = c("Row","None"),
                                                 selected = "Row")
                           ),
                           column(1,offset = 2,  checkboxGroupInput(inputId = "remove.zeros",
                                                                    label = "Invariant genes",
                                                                    choices = c("hide"),
                                                                    selected = c("hide"))),
                         ),
                         fluidRow(radioButtons(inputId = "sampleClustertype",
                                               label = "Sample Method",
                                               choices = c("Consensus","Consensus.RowScaled","Euclidean","Pearson","K.Means","Sorted"),
                                               selected = "Pearson") ),
                         fluidRow(uiOutput("ReactiveSampleSortBy") ),
                         fluidRow(radioButtons(inputId = "geneClustertype",
                                               label = "Gene Method",
                                               choices = c("Consensus","Euclidean","Pearson","Sorted","K.Means"),
                                               selected = "Pearson") ),
                         fluidRow(uiOutput("ReactiveGeneSortBy") ),
                         fluidRow(textInput(inputId = "userGeneList.1",
                                            label = "User selected genes 1",
                                            value = "COL17A1, TFF1, KRT17 EPCAM") ),
                         fluidRow(textInput(inputId = "userGeneList.2",
                                            label = "User selected genes 2",
                                            value = "PDCD1") ),
                         fluidRow(textInput(inputId = "userGeneList.3",
                                            label = "User selected genes 3",
                                            value = "CD274") )),
                  column(2, fluidRow(uiOutput("ReactiveXaxisLabel")),
                         fluidRow(uiOutput("ReactiveSignatureTrack.1")),
                         fluidRow(uiOutput("ReactiveSignatureTrack.2")),
                         fluidRow(uiOutput("ReactiveSignatureTrack.3")),
                         fluidRow(uiOutput("ReactiveSampleTracks")))
                ))

server <- function(input, output) {

  # ================================================================
  # ------------- Reactive UI objects ------------------------------
  # ================================================================
  # Which data sets are available
  output$ReactiveDatasets <- renderUI({
    radioButtons(inputId = "dataset",
                 label = "Data sets to use",
                 choices = globals$data_set_list$labels,
                 selected = "TCGA PAAD, 2017")
  })
  # ==========================================================
  # What is the source citation for the dataset currently being used?
  output$ReactiveCitation <- renderText({
    paste("<B>You are currently using data from:</B>",
          dataSet()$metadata$reference,
          "",
          "<B>Dataset accession:</B>",
          dataSet()$metadata$accession,
          "",
          sep = "<br/>")


  })


  # ==========================================================
  # Which gene sets are available
  output$GeneTitle <- renderText({
    paste("<B>Gene sets to use</B>")
  })

  shinyjs::runjs("document.getElementById('reset1').style.visibility = 'hidden';")

  observeEvent(input$bigTab, {
    choice = input$bigTab
    if(choice == "Heatmap"){
      observeEvent(input$GeneSelection, {
        shinyjs::runjs("document.getElementById('GeneSelection').style.visibility = 'hidden';")
        shinyjs::runjs("document.getElementById('reset1').style.visibility = 'visible';")
      })
      observeEvent(input$reset1, {
        shinyjs::runjs("document.getElementById('GeneSelection').style.visibility = 'visible';")
        shinyjs::runjs("document.getElementById('reset1').style.visibility = 'hidden';")
      })

      r$colorTab = "heatmap"
    }
    else {
      shinyjs::runjs("document.getElementById('GeneSelection').style.visibility = 'hidden';")
      shinyjs::runjs("document.getElementById('reset1').style.visibility = 'hidden';")
      r$colorTab = "other"
    }
  })


  output$SelectionButton <- renderUI({
    actionButton("GeneSelection", "Generate Heatmap")
  })

  observe({
    if(r$heatmapStoplight == "red"){
      output$ClearButton <- NULL

    } else{
      output$ClearButton <- renderUI({
        actionButton("reset1", "Clear")
      })

    }
  })


  output$ReactiveGenesets <- renderUI({
    meta = dataSet()$meta
    if(exists("default_selections", where = meta) &&
       exists("geneSets", where = meta$default_selections)){
      checkboxGroupInput(inputId = "genesets",
                         label = "",
                         choices = c("Most Variable 100",
                                     "Most Variable 1000",
                                     "User selected genes 1",
                                     "User selected genes 2",
                                     "User selected genes 3",
                                     gsub(pattern = ".",
                                          replacement = " ",
                                          fixed = TRUE,
                                          x = names(globals$gene_lists))),
                         selected = gsub(pattern = ".",
                                         replacement = " ",
                                         fixed = TRUE,
                                         x = meta$default_selections$geneSets))
      # shinyjs::runjs("Shiny.setInputValue(genesets.2, input$genesets);")
      # r$heatmapStoplight = "green"
    } else {
      checkboxGroupInput(inputId = "genesets",
                         label = "",
                         choices = c("Most Variable 100",
                                     "Most Variable 1000",
                                     "User selected genes 1",
                                     "User selected genes 2",
                                     "User selected genes 3",
                                     gsub(pattern = ".",
                                          replacement = " ",
                                          fixed = TRUE,
                                          x = names(globals$gene_lists))))
    }
  })

  observeEvent(input$bigTab, {
    choice = input$bigTab
    if(choice != "Heatmap"){
      output$ReactiveGenesets <- renderUI({
        checkboxGroupInput(inputId = "genesets",
                           label = "",
                           choices = c("Most Variable 100",
                                       "Most Variable 1000",
                                       "User selected genes 1",
                                       "User selected genes 2",
                                       "User selected genes 3",
                                       gsub(pattern = ".",
                                            replacement = " ",
                                            fixed = TRUE,
                                            x = names(globals$gene_lists))),
                           selected = input$genesets)
      })
    }
  })


  getGeneSets = reactive({
    pretty.selection <- input$genesets
    sacred.set <- c("Most Variable 100",
                    "Most Variable 1000",
                    "User selected genes 1",
                    "User selected genes 2",
                    "User selected genes 3")
    part1 <- intersect(pretty.selection,sacred.set)
    part2 <- setdiff(pretty.selection,sacred.set)
    machine.friendly <- gsub(x = part2,
                             pattern = " ",
                             replacement = ".",
                             fixed = TRUE)
    return(c(part1,machine.friendly))
  })

  r <- reactiveValues(heatmapStoplight = "red")

  output$text <- renderText({
    r$heatmapStoplight
  })

  observeEvent(input$GeneSelection, {
    output$ReactiveGenesets <- renderUI({
      checkboxGroupInput(inputId = "genesets.2",
                         label = "Select genesets to retain",
                         choices = c(input$genesets),
                         selected = c(input$genesets))
    })
    r$heatmapStoplight = "green"
  })


  observeEvent(input$reset1,{
    output$ReactiveGenesets <- renderUI({
      checkboxGroupInput(inputId = "genesets",
                         label = "",
                         choices = c("Most Variable 100",
                                     "Most Variable 1000",
                                     "User selected genes 1",
                                     "User selected genes 2",
                                     "User selected genes 3",
                                     gsub(pattern = ".",
                                          replacement = " ",
                                          fixed = TRUE,
                                          x = names(globals$gene_lists))),
                         selected = input$genesets)
    })
    r$heatmapStoplight = "red"
    output$heatmap <- NULL
  })



  # ==========================================================
  # What survival data is in the dataset
  output$ReactiveSurvivaltwo <- renderText({

    if(dataSet()$metadata$survivalA %in% "None"){
      paste("** This dataset does not have survival data **")
    } else {
      paste("You are currently using < ",
            dataSet()$metadata[which(dataSet()$metadata %in% input$survivalDays)],
            " > for your survival analysis",
            sep = "")
    }

  })


  output$ReactiveSurvival <- renderUI({

    if((!dataSet()$metadata$survivalA %in% "None") &&
       (dataSet()$metadata$survivalB %in% "None")){
      radioButtons(inputId = "survivalDays",
                   label = "Survival Factors",
                   choices = c(dataSet()$metadata[grep(dataSet()$metadata,
                                                       pattern = "survival")[1]]))
    } else if((!dataSet()$metadata$survivalA %in% "None") &&
              (!dataSet()$metadata$survivalB %in% "None")) {
      radioButtons(inputId = "survivalDays",
                   label = "Survival Factors",
                   choices = c(dataSet()$metadata[grep(dataSet()$metadata,
                                                       pattern = "survival")]))
    } else {

    }
  })
  # ==========================================================
  # How continuous variables should be split for survival analysis

  output$ReactiveContinuousSurv.1 <- renderUI({
    validate(need(length(getSampleTracks()) > 0, "Please select a first sample track"))

    x <- dataSet()
    s = sampleSet()
    x$sampInfo = x$sampInfo[s,,drop=T]
    tmp.info <- getSignatures()
    x$sampInfo <- cbind(x$sampInfo,tmp.info)

    if((length(levels(as.factor(x$sampInfo[, getSampleTracks()[1]]))) >= 5) &&
       (is.numeric(x$sampInfo[, getSampleTracks()[1]]))){
      sliderInput(inputId = "factorized1",
                  label = paste(getSampleTracks()[1]),
                  min = 0,
                  max = 1,
                  value = 0.50,
                  step = 0.05)
    } else{
      radioButtons(inputId = "placeholder1",
                   label = paste(getSampleTracks()[1]),
                   choices = "This factor is not continuous")
    }
  })

  output$ReactiveContinuousSurv.2 <- renderUI({
    validate(need((length(getSampleTracks()) > 1) &&
                    (length(getSampleTracks()) < 3), "Please select no more than a second sample track"))
    x <- dataSet()
    s = sampleSet()
    x$sampInfo = x$sampInfo[s,,drop=T]
    tmp.info <- getSignatures()
    x$sampInfo <- cbind(x$sampInfo,tmp.info)

    if((length(levels(as.factor(x$sampInfo[, getSampleTracks()[2]]))) >= 5) &&
       (is.numeric(x$sampInfo[, getSampleTracks()[2]]))){
      sliderInput(inputId = "factorized2",
                  label = paste(getSampleTracks()[2]),
                  min = 0,
                  max = 1,
                  value = 0.50,
                  step = 0.05)
    } else{
      radioButtons(inputId = "placeholder2",
                   label = paste(getSampleTracks()[2]),
                   choices = "This factor is not continuous")
    }
  })
  # ==========================================================
  # Which sample filter categories should be available
  output$ReactiveSampleFilters <- renderUI({
    meta <- dataSet()$metadata
    possibleTracks <- gsub(pattern = ".", replacement = " ",fixed = TRUE, x = names(dataSet()$sampInfo))
    if(exists("default_selections", where = meta) && input$Defaults == "Yes"){
      checkboxGroupInput(inputId = "sampleFilters",
                         label = "Filter samples by",
                         choices = possibleTracks,
                         selected = meta$default_selections$filter_column)
    } else {
      checkboxGroupInput(inputId = "sampleFilters",
                         label = "Filter samples by",
                         choices = possibleTracks)
    }
  })

  # ==========================================================
  # Which sample sub-categories are actually used?
  output$ReactiveSampleFilterSpecifics <- renderUI({
    sampInfo <- dataSet()$sampInfo
    meta <- dataSet()$metadata
    possibleTracks <- NULL
    for(i in gsub(pattern = " ", replacement = ".",fixed = TRUE, x = input$sampleFilters)){
      if(length(levels(as.factor(sampInfo[,names(sampInfo)==i])))>1){
        possibleTracks <- c(possibleTracks,paste(i,levels(as.factor(sampInfo[,names(sampInfo)==i])), sep=":"))
      }
    }

    if(!is.null(input$sampleFilters) & class(sampInfo[,names(sampInfo)==i]) != "numeric"){
      if(exists("default_selections", where = meta) && input$Defaults == "Yes"){
        # print("Filter by the following")
        # print(meta$default_selections$filter_levels)
        checkboxGroupInput(inputId = "sampleFilterSpecifics",
                           label = "Select to Remove",
                           choices = possibleTracks,
                           selected = meta$default_selections$filter_levels)
      } else {
        # print("Filter options are")
        # print(possibleTracks)
        checkboxGroupInput(inputId = "sampleFilterSpecifics",
                           label = "Select to Remove",
                           choices = possibleTracks)
      }
    }
    else if(class(sampInfo[,names(sampInfo)==i]) == "numeric"){
      radioButtons(inputId = "placeholder3",
                   label = "Select to Remove",
                   choices = "Cannot filter by continuous variables")
    }
  })

  # ==========================================================
  # Which sample tracks should be selectable for visualization
  output$ReactiveSampleTracks <- renderUI({
    meta = dataSet()$metadata
    possibleTracks <- names(dataSet()$sampInfo)
    if(exists("default_selections", where = meta) &&
       exists("sampleTracks", where = meta$default_selections)){
      checkboxGroupInput(inputId = "sampleTracks",
                         label = "Sample Tracks",
                         choices =  gsub(pattern = ".",
                                         replacement = " ",
                                         fixed = TRUE,
                                         x =
                                           c("Expression.signature.1",
                                             "Expression.signature.2",
                                             "Expression.signature.3",
                                             as.character(globals$classifier_list$labels),
                                             "purIST_2019_Call",
                                             "molgrad_PDX",
                                             "molgrad_Puleo",
                                             "molgrad_ICGCarray",
                                             "molgrad_ICGCrnaseq",
                                             possibleTracks)),
                         selected = gsub(pattern = ".",
                                         replacement = " ",
                                         fixed = TRUE,
                                         x =meta$default_selections$sampleTracks))
    } else {
      checkboxGroupInput(inputId = "sampleTracks",
                         label = "Sample Tracks",
                         choices =  gsub(pattern = ".",
                                         replacement = " ",
                                         fixed = TRUE,
                                         x =
                                           c("Expression.signature.1",
                                             "Expression.signature.2",
                                             "Expression.signature.3",
                                             as.character(globals$classifier_list$labels),
                                             "purIST_2019_Call",
                                             "molgrad_PDX",
                                             "molgrad_Puleo",
                                             "molgrad_ICGCarray",
                                             "molgrad_ICGCrnaseq",
                                             possibleTracks)))
    }
  })
  getSampleTracks = reactive({
    pretty.selection <- input$sampleTracks
    machine.friendly <- gsub(x = pretty.selection,
                             pattern = " ",
                             replacement = ".",
                             fixed = TRUE)
    return(machine.friendly)
  })
  # ==========================================================
  # Which tracks to use for x axis label
  output$ReactiveXaxisLabel <- renderUI({
    possibleTracks <- names(dataSet()$sampInfo)
    selectInput(inputId = "XaxisLabel",
                label = "X axis label (heatmap) -or- Color (cartesian)",
                choices = possibleTracks)
  })
  # ==========================================================
  # Which tracks to use for Barplot coloring
  output$ReactiveBarplotcolor <- renderUI({
    possibleTracks <- names(dataSet()$sampInfo)
    selectInput(inputId = "Barplotcolor",
                label = "Box/Scatter Color",
                choices = c("No Color", possibleTracks),
                selected = "No Color")
  })
  output$ReactiveBarplotcolor2 <- renderUI({
    possibleTracks <- names(dataSet()$sampInfo)
    selectInput(inputId = "Barplotcolor2",
                label = "Box/Scatter Color",
                choices = c("No Color", possibleTracks),
                selected = "No Color")
  })

  getColor = reactive({
    if(r$colorTab == "heatmap"){
      pretty.selection <- input$Barplotcolor
    }
    else if(r$colorTab == "other") {
      pretty.selection <- input$Barplotcolor2
    }
    machine.friendly <- gsub(x = pretty.selection,
                             pattern = " ",
                             replacement = ".",
                             fixed = TRUE)
    return(machine.friendly)
  })
  # ==========================================================
  # Expression scores as tracks and colors
  output$ReactiveSignatureTrack.1 <- renderUI({
    possibleTracks <- c("User selected genes 1","User selected genes 2","User selected genes 3",names(globals$gene_lists))
    selectInput(inputId = "SignatureTrack.1",
                label ="Expression signature 1 (Red)",
                choices = possibleTracks,
                selected = "User selected genes 1")
  })
  output$ReactiveSignatureTrack.2 <- renderUI({
    possibleTracks <- c("User selected genes 1","User selected genes 2","User selected genes 3",names(globals$gene_lists))
    selectInput(inputId = "SignatureTrack.2",
                label ="Expression signature 2 (Green)",
                choices = possibleTracks,
                selected = "User selected genes 2")
  })
  output$ReactiveSignatureTrack.3 <- renderUI({
    possibleTracks <- c("User selected genes 1","User selected genes 2","User selected genes 3",names(globals$gene_lists))
    selectInput(inputId = "SignatureTrack.3",
                label ="Expression signature 3 (Blue)",
                choices = possibleTracks,
                selected = "User selected genes 3")
  })
  # ==========================================================
  # What can we sort samples by
  output$ReactiveSampleSortBy <- renderUI({
    if(input$sampleClustertype == "Sorted"){
      possibleTracks <- c(names(getSignatures()),names(dataSet()$sampInfo))
      selectInput(inputId = "sampleSortBy",
                  label = NULL,
                  choices = c("expression",possibleTracks),
                  selected = "expression")
    } else if (input$sampleClustertype %in% c("Consensus.RowScaled","Consensus","K.Means")){
      sliderInput(inputId = "sampleSortBy",
                  label = "Number of Clusters",
                  min=2,max=10,value=2)
    }
  })

  # ==========================================================
  # What can we sort genes by
  output$ReactiveGeneSortBy <- renderUI({
    if(input$geneClustertype == "Sorted"){
      selectInput(inputId = "geneSortBy",
                  label = NULL,
                  choices = c("expression","gene sets"),
                  selected = "gene sets")
    } else if(input$geneClustertype %in% c("Consensus.RowScaled","Consensus","K.Means")){
      sliderInput(inputId = "geneSortBy",
                  label = "Number of Clusters",
                  min=2,max=10,value=2)
    }
  })
  # ==========================================================
  # Radio buttons for DESeq tab
  output$DESeqcontrastA <- renderUI({
    validate(
      need(length(input$dataset) > 0, "Please select a dataset, one sample track, one geneset, and an experiment type")
    )
    validate(
      need((2 > length(getSampleTracks())) && length(getSampleTracks()) > 0, "Please select one sample track from below")
    )

    sampInfo <- dataSet()$sampInfo
    possibleTracks.d <- sampInfo[[which(names(sampInfo) %in% getSampleTracks())]]

    if (class(possibleTracks.d) != "factor"){
      possibleTracks.d <- as.factor(possibleTracks.d)
    }
    validate(
      need(length(levels(possibleTracks.d)) > 0, "Sorry, no options. Please select a different track")
    )
    radioButtons(inputId = "contrastA",
                 label = "A (numerator)",
                 choices = levels(possibleTracks.d))
  })
  output$DESeqcontrastB <- renderUI({
    validate(
      need((2 > length(getSampleTracks())) && length(getSampleTracks()) > 0, "")
    )
    sampInfo <- dataSet()$sampInfo
    possibleTracks.d <- sampInfo[[which(names(sampInfo) %in% getSampleTracks())]]

    if (class(possibleTracks.d) != "factor"){
      possibleTracks.d <- as.factor(possibleTracks.d)
    }
    validate(
      need(length(levels(possibleTracks.d)) > 0, "")
    )

    radioButtons(inputId = "contrastB",
                 label = "B (denominator)",
                 choices = levels(possibleTracks.d))
  })
  # New button for exp type
  output$Experiment_Type <- renderUI({
    validate(
      need(length(input$dataset) > 0, "")
    )

    meta <- getX()$metadata

    if(!(exists("exp.type", where = meta))){
      radioButtons(inputId = "exp.type",
                   label = "Experiment Type",
                   choices = c("scRNA",
                               "RNAseq",
                               "Array"),
                   selected = character(0))
    }
    else{
      radioButtons(inputId = "exp.type",
                   label = "Experiment Type",
                   choices = c("scRNA",
                               "RNAseq",
                               "Array"),
                   selected = meta$exp.type)
    }
  })
  # ================================================================
  # ------------- Reactive Values -------------------------------
  # ================================================================
  ## Dataset location
  globals <- reactiveValues(data_set_list = readRDS("/data/data_set_list.rds"),
                            gene_lists = readRDS("/data/gene_lists.rds"),
                            classifier_list = pdacR::classifier_list,
                            res = NULL,
                            genes2color = NULL)

  observeEvent(input$Species, {
    this.species <- isolate(input$Species)
    if(this.species == "Mouse"){
      globals$data_set_list <- readRDS("/data/mouse_data_set_list.rds")
      globals$gene_lists <- pdacR::mouse_gene_lists
      globals$classifier_list = NULL
    }
    else if (this.species == "Human"){
      globals$data_set_list = readRDS("/data/data_set_list.rds")
      globals$gene_lists = readRDS("/data/gene_lists.rds")
      globals$classifier_list = pdacR::classifier_list
    }
  }
  )

  # ================================================================
  # ------------- Event Observations -------------------------------
  # ================================================================

  observeEvent(input$addPackageButton, {
    new.package <- isolate(input$addPackageText)
    loading.result <- require(new.package,character.only=TRUE)
    if(loading.result==TRUE){
      showNotification(ui = paste("Successfully loaded",new.package),
                       type="message")
      #------------
      # need to avoid collision of variable names when loading new gene lists
      get.new.gene_list <- local( function(new.package = new.package){
        #gene_lists <- NULL
        data(package = new.package,
             list = "gene_lists")
        new.gene_lists <- gene_lists
        return(new.gene_lists)
      })
      new.gene_lists <- get.new.gene_list(new.package)
      # need to avoid collision of variable names when loading new data sets
      get.new.data_set_list <- local( function(new.package = new.package){
        #data_set_list <- NULL
        data(package = new.package,
             list = as.character("data_set_list"))
        new.data_set_list <- data_set_list
        return(new.data_set_list)
      })
      new.data_set_list <- get.new.data_set_list(new.package)
      data(list = as.character(new.data_set_list[,2]))
      #------------
      cat('before action')
      if(!is.null(new.gene_lists)){
        if(length(new.gene_lists)>0){
          globals$gene_lists <- c(globals$gene_lists, new.gene_lists)
          globals$new.gene_lists <- new.gene_lists
          globals$gene_lists <- globals$gene_lists[!duplicated(globals$gene_lists)]
        }
      }
      if(!is.null(new.data_set_list)){
        if(length(new.data_set_list)>0){
          globals$data_set_list <- rbind(globals$data_set_list,new.data_set_list)
          globals$mouse_data_set_list <- rbind(globals$mouse_data_set_list,new.data_set_list)
          globals$data_set_list <- globals$data_set_list[!duplicated.data.frame(globals$data_set_list),]
        }
      }
    }
    if(loading.result==FALSE){
      showNotification(ui = paste("Could not find a package called :",new.package),
                       type="error")
    }
  })

  # observe go button to run differential expression
  observeEvent(input$go, {
    validate(
      need((length(getSampleTracks() == 1)) && (length(getGeneSets() == 1)), "Please select 1 sample track and 1 gene set from below"))

    showNotification(ui = "Running Differential expression, this may take a while",
                     duration = NULL,
                     closeButton = F,
                     type = "message",
                     id = "DESeq_notification")

    tracks <- getSampleTracks()

    x <- isolate(dataSet())
    coldata <- x$sampInfo
    tracks <- coldata[names(coldata) == tracks]

    cts = x$ex
    if (anyNA(tracks[[1]])){
      print(". . . There are NA's in this track, trimming . . .")
      not.NA <- !is.na(tracks[[1]])
      print("dim coldata before and after")
      print(dim(coldata))
      coldata <- coldata[not.NA,]
      print(dim(coldata))
      print("dim cts before and after")
      print(dim(cts))
      cts <- cts[,not.NA]
      print(dim(cts))
    } else {print(". . . There are no NA's, moving to Diff Expr . . .")}


    # make coldata compatible with DESeq design format
    replace <- which(colnames(coldata) %in% getSampleTracks())
    colnames(coldata)[replace] <- "use"

    ## Remove samples not within the A vs B comparison
    comparisons = which(coldata$use %in% c(input$contrastA,input$contrastB))
    coldata = droplevels(coldata[comparisons,])
    cts = cts[,comparisons]

    if(input$exp.type == "scRNA"){
      # quick and dirty b/c its exploratory
      change <- cbind(numerator = rowMeans(cts[,which(coldata$use == input$contrastA)],na.rm=T),
                      denominator = rowMeans(cts[,which(coldata$use == input$contrastB)],na.rm=T))

      test.res <- numeric(length = nrow(cts))
      for(i in 1:nrow(cts)){
        test.res[i] <- t.test(cts[i,which(coldata$use == input$contrastA)],
                              cts[i,which(coldata$use == input$contrastB)])$p.value
      }

      globals$res <- data.frame(labels = x$featInfo$SYMBOL,
                                padj = test.res,
                                log2FoldChange = log2(change$numerator/change$denominator),
                                row.names = x$featInfo$SYMBOL
      )
    }
    else if(input$exp.type == "Array"){
      ## Expects log normalized data, no change
      rownames(cts) <- x$featInfo$SYMBOL

      ## Relevel to make contrastB the intercept
      coldata$use = relevel(coldata$use, ref = input$contrastB)
      design <- model.matrix(~ coldata$use)
      fit <- lmFit(cts, design)
      fit <- eBayes(fit)
      globals$res <- topTable(fit, n = nrow(fit))
      print(summary(globals$res))
      globals$res = globals$res[,c("logFC","adj.P.Val")]
      colnames(globals$res) <- c("log2FoldChange", "padj")
      globals$res$labels <- x$featInfo$SYMBOL[as.numeric(rownames(globals$res))]

    } else {
      ## Undo lognormalization below, easier than replicating parsing (for now)
      cts <- 2^(cts)-1
      cts = apply(cts,2,FUN = function(x){round(x) %>% as.integer()})

      dds <- DESeqDataSetFromMatrix(countData = cts,
                                    colData = coldata,
                                    design = ~ use)
      dds <- DESeq(dds)

      globals$res = as.data.frame(results(dds, contrast = c("use",input$contrastA,input$contrastB)))
      globals$res$labels = x$featInfo$SYMBOL[as.numeric(rownames(globals$res))]
      #print(summary(globals$res))
    }



    this.set <- isolate(input$genesets)
    if("User selected genes 1" %in% input$genesets){
      set <- strsplit(input$userGeneList.1,"[,;[:space:]]+")[[1]]
      print(paste0(". . . Searching for ", set, " in dataset . . ."))
    }else if("User selected genes 2" %in% input$geneset){
      set <- strsplit(input$userGeneList.2,"[,;[:space:]]+")[[1]]
      print(paste0(". . . Searching for ", set, " in dataset . . ."))
    }else if("User selected genes 3" %in% input$genesets){
      set <- strsplit(input$userGeneList.3,"[,;[:space:]]+")[[1]]
      print(paste0(". . . Searching for ", set, " in dataset . . ."))
    }else{
      print(paste0(". . . Searching for ", this.set, " in dataset . . ."))
      set <- gsub(x = this.set,
                  pattern = " ",
                  replacement = ".",
                  fixed = TRUE)
      set <- globals$gene_lists[names(globals$gene_lists) == set]
      set <- set[[1]]
    }

    globals$genes2color = NULL

    for (i in globals$res$labels){
      # print(i)
      if (i %in% set){
        print(paste(i, "is present"))
        globals$genes2color <- c(globals$genes2color,i)
      }
    }
    print("--------- Genes in geneset and data ----------")
    print(globals$genes2color)
  })

  #observe genesets for coloring volcano
  observeEvent(input$genesets, {
    validate(
      need(!is.null(globals$res), "Please run Diff Expr first -- Click the button"))
    print(". . . Updating Gene lists . . .")
    this.set <-isolate(input$genesets)
    print(". . . this.set . . .")
    print(this.set)
    # =============================================
    # add an if here for "user selected genes"
    # =============================================
    if("User selected genes 1" %in% input$genesets){
      set <- strsplit(input$userGeneList.1,"[,;[:space:]]+")[[1]]
      print(paste0(". . . Searching for ", set, " in dataset . . ."))
      print(set)
    }else if("User selected genes 2" %in% input$geneset){
      set <- strsplit(input$userGeneList.2,"[,;[:space:]]+")[[1]]
      print(paste0(". . . Searching for ", set, " in dataset . . ."))
    }else if("User selected genes 3" %in% input$genesets){
      set <- strsplit(input$userGeneList.3,"[,;[:space:]]+")[[1]]
      print(paste0(". . . Searching for ", set, " in dataset . . ."))
    }else{
      print(paste0(". . . Searching for ", this.set, " in dataset . . ."))
      set <- gsub(x = this.set,
                  pattern = " ",
                  replacement = ".",
                  fixed = TRUE)
      set <- globals$gene_lists[names(globals$gene_lists) == set]
      set <- set[[1]]
      print(set)
    }

    globals$genes2color = NULL
    for (i in globals$res$labels){
      if (i %in% set){
        print(paste(i, "is present"))
        globals$genes2color[i] <- i
      }
    }
    print("--------- Genes in geneset and data ----------")
    print(globals$genes2color)
  })


  # ================================================================
  # ------------- Reactive data ------------------------------------
  # ================================================================

  # load full data sets
  dataSet <- reactive({
    writeLines('----------')
    print(". . . Running dataSet . . .")
    validate(
      need(length(input$dataset) > 0, "Please select a data set")
    )
    x <- NULL
    for(selectedset in input$dataset){
      selectedvariable <- globals$data_set_list$variablenames[globals$data_set_list$labels %in% selectedset]
      filename = paste0("/data/",selectedvariable,".rds")
      y = readRDS(file = filename)
      y$sampInfo$source <- selectedset
      do.a.log.transform = TRUE # the default behavior

      if(typeof(y$ex) %in% "S4"){
        #print("Converted from S4 to a regular matrix?")
        #print(typeof(y$ex))
        y$ex <- as.matrix(y$ex)
        #print(typeof(y$ex))
        #print(dim(y$ex))
      }
      if(exists('metadata',where=y)){
        if(exists('log.transformed',where=y$metadata)){
          if(y$metadata$log.transformed==TRUE){
            do.a.log.transform = FALSE
          }
        }
      }
      if(do.a.log.transform){
        y$ex <- log2(1+y$ex)
      }
      x <- mergeDataSets(x,y)
    }
    if(object.size(x)<(10^8.5)){
      x$ex <- limma::normalizeQuantiles(as.matrix(x$ex))
    }
    return(x)
  })

  # ==========================================================
  # load and merge feature sets
  featureSet <- reactive({
    print(". . . Running featureSet . . .")
    validate(
      need(length(getGeneSets()) > 0, "Please select a gene set")
    )
    x <- dataSet()
    s <- sampleSet()
    f <- as.numeric()
    g <- data.frame(SYMBOL=x$featInfo$SYMBOL)

    if("Most Variable 100" %in% getGeneSets()){
      Vs <- apply(X=x$ex[,s],MARGIN=1,FUN = function(x){var(x)})
      Fs <-  which(Vs>quantile(Vs,1-100/length(Vs)))
      f <- c(f,Fs); g$Most.Variable.100 <- "";  g$Most.Variable.100[Fs] <- "variable" ; g$Most.Variable.100 <- as.factor(g$Most.Variable.100)}
    if("Most Variable 1000" %in% getGeneSets()){
      Vs <- apply(X=x$ex[,s],MARGIN=1,FUN = function(x){var(x)})
      Fs <-  which(Vs>quantile(Vs,1-1000/length(Vs)))
      f <- c(f,Fs); g$Most.Variable.1000 <- "";  g$Most.Variable.1000[Fs] <- "Most" ; g$Most.Variable.1000 <- as.factor(g$Most.Variable.1000)}

    for(this_gene_list in intersect(getGeneSets(),c(names(globals$gene_lists), names(globals$mouse_gene_lists)))){
      gene_list_obj <- globals$gene_lists[[which(names(globals$gene_lists)==this_gene_list)]]
      if(is.data.frame(gene_list_obj)){
        # fancy gene list with attributes
        symbol_overlap <- intersect(x$featInfo$SYMBOL,gene_list_obj[,1])
        Fs <- match(table = x$featInfo$SYMBOL , x = symbol_overlap)
        Freverse <- match(table = gene_list_obj[,1] , symbol_overlap)
        f <- c(f,Fs)
        g$new <- ""
        g[Fs,length(g)] <- as.character(gene_list_obj[Freverse,2] )
        g$new <- as.factor(g$new)
        names(g)[length(g)] <- this_gene_list
      } else{
        # just a list of gene symbols
        # print(this_gene_list)
        # print(gene_list_obj)
        Fs <- which(x$featInfo$SYMBOL %in% gene_list_obj)
        #print(Fs)
        f <- c(f,Fs)
        #print(f)
        #print("------------")
        g$new <- ""
        g[Fs,length(g)] <- this_gene_list
        g$new <- as.factor(g$new)
        names(g)[length(g)] <- this_gene_list
        #print(head(g, n = 10))
        #levels(g$new)
      }

    }
    if("User selected genes 1" %in% getGeneSets()){
      userGenes <- strsplit(input$userGeneList.1,"[,;[:space:]]+")[[1]]
      Fs <- which(x$featInfo$SYMBOL %in% userGenes)
      f <- c(f,Fs); g$User.selected.genes.1 <- ""; g$User.selected.genes.1[Fs] <- "user.1" ; g$User.selected.genes.1 <- as.factor(g$User.selected.genes.1) }
    if("User selected genes 2" %in% getGeneSets()){
      userGenes <- strsplit(input$userGeneList.2,"[,;[:space:]]+")[[1]]
      Fs <- which(x$featInfo$SYMBOL %in% userGenes)
      f <- c(f,Fs); g$User.selected.genes.2 <- ""; g$User.selected.genes.2[Fs] <- "user.2" ; g$User.selected.genes.2 <- as.factor(g$User.selected.genes.2) }
    if("User selected genes 3" %in% getGeneSets()){
      userGenes <- strsplit(input$userGeneList.3,"[,;[:space:]]+")[[1]]
      Fs <- which(x$featInfo$SYMBOL %in% userGenes)
      f <- c(f,Fs); g$User.selected.genes.3 <- ""; g$User.selected.genes.3[Fs] <- "user.3" ; g$User.selected.genes.3 <- as.factor(g$User.selected.genes.3) }
    f <- f[!duplicated(f)]
    validate(
      need(length(f) > 0, "No matches found for selected genes")
    )
    print("returning now!")
    return(list(f=f,g=g[f,!(names(g)=="SYMBOL"),drop=FALSE]))
  })

  mergeDataSets <- function(setA,setB){
    if(is.null(setA)){return(setB)}
    if(is.null(setB)){return(setA)}
    commonGenes <- intersect(setA$featInfo$SYMBOL,setB$featInfo$SYMBOL)
    commonGenes <- setdiff(commonGenes,NA)
    subInA <- which(setA$featInfo$SYMBOL %in% commonGenes)
    setA$featInfo <- setA$featInfo[subInA,,drop=FALSE]
    setA$ex <- setA$ex[subInA,,drop=FALSE]
    subInB <- which(setB$featInfo$SYMBOL %in% commonGenes)
    setB$featInfo <- setB$featInfo[subInB,,drop=FALSE]
    setB$ex <- setB$ex[subInB,,drop=FALSE]
    subInA <- order(setA$featInfo$SYMBOL)
    setA$featInfo <- setA$featInfo[subInA,,drop=FALSE]
    setA$ex <- setA$ex[subInA,,drop=FALSE]
    if(anyDuplicated(setA$featInfo$SYMBOL)){
      subInA <- which(!duplicated(setA$featInfo$SYMBOL))
      setA$ex <- aggregate(setA$ex, by = list(setA$featInfo$SYMBOL), FUN = mean)
      setA$ex <- setA$ex[,-1,drop=FALSE]
      setA$featInfo <- setA$featInfo[subInA,,drop=FALSE]
      rownames(setA$ex) <- as.character(setA$featInfo$SYMBOL)
    }
    subInB <- order(setB$featInfo$SYMBOL)
    setB$featInfo <- setB$featInfo[subInB,,drop=FALSE]
    setB$ex <- setB$ex[subInB,,drop=FALSE]
    if(anyDuplicated(setB$featInfo$SYMBOL)){
      subInB <- which(!duplicated(setB$featInfo$SYMBOL))
      setB$ex <- aggregate(setB$ex, by = list(setB$featInfo$SYMBOL), FUN = mean)
      setB$ex <- setB$ex[,-1,drop=FALSE]
      setB$featInfo <- setB$featInfo[subInB,,drop=FALSE]
      rownames(setB$ex) <- as.character(setB$featInfo$SYMBOL)
    }
    combined <- list()
    combined$ex <- cbind(setA$ex,setB$ex)
    if(ncol(setA$featInfo) == ncol(setB$featInfo)){
      if(nrow(setA$featInfo) == nrow(setB$featInfo)){
        combined$featInfo <- setA$featInfo
      }
    }
    else{combined$featInfo <- cbind(setA$featInfo,setB$featInfo)
    combined$featInfo <- combined$featInfo[,!duplicated(names(combined$featInfo))]}
    combined$sampInfo   <- plyr::rbind.fill(setA$sampInfo,setB$sampInfo)
    combined$sampInfo <- data.frame(
      lapply(X = combined$sampInfo,FUN = function(x){
        if(is.factor(x)){
          levels(x) <- c(levels(x),"NotAvailable")
          x[is.na(x)] <- "NotAvailable"
        }
        if(is.numeric(x)){
          x[is.na(x)] <- NaN
        }
        if(is.character(x)){
          x[is.na(x)] <- "NotAvailable"
        }
        return(x)})
    )
    # print(class(combined$featInfo))
    # print(names(combined$featInfo))
    # print(dim(combined$featInfo))
    # print("Done merging.")
    return(combined)
  }

  # ==========================================================
  # select only a subset of samples
  sampleSet <- reactive({
    print(". . . Running sampleSet . . .")
    sampInfo <- dataSet()$sampInfo
    sampleCategories <- input$sampleFilters
    sampleCategories <- gsub(x = sampleCategories,
                             pattern = " ",
                             replacement = ".",
                             fixed = T)
    print("========== sampleCategories =========")
    print(sampleCategories)

    selectedFilters <- input$sampleFilterSpecifics
    print("========== selectedFilters =========")
    print(selectedFilters)

    selectedSamples <- 1:(nrow(sampInfo))
    if(length(selectedFilters)<1){
      return(selectedSamples) }

    for(i in sampleCategories){
      withinCategoryData <- sampInfo[,names(sampInfo)==i]
      print(withinCategoryData)
      withinCategoryToRemove <- stringr::str_split_fixed(string = selectedFilters, pattern = ":", n = 2)[,2]
      selectedSamples <- intersect(selectedSamples,
                                   which(!(withinCategoryData %in% withinCategoryToRemove)))
    }
    return(selectedSamples)
  })

  # ==========================================================
  # pull samples and features from loaded data
  getX <- reactive({
    print(". . . Running getX . . .")
    x <- dataSet()
    s <- sampleSet()
    f <- featureSet()

    tmp.ex <- x$ex[f$f,s,drop=FALSE]
    x$Vs <- apply(X=tmp.ex,MARGIN=1,FUN = function(x){var(x)})
    x$Vs[is.na(x$Vs)] <- 0

    if("hide" %in% input$remove.zeros){
      f$f <- subset(f$f,!(x$Vs == 0))
      f$g <- subset(f$g,!(x$Vs == 0))
    }


    x$featInfo <-  x$featInfo[f$f,,drop=FALSE]
    x$ex <- x$ex[f$f,s,drop=FALSE]
    x$sampInfo <- x$sampInfo[s,,drop=FALSE]
    return(x)
  })

  getSignatures <- reactive({
    print(". . . Running getSignatures . . .")
    x <- dataSet()
    s <- sampleSet()
    x$ex <- x$ex[,s,drop=FALSE]
    tmp.info <- NULL

    print(paste0("SignatureTrack.1 = ",input$SignatureTrack.1))
    print(paste0("SignatureTrack.2 = ",input$SignatureTrack.2))
    print(paste0("SignatureTrack.3 = ",input$SignatureTrack.3))
    Expression.signature <- calculate.signature(dataset = x , track = input$SignatureTrack.1)
    tmp.info$Expression.signature.1 <- Expression.signature

    Expression.signature <- calculate.signature(dataset = x , track = input$SignatureTrack.2)
    tmp.info$Expression.signature.2 <- Expression.signature

    Expression.signature <- calculate.signature(dataset = x , track = input$SignatureTrack.3)
    tmp.info$Expression.signature.3 <- Expression.signature

    rownames(x$ex) <- make.names(x$featInfo$SYMBOL, unique=TRUE)

    if(!"Mouse" %in% input$Species){
      for(i in 1:length(globals$classifier_list[[1]])){
        this.classifier <-
          get(
            as.character(
              globals$classifier_list$variablenames[i]))
        tmp.info[[as.character(globals$classifier_list$labels[i])]] <-
          as.numeric(
            create.classif(
              dat=x$ex,
              fit=this.classifier$fit,
              classifier=this.classifier)$predprob)
      }
      ## Add PAMG to classifier obj
      tmp.info = c(tmp.info, implement_PAMG(x))

      ## Bin purIST outputs
      tmp.info$purIST_2019_Call = factor(ifelse(tmp.info$puRIST_2019 >= .5,"Basal-like","Classical"))

      return(tmp.info)
    }
  })

  calculate.signature <- function(dataset , track){
    print(paste0("Track = ",track))
    if("User selected genes 1" %in% track){
      gene_list_obj <- strsplit(input$userGeneList.1,"[,;[:space:]]+")[[1]]
    }else if("User selected genes 2" %in% track){
      gene_list_obj <- strsplit(input$userGeneList.2,"[,;[:space:]]+")[[1]]
    }else if("User selected genes 3" %in% track){
      gene_list_obj <- strsplit(input$userGeneList.3,"[,;[:space:]]+")[[1]]
    }else{
      gene_list_obj <- globals$gene_lists[[which(names(globals$gene_lists)==track)]]
    }

    if(is.data.frame(gene_list_obj)){gene_list_obj <- gene_list_obj[[1]]}
    symbol_overlap <- intersect(dataset$featInfo$SYMBOL,as.character(gene_list_obj))
    Fs <- match(table = dataset$featInfo$SYMBOL ,
                x     = symbol_overlap)
    y <- as.matrix(dataset$ex[Fs,,drop=FALSE])
    Expression.signature <- colMeans(x = y,
                                     na.rm = TRUE)
    return(Expression.signature)
  }

  # ==========================================================
  # do column clustering
  getColv <- reactive({
    print(". . . Running column clustering . . .")
    x <- getX()
    # note that x already includes log scaled expression at this point.
    if(length(x$ex[,1]) > 1){
      if("Consensus.RowScaled" %in% input$sampleClustertype){
        ex2 <- t(scale(t(x$ex),center=TRUE,scale=TRUE))
        k <- input$sampleSortBy
        Colv <- as.dendrogram(
          ConsensusClusterPlus::ConsensusClusterPlus(d = (as.matrix(ex2)),
                                                     seed = 1234,
                                                     maxK = k+1,
                                                     reps=50,
                                                     distance="pearson",
                                                     clusterAlg="pam")[[k]]$consensusTree)
      } else if ("Consensus" %in% input$sampleClustertype){
        k <- input$sampleSortBy
        Colv <- as.dendrogram(
          ConsensusClusterPlus::ConsensusClusterPlus(d = (as.matrix(x$ex)),
                                                     seed = 1234,
                                                     maxK = k+1,
                                                     reps = 50,
                                                     distance = "pearson",
                                                     clusterAlg = "pam")[[k]]$consensusTree)
      } else if("Pearson" %in% input$sampleClustertype){
        Colv <- as.dendrogram( hclust( bioDist::cor.dist(x = t(as.matrix(x$ex)) ) ) )
      } else if("Euclidean" %in% input$sampleClustertype){
        Colv <- as.dendrogram( hclust( dist(x = t(as.matrix(x$ex)),
                                            method="euclidean")))
      } else if("K.Means" %in% input$sampleClustertype){
        k <- input$sampleSortBy
        set.seed(1234)
        ProcessedKmeans <-(kmeans((x = t(as.matrix(x$ex))),
                                  centers = k,
                                  iter.max = 20,
                                  nstart = 1))
        Colv <- convert_kmeans_to_dendrogram(ProcessedKmeans$cluster)
      }
    } else { Colv <- FALSE }
    if("Sorted" %in% input$sampleClustertype){
      if(input$sampleSortBy == "expression"){
        sampleOrder <- order(colMeans(x$ex))
      } else {
        tmp.info <- getSignatures()
        if(!"Mouse" %in% input$Species){       x$sampInfo <- cbind(x$sampInfo,tmp.info)     }
        sampleOrder <- order(x$sampInfo[,input$sampleSortBy])
      }
      Colv <- convert_order_to_dendrogram(sampleOrder)
    }
    return(Colv)
  })

  # ==========================================================
  # do row clustering
  getRowv <- reactive({
    print(". . . Running row clustering . . .")
    x <- getX()
    f <- featureSet()

    if("hide" %in% input$remove.zeros){
      f$f <- subset(f$f,!(x$Vs == 0))
      f$g <- subset(f$g,!(x$Vs == 0))
    }

    friendly.filter <- character()
    for(i in input$genesets.2){
      friendly.filter <-  c(friendly.filter, gsub(x = i,
                                                  pattern = " ",
                                                  replacement = ".",
                                                  fixed = TRUE))
    }

    featureSetFilter <- which(names(f$g) %in% friendly.filter)

    print(featureSetFilter)

    genesToPlot <- numeric()
    for(i in featureSetFilter){
      genesToPlot <- union(genesToPlot, which(!(f$g[,i] %in% "")))
    }


    x$featInfo <- x$featInfo[genesToPlot,,drop = FALSE]
    x$ex <- x$ex[genesToPlot,,drop = FALSE]
    x$Vs <- x$Vs[genesToPlot]
    f$f <- f$f[genesToPlot]
    f$g <- f$g[genesToPlot,,drop = FALSE]

    ex2 <- t(scale(t(x$ex),center=TRUE,scale=TRUE))
    if(length(x$ex[,1])<2){
      Rowv <- FALSE
    } else if("Consensus" %in% input$geneClustertype){
      k <- input$geneSortBy
      Rowv <- as.dendrogram(
        ConsensusClusterPlus::ConsensusClusterPlus(d = t(as.matrix(ex2)),
                                                   seed = 1234,
                                                   maxK = k+1,
                                                   reps=50,
                                                   distance="pearson",
                                                   clusterAlg="pam")[[k]]$consensusTree)

    } else if("Pearson" %in% input$geneClustertype){
      Rowv <- as.dendrogram( hclust( bioDist::cor.dist(x = (as.matrix(ex2)) ) ) )


    } else if("Euclidean" %in% input$geneClustertype){
      Rowv <- as.dendrogram( hclust( dist(x = (as.matrix(ex2)),
                                          method="euclidean")))

    } else if("K.Means" %in% input$geneClustertype){
      k <- input$geneSortBy
      set.seed(1234)
      ProcessedKmeans <-(kmeans((x = (as.matrix(ex2))),
                                centers = k,
                                iter.max = 25,
                                nstart = 1))
      Rowv <- convert_kmeans_to_dendrogram(ProcessedKmeans$cluster)


    }    else if("Sorted" %in% input$geneClustertype){
      if(input$geneSortBy == "expression"){
        geneOrder <- order(rowMeans(x$ex))
      } else if(input$geneSortBy == "gene sets"){
        geneOrder <- do.call(order, as.list(f$g))
      }
      Rowv <- convert_order_to_dendrogram(geneOrder)
    }

    return(Rowv)
  })

  # ==========================================================
  # do t-SNE projection
  getTSNE <- reactive({
    x <- getX()
    samples.to.remove <- which(duplicated(t(x$ex)))
    if(length(samples.to.remove)>0){
      x$ex <- x$ex[,-samples.to.remove]
    }

    if (input$projection %in% "tSNE") {
      print(". . . Running tSNE . . .")
      set.seed(1234)
      tsne <- NA
      tsne <- Rtsne(X = t(x$ex),
                    dims = 2,
                    initial_dims = 20, # 50
                    perplexity = 10, # 30
                    theta = 0.5, # 0.5
                    max_iter = 1000, # 1000
                    verbose = FALSE,
                    momentum = 0.5, # 0.5
                    final_momentum = 0.8,  # 0.8
                    eta = 200, # 200
                    exaggeration_factor = 12 # 12
      )
    }
    # if (input$projection %in% "UMAP") {
    #   print(". . . Running UMAP . . .")
    #
    #   tsne <- NA
    #   tsne <- umap(t(x$ex),init="spectral",metric="cosine",random_state = 1234L)
    #   tsne<-list(Y=cbind(tsne$UMAP1,tsne$UMAP2))
    #
    # }

    else {
      pcaData  <- prcomp(x$ex)
      tsne <- list(Y = cbind(pcaData$rotation[,1], pcaData$rotation[,2]))
    }
    tsne$samples.to.remove <- samples.to.remove
    return(tsne)
  })


  # ==========================================================
  #
  getDendrogramParam <- reactive({
    Colv <- getColv()
    Rowv <- getRowv()
    if(((Rowv %in% FALSE) && (Colv %in% FALSE))[1]){
      dendrogramParam <- "none"
    } else if((Rowv %in% FALSE)[1]) {
      dendrogramParam <- "column"
    } else if((Colv %in% FALSE)[1]) {
      dendrogramParam <- "row"
    } else {
      dendrogramParam <- "both" }
    return(dendrogramParam)
  })
  # ==========================================================
  # Get survival data from ReactiveSurvival

  getsurvivalParam <- reactive({
    tmp <- names(dataSet()$metadata)[which(dataSet()$metadata %in% input$survivalDays)]
    return(tmp)
  })

  # ================================================================
  # ------------- Reactive outputs ---------------------------------
  # ================================================================

  # ================================================================
  # write heatmap to table
  observe({
    if(r$heatmapStoplight == "green"){
      writeHeatmapToTable <- function() {
        #  ------- retrieve data set and clustering results  -------
        x <- getX()

        tmp.info <- getSignatures()
        if(!"Mouse" %in% input$Species){       x$sampInfo <- cbind(x$sampInfo,tmp.info)     }

        f <- featureSet()

        if("hide" %in% input$remove.zeros){
          f$f <- subset(f$f,!(x$Vs == 0))
          f$g <- subset(f$g,!(x$Vs == 0))
        }

        friendly.filter <- character()
        for(i in input$genesets.2){
          friendly.filter <-  c(friendly.filter, gsub(x = i,
                                                      pattern = " ",
                                                      replacement = ".",
                                                      fixed = TRUE))
        }

        featureSetFilter <- which(names(f$g) %in% friendly.filter)

        print(featureSetFilter)

        genesToPlot <- numeric()
        for(i in featureSetFilter){
          genesToPlot <- union(genesToPlot, which(!(f$g[,i] %in% "")))
          #genesToPlot <- c(genesToPlot, which(f$g[,i] %in% friendly.filter))
        }


        x$featInfo <- x$featInfo[genesToPlot,, drop = FALSE]
        x$ex <- x$ex[genesToPlot,,drop = FALSE]
        x$Vs <- x$Vs[genesToPlot]
        f$f <- f$f[genesToPlot]
        f$g <- f$g[genesToPlot,,drop = FALSE]

        Colv <- getColv()
        Rowv <- getRowv()
        sample.order <- order.dendrogram(Colv)
        feat.order <- rev(order.dendrogram(Rowv))
        # ------------- prepare top and left margin data ---------------
        featInfoToPlot <- as.matrix(f$g[feat.order,,drop=FALSE])
        featInfoNames <- t(as.matrix(names(f$g)))
        sampInfoToPlot <- t(as.matrix(x$sampInfo[sample.order,
                                                 names(x$sampInfo) %in% getSampleTracks(),
                                                 drop=FALSE]))
        sampInfoNames <- as.matrix(names(x$sampInfo)[names(x$sampInfo) %in% getSampleTracks(),
                                                     drop=FALSE])
        # ------------- construct output table -------------------------
        table.data <- cbind(as.matrix(x$featInfo$SYMBOL[feat.order]),as.matrix(x$ex)[feat.order,sample.order])
        table.data <- rbind(cbind(sampInfoNames,sampInfoToPlot),
                            table.data)
        if(!(input$sampleClustertype %in% "Sorted")){
          for(k in 2:min(5,length(sample.order))){
            cluster.labels = t(as.matrix(cutree(as.hclust(Colv), k = k))[sample.order])
            table.data <- rbind(cbind(paste("k",k,sep="."),cluster.labels),table.data)
          }
        }
        row.padding <- (dim(table.data)[1]-dim(featInfoToPlot)[1]) - 1
        if(row.padding<1){
          table.data <- rbind("",table.data)
          row.padding <- (dim(table.data)[1]-dim(featInfoToPlot)[1]) - 1
        }
        col.padding <- dim(featInfoToPlot)[2]
        table.data <- cbind(rbind(matrix("",nrow = row.padding,ncol = col.padding),
                                  featInfoNames,
                                  featInfoToPlot),
                            table.data)
        # ------------- write output table -------------------------
        write.table(x = table.data,
                    row.names = FALSE,
                    col.names = FALSE,
                    sep = "\t",
                    file = "mytable.tsv")
      }
    }
  })


  # ================================================================
  # plot heatmap to pdf file
  observe({
    if(r$heatmapStoplight == "green"){
      writeHeatmapToPdf <- function() {
        #  ------- retrieve data set and clustering results  -------
        x <- getX()

        tmp.info <- getSignatures()
        if(!"Mouse" %in% input$Species){       x$sampInfo <- cbind(x$sampInfo,tmp.info)     }

        f <- featureSet()

        if("hide" %in% input$remove.zeros){
          f$f <- subset(f$f,!(x$Vs == 0))
          f$g <- subset(f$g,!(x$Vs == 0))
        }

        friendly.filter <- character()
        for(i in input$genesets.2){
          friendly.filter <-  c(friendly.filter, gsub(x = i,
                                                      pattern = " ",
                                                      replacement = ".",
                                                      fixed = TRUE))
        }

        featureSetFilter <- which(names(f$g) %in% friendly.filter)

        print(featureSetFilter)

        genesToPlot <- numeric()
        for(i in featureSetFilter){
          genesToPlot <- union(genesToPlot, which(!(f$g[,i] %in% "")))
          #genesToPlot <- c(genesToPlot, which(f$g[,i] %in% friendly.filter))
        }

        x$featInfo <- x$featInfo[genesToPlot,, drop = FALSE]
        x$ex <- x$ex[genesToPlot,,drop = FALSE]
        x$Vs <- x$Vs[genesToPlot]
        f$f <- f$f[genesToPlot]
        f$g <- f$g[genesToPlot,,drop = FALSE]


        Colv <- getColv()
        Rowv <- getRowv()
        dendrogramParam <- getDendrogramParam()

        # ------------- visualization options ----------------------
        if("Row" %in% input$scaling){
          scale <- "row"
          breaks <- (((0:299)/299)-0.5)*5
        } else if("None" %in% input$scaling){
          scale <- "none"
          breaks <- (((0:299)/299))*quantile(c(as.matrix(x$ex)),0.99,na.rm=TRUE)
        }
        myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

        # ------------- prepare for plotting ----------------------
        RowSideColors <- getSideColors(f$g,names(f$g),drop.levels = TRUE)
        colorlists <- rep(list(c("gray94", "blue", "green",
                                 "yellow", "orange", "red","black")),
                          length(getSampleTracks()))
        colorlists[getSampleTracks() %in% classifier_list$labels] <-
          classifier_list$colors[classifier_list$labels %in% getSampleTracks()]
        ColSideColors <- getSideColors(sampInfo = x$sampInfo,
                                       sampleTracks = getSampleTracks(),
                                       colorlists = colorlists,
                                       drop.levels = TRUE)
        if(length(x$ex[,1])==1){
          x$ex <- rbind(x$ex,x$ex)
          RowSideColors$SideColors <- rbind(RowSideColors$SideColors,RowSideColors$SideColors)
        }
        # ------------- plot pdf ---------------------------------------
        scaleFactor <- sqrt(1+length(x$featInfo$SYMBOL))
        pdf("myplot.pdf", width=12, height=2.5+scaleFactor)
        suppressWarnings(
          heatmap.3(x = as.matrix(x$ex)
                    ,main = paste(input$dataset,collapse=" and ")
                    ,Rowv = Rowv
                    ,Colv = Colv
                    ,dendrogram = dendrogramParam
                    ,trace="none"
                    ,scale=scale
                    ,labRow = x$featInfo$SYMBOL
                    ,col=myPalette,breaks=breaks
                    ,ColSideColors = ColSideColors$SideColors
                    ,ColSideColorsSize = dim(ColSideColors$SideColors)[2]+3/scaleFactor
                    ,RowSideColors = t(RowSideColors$SideColors)
                    ,RowSideColorsSize = dim(t(RowSideColors$SideColors))[1]*1.5
                    ,lwid = c(1,5),lhei = c(.5+3/scaleFactor,5)
                    ,margins =c(10,20)
                    ,labCol = x$sampInfo[,names(x$sampInfo) %in% input$XaxisLabel]
          )
        )
        legend(xy.coords(x=.86,y=.92),
               legend=c("Gene Labels","",RowSideColors$text,"Sample Labels","",ColSideColors$text),
               fill=c("white","white",RowSideColors$colors,"white","white",ColSideColors$colors),
               border=FALSE, bty="n",
               y.intersp = 0.9, cex=0.5)
        dev.off()
      }
    }
  })

  # ================================================================
  # Plot heatmaps to screen
  observe({
    if(r$heatmapStoplight == "green"){
      output$ReactiveGenesets <- renderUI({
        checkboxGroupInput(inputId = "genesets.2",
                           label = "Select genesets to retain",
                           choices = c(input$genesets),
                           selected = c(input$genesets))
      })
      output$heatmap <- renderPlot({
        #  ------- retrieve data set and clustering results  -------
        x <- getX()

        tmp.info <- getSignatures()
        if(!"Mouse" %in% input$Species){       x$sampInfo <- cbind(x$sampInfo,tmp.info)     }

        f <- featureSet()

        if("hide" %in% input$remove.zeros){
          f$f <- subset(f$f,!(x$Vs == 0))
          f$g <- subset(f$g,!(x$Vs == 0))
        }

        friendly.filter <- character()
        for(i in input$genesets.2){
          friendly.filter <-  c(friendly.filter, gsub(x = i,
                                                      pattern = " ",
                                                      replacement = ".",
                                                      fixed = TRUE))
        }


        print("The genesets.2")
        print(input$genesets.2)

        print("The friendly.filter")
        print(friendly.filter)

        featureSetFilter <- which(names(f$g) %in% friendly.filter)

        print("The featureSet Filter")
        print(featureSetFilter)

        genesToPlot <- numeric()
        for(i in featureSetFilter){
          genesToPlot <- union(genesToPlot, which(!(f$g[,i] %in% "")))
          #genesToPlot <- c(genesToPlot, which(f$g[,i] %in% friendly.filter))
        }

        x$featInfo <- x$featInfo[genesToPlot,, drop = FALSE]
        x$ex <- x$ex[genesToPlot,,drop = FALSE]
        x$Vs <- x$Vs[genesToPlot]
        f$f <- f$f[genesToPlot]
        f$g <- f$g[genesToPlot,,drop = FALSE]

        Colv <- getColv()
        Rowv <- getRowv()

        dendrogramParam <- getDendrogramParam()
        # ------------- visualization options ----------------------
        if("Row" %in% input$scaling){
          scale <- "row"
          breaks <- (((0:299)/299)-0.5)*5
        } else if("None" %in% input$scaling){
          scale <- "none"
          breaks <- (((0:299)/299))*quantile(c(as.matrix(x$ex)),0.99,na.rm=TRUE)
        }
        myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

        # ------------- prepare for plotting ----------------------
        ## Colorlist is always rainbow regardless. Want to make continuous variables blue-white-red
        RowSideColors <- getSideColors(f$g,names(f$g),drop.levels = FALSE)

        # colorlists <- rep(list(c("gray94", "blue", "green",
        #                          "yellow", "orange", "red","black")),
        #                   length(getSampleTracks()))
        # print("THESE ARE THE SAMPLETRACKS()")
        # print(getSampleTracks())
        # print(" THOSE were the tracks ")

        colorlists = list()
        for(colname in getSampleTracks()){
          if(class(x$sampInfo[,colname]) == "numeric"){
            colorlists = c(colorlists,
                           list(c("blue","white","red")))
          }
          else{
            colorlists = c(colorlists,
                           list(c("gray94", "blue", "green",
                                  "yellow", "orange", "red","black")))
          }
        }
        colorlists[getSampleTracks() %in% classifier_list$labels] <-
          classifier_list$colors[classifier_list$labels %in% getSampleTracks()]
        ColSideColors <- getSideColors(sampInfo = x$sampInfo,
                                       sampleTracks = getSampleTracks(),
                                       colorlists = colorlists,
                                       drop.levels = TRUE)
        if(length(x$ex[,1])==1){
          x$ex <- rbind(x$ex,x$ex)
          RowSideColors$SideColors <- rbind(RowSideColors$SideColors,RowSideColors$SideColors)
        }

        # ------------- plot to screen ---------------------------------------
        print(x$ex)
        suppressWarnings(
          heatmap.3(x = as.matrix(x$ex)
                    ,Rowv = Rowv
                    ,Colv = Colv
                    ,dendrogram = dendrogramParam
                    ,trace="none"
                    ,scale=scale
                    ,labRow = x$featInfo$SYMBOL
                    ,col=myPalette,breaks=breaks
                    ,ColSideColors = ColSideColors$SideColors
                    ,ColSideColorsSize = dim(ColSideColors$SideColors)[2]
                    ,RowSideColors = t(RowSideColors$SideColors)
                    ,RowSideColorsSize = dim(t(RowSideColors$SideColors))[1]*1.5
                    ,lwid = c(1,5),lhei = c(1,5)
                    ,margins =c(10,30)
                    ,labCol = x$sampInfo[,names(x$sampInfo) %in% input$XaxisLabel]
          )
        )

        legend.scale <- 0.90
        if(length(RowSideColors$text) + length(ColSideColors$text)<18){legend.scale <- 1.2}
        if(length(RowSideColors$text) + length(ColSideColors$text)>32){legend.scale <- 0.67}
        legend(xy.coords(x=.92,y=.92),
               legend=c("Gene Labels","",RowSideColors$text,"Sample Labels","",ColSideColors$text),
               fill=c("white","white",RowSideColors$colors,"white","white",ColSideColors$colors),
               border=FALSE, bty="n",
               y.intersp = legend.scale , cex = legend.scale)
      })
    }
  })

  # ==========================================================
  output$cartesian <- renderPlot({
    print("Rendering cartesian plot")
    x <- getX()

    tmp.info <- getSignatures()
    if(!"Mouse" %in% input$Species){       x$sampInfo <- cbind(x$sampInfo,tmp.info)     }

    tsne <- getTSNE()

    if(length(tsne$samples.to.remove)>0){
      x$ex <- x$ex[,-tsne$samples.to.remove]
      x$sampInfo <- x$sampInfo[-tsne$samples.to.remove,]
    }

    df <- data.frame(x = tsne$Y[,1],
                     y = tsne$Y[,2])

    if(input$painting %in% "Signatures (R+G+B)"){
      df$R <- as.hexmode(sapply(round(400*x$sampInfo$Expression.signature.1/
                                        max(x$sampInfo$Expression.signature.1)),
                                FUN = function(x){min(255,x)}))

      df$G <- as.hexmode(sapply(round(400*x$sampInfo$Expression.signature.2/
                                        max(x$sampInfo$Expression.signature.2)),
                                FUN = function(x){min(255,x)}))

      df$B <- as.hexmode(sapply(round(400*x$sampInfo$Expression.signature.3/
                                        max(x$sampInfo$Expression.signature.3)),
                                FUN = function(x){min(255,x)}))

      df$RGB <- paste("#",df$R,df$G,df$B,sep="")
    } else if(input$painting %in% "Categorical"){
      SideColors <- getSideColors(sampInfo = x$sampInfo,
                                  sampleTracks = input$XaxisLabel,
                                  colorlists = c("blue", "green", "yellow", "orange", "red"),
                                  drop.levels = TRUE)
      df$RGB <- SideColors$SideColors
    }


    if(input$painting %in% "Signatures (R+G+B)"){
      df$R <- as.hexmode(sapply(round(400*x$sampInfo$Expression.signature.1/
                                        max(x$sampInfo$Expression.signature.1)),
                                FUN = function(x){min(255,x)}))

      df$G <- as.hexmode(sapply(round(400*x$sampInfo$Expression.signature.2/
                                        max(x$sampInfo$Expression.signature.2)),
                                FUN = function(x){min(255,x)}))

      df$B <- as.hexmode(sapply(round(400*x$sampInfo$Expression.signature.3/
                                        max(x$sampInfo$Expression.signature.3)),
                                FUN = function(x){min(255,x)}))

      df$RGB <- paste("#",df$R,df$G,df$B,sep="")
    } else if(input$painting %in% "Categorical"){
      SideColors <- getSideColors(sampInfo = x$sampInfo,
                                  sampleTracks = input$XaxisLabel,
                                  colorlists = c("blue", "green", "yellow", "orange", "red"),
                                  drop.levels = TRUE)
      df$RGB <- SideColors$SideColors
    }


    #=================================================================
    p <- ggplot(data = df,
                aes(x = x,
                    y = y))
    p <- p + geom_point(fill = df$RGB,
                        color = "gray40",
                        size = 4,
                        alpha = 0.7,
                        pch =21)
    p <- p + labs(title = paste(input$projection))
    #p <- p + theme_dark()
    p <- p + theme_classic()
    # p <- p + theme(panel.background = element_rect(fill = "black"),
    #                panel.grid.major = element_line(colour= "gray25"),
    #                panel.grid.minor = element_line(colour= "gray30"))
    if(input$painting %in% "Categorical"){
      SideColors$text <- SideColors$text[c(-1,-2,-length(SideColors$text))]
      SideColors$colors <- SideColors$colors[c(-1,-2,-length(SideColors$colors))]
      p <- p + geom_point(data = data.frame(x=0,y=0,
                                            t=SideColors$text,
                                            c=SideColors$colors),
                          aes(x=x,y=y,color=t),size = 5, alpha = 0)
      p <- p + scale_color_manual(labels = SideColors$text,
                                  values = SideColors$colors)
      p <- p + guides(colour = guide_legend(override.aes = list(alpha=1)))
    }
    plot(p)


  })
  # ==========================================================
  output$barplot <- output$barplot2 <- renderPlot({
    validate(
      need(length(getSampleTracks()) > 0, "Please select a sample track")
    )
    x <- getX()
    tmp.info <- getSignatures()

    if(!"Mouse" %in% input$Species){       x$sampInfo <- cbind(x$sampInfo,tmp.info)     }
    f <- featureSet()

    if("hide" %in% input$remove.zeros){
      f$f <- subset(f$f,!(x$Vs == 0))
      f$g <- subset(f$g,!(x$Vs == 0))
    }

    #print(summary(f$g))

    tmp.ex <- data.frame(matrix(nrow = length(x$sampInfo[[1]]), ncol=0))
    for(gs.name in names(f$g)){
      tmp.ex[[gs.name]] <- colMeans(x$ex[!(f$g[[gs.name]] %in% ""),,drop = FALSE])
    }
    print(". . . tmp.ex . . . ")
    print(summary(tmp.ex))
    color_by <- getColor()

    if(color_by == "No.Color"){
      tmp.si <- as.data.frame(x$sampInfo[,getSampleTracks()])
      names(tmp.si) <- as.character(getSampleTracks())
    }
    else{
      tmp.si <- as.data.frame(x$sampInfo[,c(getSampleTracks(),color_by)])
      #names(tmp.si) <- as.character(c(getSampleTracks(), color_by))
    }


    print(summary(tmp.si))
    tmp.df <- cbind(tmp.ex,tmp.si)

    print(summary(tmp.df))

    melted.df <- melt(data = tmp.df,
                      id.vars = names(tmp.si) ,
                      variable.names = names(tmp.ex))

    names(melted.df)[names(melted.df) %in% "variable"] <- "SYMBOL"
    names(melted.df)[names(melted.df) %in% "value"] <- "expression"
    names(melted.df)[names(melted.df) %in% color_by] <- "color"

    if(color_by == "No.Color"){
      # They dont want to color

      melted.df <- melt(data = melted.df,
                        id.vars = c("SYMBOL","expression") ,
                        variable.names = names(tmp.si))
      #melted.df <- melted.df[complete.cases(melted.df),]
      # print("df for dotplot")
      # print(head(melted.df))
      # print(class(melted.df$value))
      ggplot(data =  melted.df) +
        facet_grid(SYMBOL ~ variable, scales="free") +
        geom_dotplot(data =  subset(melted.df,
                                    (is.na(as.numeric(as.character(melted.df$value)))) &
                                      !is.na(melted.df$value)),
                     aes(x=value,y=expression),
                     binaxis='y',
                     stackdir='center',
                     dotsize=.3) +
        geom_point(data = subset(melted.df,!(is.na(as.numeric(as.character(melted.df$value))))),
                   aes(x=as.numeric(as.character(value)),y=expression),
                   size = 2#,clip = "off"
        ) +
        stat_cor(data = subset(melted.df,!(is.na(as.numeric(as.character(melted.df$value))))),
                 aes(x=as.numeric(as.character(value)),y=expression),
                 method = "spearman",
                 label.y.npc = "top") +
        geom_smooth(data = subset(melted.df,!(is.na(as.numeric(as.character(melted.df$value))))),
                    aes(x=as.numeric(as.character(value)),y=expression),
                    method = "lm",
                    se = F,
                    color = "Black",
                    lty = 2) + theme_pubclean() +
        theme(axis.text.x = element_text(angle=20,margin = margin(t = 15)))
    }
    else {
      # They want to color by something new
      melted.df <- melt(data = melted.df,
                        id.vars = c("SYMBOL","expression", "color") ,
                        variable.names = names(tmp.si))
      #melted.df <- melted.df[complete.cases(melted.df),]

      ggplot(data =  melted.df) +
        facet_grid(SYMBOL ~ variable, scales="free") +
        geom_dotplot(data =  subset(melted.df,
                                    (is.na(as.numeric(as.character(melted.df$value)))) &
                                      !is.na(melted.df$value)),
                     aes(x=value,y=expression, fill = color, color = color),
                     binaxis='y',
                     stackdir='center',
                     dotsize=.3) +
        geom_point(data =  subset(melted.df,
                                  !(is.na(as.numeric(as.character(melted.df$value))))),
                   aes(x=as.numeric(as.character(value)),y=expression, color = color),
                   size = 2#,clip = "off"
        ) +
        theme_pubclean() + theme(axis.text.x = element_text(angle=20,margin = margin(t = 15)))
    }
  })
  # ==========================================================
  output$survival <- renderPlot({
    validate(need(!(dataSet()$metadata$survivalA %in% "None"), "This dataset does not have survival data, please pick another set"))
    validate(need((3 > length(getSampleTracks())) && (length(getSampleTracks()) > 0), "Please select 1-2 sample tracks"))
    print("..Rendering survival plot..")
    print(input$survivalType)

    x <- dataSet()
    s = sampleSet()
    x$sampInfo = x$sampInfo[s,,drop=T]
    # x = getX()
    tmp.info <- getSignatures()
    x$sampInfo <- cbind(x$sampInfo,tmp.info)
    survivalParam <- getsurvivalParam()

    if(input$survivalType == "Kaplan Meier"){
      print("..Running KM Analysis..")

      tmp.df <- na.omit(x$sampInfo[, c(getSampleTracks(),
                                       survivalParam,
                                       "censorA.0yes.1no")])
      tmp.df$censorA.0yes.1no <- as.integer(as.character(tmp.df$censorA.0yes.1no))

      if(ncol(tmp.df) == 3){
        names(tmp.df) <- c("factor1",
                           "survivalParam",
                           "censorA.0yes.1no")
        # ================================
        # This is for continuous sampInfo
        # ================================
        if(length(levels(as.factor(tmp.df$factor1))) >= 5 && is.numeric(tmp.df$factor1)){

          quants <- quantile(tmp.df$factor,
                             probs = seq(0,1,.05),
                             na.rm = T)

          tmp.df$comp <- as.factor(cut(tmp.df$factor1,
                                       breaks = c(quants[[1]], quants[names(quants) %in% percent(input$factorized1)], quants[[length(quants)]]),
                                       labels = c("Low", "High"),
                                       include.lowest = T))

          tmp.df$SurvObj <- with(tmp.df, Surv(time = survivalParam,
                                              event = censorA.0yes.1no))
          fit <- survfit(SurvObj ~ comp,
                         data = tmp.df,
                         type = "kaplan-meier")
          ggsurvplot(fit,
                     data = tmp.df,
                     pval = TRUE,
                     censor = TRUE,
                     surv.median.line = "hv",
                     legend.title = paste(getSampleTracks()),
                     xlab = "Time (days)",
                     risk.table = TRUE,
                     title = "Kaplan Meier")
        } else{
          # ================================
          # This is for categorical sampInfo
          # ================================
          tmp.df$SurvObj <- with(tmp.df, Surv(time = survivalParam,
                                              event = censorA.0yes.1no))
          fit <- survfit(SurvObj ~ factor1,
                         data = tmp.df,
                         type = "kaplan-meier")
          ggsurvplot(fit,
                     data = tmp.df,
                     pval = TRUE,
                     censor = TRUE,
                     surv.median.line = "hv",
                     legend.title = paste(getSampleTracks()),
                     xlab = "Time (days)",
                     risk.table = TRUE,
                     title = "Kaplan Meier")
        }
      } else if(ncol(tmp.df) == 4){
        names(tmp.df) <- c("factor1",
                           "factor2",
                           "survivalParam",
                           "censorA.0yes.1no")

        if(length(levels(as.factor(tmp.df$factor1))) >= 5 && is.numeric(tmp.df$factor1)){
          quants1 <- quantile(tmp.df$factor1,
                              probs = seq(0,1,.05),
                              na.rm = T)

          tmp.df$comp1 <- cut(tmp.df$factor1,
                              breaks = c(quants1[[1]], quants1[names(quants1) %in% percent(input$factorized1)], quants1[[length(quants1)]]),
                              labels = c("Low1", "High1"),
                              include.lowest = T)

          if((length(levels(as.factor(tmp.df$factor2 >= 5)))) && is.numeric(tmp.df$factor2)){
            quants2 <- quantile(tmp.df$factor2,
                                probs = seq(0,1,.05),
                                na.rm = T)

            tmp.df$comp2 <- cut(tmp.df$factor2, breaks = c(quants2[[1]], quants2[names(quants2) %in% percent(input$factorized2)], quants2[[length(quants2)]]),
                                labels = c("Low2", "High2"),
                                include.lowest = T)

            # if both factors are continuous
            tmp.df$SurvObj <- with(tmp.df, Surv(time = survivalParam,
                                                event = censorA.0yes.1no))
            fit <- survfit(SurvObj ~ comp1 + comp2,
                           data = tmp.df,
                           type = "kaplan-meier")
            ggsurvplot(fit,
                       data = tmp.df,
                       pval = TRUE,
                       censor = TRUE,
                       surv.median.line = "hv",
                       legend.title = paste(getSampleTracks()),
                       xlab = "Time (days)",
                       risk.table = TRUE,
                       title = "Kaplan Meier")
          }
          else{
            # if only factor1 is continuous
            tmp.df$SurvObj <- with(tmp.df, Surv(time = survivalParam,
                                                event = censorA.0yes.1no))
            fit <- survfit(SurvObj ~ comp1 + factor2,
                           data = tmp.df,
                           type = "kaplan-meier")
            ggsurvplot(fit,
                       data = tmp.df,
                       pval = TRUE,
                       censor = TRUE,
                       surv.median.line = "hv",
                       legend.title = paste(getSampleTracks()),
                       xlab = "Time (days)",
                       risk.table = TRUE,
                       title = "Kaplan Meier")
          }
        }
        else{
          if(length(levels(as.factor(tmp.df$factor2))) >= 5 && is.numeric(tmp.df$factor2)){
            quants2 <- quantile(tmp.df$factor2,
                                probs = seq(0,1,.05),
                                na.rm = T)

            tmp.df$comp2 <- cut(tmp.df$factor2, breaks = c(quants2[[1]], quants2[names(quants2 %in% percent(input$factorized2))], quants2[[length(quants2)]]),
                                labels = c("Low", "High"),
                                include.lowest = T)

            # if only factor2 is continuous
            tmp.df$SurvObj <- with(tmp.df, Surv(time = survivalParam,
                                                event = censorA.0yes.1no))
            fit <- survfit(SurvObj ~ factor1 + comp2,
                           data = tmp.df,
                           type = "kaplan-meier")
            ggsurvplot(fit,
                       data = tmp.df,
                       pval = TRUE,
                       censor = TRUE,
                       surv.median.line = "hv",
                       legend.title = paste(getSampleTracks()),
                       xlab = "Time (days)",
                       risk.table = TRUE,
                       title = "Kaplan Meier")
          }
          # if both factors are categorical
          tmp.df$SurvObj <- with(tmp.df, Surv(time = survivalParam,
                                              event = censorA.0yes.1no))
          fit <- survfit(SurvObj ~ factor1 + factor2,
                         data = tmp.df,
                         type = "kaplan-meier")
          ggsurvplot(fit,
                     data = tmp.df,
                     pval = TRUE,
                     censor = TRUE,
                     surv.median.line = "hv",
                     legend.title = paste(getSampleTracks()[1],
                                          "+",
                                          getSampleTracks()[2]),
                     xlab = "Time (days)",
                     risk.table = TRUE,
                     title = "Kaplan Meier")
        }

      }
    }



    else if(input$survivalType %in% "Cox regression"){
      print("..Running Cox Analysis..")
      plots <- list()

      for(i in 1:length(getSampleTracks())){
        this.track <- getSampleTracks()[i]

        tmp.df <- na.omit(x$sampInfo[, c(this.track,
                                         survivalParam,
                                         "censorA.0yes.1no")])
        tmp.df$censorA.0yes.1no <- as.integer(as.character(tmp.df$censorA.0yes.1no))

        names(tmp.df) <- c("factor",
                           "survivalParam",
                           "censorA.0yes.1no")

        # check if factor is continuous
        if(length(levels(as.factor(tmp.df$factor))) >= 5 && is.numeric(tmp.df$factor)){
          quants <- quantile(tmp.df$factor,
                             probs = seq(0,1,.05),
                             na.rm = T)

          # which slider bar to use to quantile
          if(i == 1){
            tmp.df$factor <- as.factor(cut(tmp.df$factor,
                                           breaks = c(quants[[1]], quants[names(quants) %in% percent(input$factorized1)], quants[[length(quants)]]),
                                           labels = c("Low", "High"),
                                           include.lowest = T))
          }
          else{
            tmp.df$factor <- as.factor(cut(tmp.df$factor, breaks = c(quants[[1]], quants[names(quants) %in% percent(input$factorized2)], quants[[length(quants)]]),
                                           labels = c("Low", "High"),
                                           include.lowest = T))
          }

          print(head(tmp.df))
          print(levels(tmp.df$factor))

          cox <- coxph(Surv(survivalParam,
                            censorA.0yes.1no) ~ factor, data = tmp.df)

          new_df <- with(tmp.df,
                         data.frame(factor = levels(tmp.df$factor)))
        }
        else{
          cox <- coxph(Surv(survivalParam,
                            censorA.0yes.1no) ~ factor, data = tmp.df)

          new_df <- with(tmp.df,
                         data.frame(factor = levels(tmp.df$factor)))
        }

        fit <- survfit(cox, newdata = new_df)
        plots[[this.track]] <- ggsurvplot(fit,
                                          data = new_df,
                                          pval = formatC(summary(cox)$logtest[['pvalue']]),
                                          censor = TRUE,
                                          surv.median.line = "hv",
                                          legend.title = paste(this.track),
                                          legend.labs = levels(tmp.df$factor),
                                          xlab = "Time (days)",
                                          risk.table = FALSE,
                                          title = paste0("Cox Regression - ", this.track))
      }
      arrange_ggsurvplots(plots,
                          print = TRUE)


    }


    # for(i in getSampleTracks()){
    #   # make the tmp.df
    #   tmp.df <- x$sampInfo[, c(i,
    #                            survivalParam,
    #                            "censorA.0yes.1no")]
    #   tmp.df$censorA.0yes.1no <- as.integer(as.character(tmp.df$censorA.0yes.1no))
    #
    #
    #   names(tmp.df) <- c("factor",
    #                      "survivalParam",
    #                      "censorA.0yes.1no")
    #
    #   # check if factor is continuous
    #   if(length(levels(as.factor(tmp.df$factor))) >= 5 && is.numeric(tmp.df$factor)){
    #     quants <- quantile(tmp.df$factor,
    #                        na.rm = T)
    #     print(quants)
    #     names(quants) <- c("0","0.25","0.5","0.75","1")
    #
    #
    #     tmp.df$factor <- as.factor(cut(tmp.df$factor, breaks = c(quants[[1]], quants[names(quants) %in% as.character(input$factorized1)], quants[[5]]),
    #                                    labels = c("Low", "High"),
    #                                    include.lowest = T))
    #
    #     print(head(tmp.df))
    #
    #     cox <- coxph(Surv(survivalParam,
    #                       censorA.0yes.1no) ~ factor, data = tmp.df)
    #
    #     new_df <- with(tmp.df,
    #                    data.frame(factor = levels(tmp.df$factor)))
    #   }
    #   else{
    #     cox <- coxph(Surv(survivalParam,
    #                       censorA.0yes.1no) ~ factor, data = tmp.df)
    #
    #     new_df <- with(tmp.df,
    #                    data.frame(factor = levels(tmp.df$factor)))
    #   }
    #
    #   fit <- survfit(cox, newdata = new_df)
    #
    #   plots[[i]] <- ggsurvplot(fit,
    #                            data = new_df,
    #                            pval = formatC(summary(cox)$logtest[['pvalue']]),
    #                            censor = TRUE,
    #                            surv.median.line = "hv",
    #                            legend.title = paste(i),
    #                            xlab = "Time (days)",
    #                            risk.table = FALSE,
    #                            title = "Cox Regression")
    # }
    # arrange_ggsurvplots(plots,
    #                     print = TRUE)

    # }


  })
  # ================================================================
  # Plot DESeq to screen
  output$volcano <- renderPlot({
    validate(need(!is.null(globals$res), message = "To run differential expression, please click the button in top right corner"))
    print("...Rendering standard volcano plot...")
    #  ------- retrieve data set and clustering results  ------
    # ------------- Run DESeq --------------

    # ------------- Standard Volcano ----------------------
    globals$res$subtype <- factor(ifelse(globals$res$labels %in% globals$genes2color, "present", "absent"))

    by_subtype <- globals$res[order(globals$res$subtype),]
    # by_subtype$log2FoldChange[by_subtype$log2FoldChange > 20] = 20
    # by_subtype$log2FoldChange[by_subtype$log2FoldChange < -20] = -20
    # by_subtype$padj[by_subtype$padj < .000000000000001] = .000000000000001
    by_subtype$both_cutoff <- by_subtype$padj < .05 & abs(by_subtype$log2FoldChange) > 1

    e <- ggplot(by_subtype, aes(x = log2FoldChange, y = -log10(padj), label = labels))
    e <- e + geom_point(aes(color = subtype))
    e <- e + geom_label(data = by_subtype[which(by_subtype$both_cutoff == T & by_subtype$subtype == "present"),],
                        hjust = 0, vjust = 0)
    e <- e + scale_color_manual(values = c("grey", "red"))
    e <- e + theme_classic() + labs(color = paste0("Present in ", input$genesets))
    e <- e + theme(legend.title = element_text(size = 16),
                   legend.text = element_text(size = 14),
                   legend.position = "bottom",
                   plot.title = element_text(size = 20, face = "bold", hjust = .5),
                   legend.background = element_rect(color = "steelblue"),
                   axis.text.x = element_text(face = "bold", size = 16),
                   axis.text.y = element_text(face = "bold", size = 16),
                   axis.title = element_text(face = "bold", size = 20),
                   plot.margin = unit(c(1,3,3,1), "cm"))
    e <- e + labs(y = expression('-log'[10]*'pval'),
                  x = expression('log'[2]*'fold change'),
                  title = (paste0("Differential Expression of genes comparing ",
                                  input$contrastA,"/", input$contrastB)))
    #e <- e + xlim(-15,15) + ylim(0,15)
    e <- e + geom_hline(yintercept = 1.3, linetype = "dotted") + geom_vline(xintercept=c(-2,2), linetype="dotted")
    removeNotification(id = "DESeq_notification")

    e
  })

  # ==========================================================
  # output link to save rendered pdf heatmaps
  output$pdflink <- downloadHandler(
    filename = "myplot.pdf",
    content = function(file) {
      writeHeatmapToPdf()
      file.copy("myplot.pdf", file)
    })
  # ==========================================================
  # output link to save table
  output$tablelink <- downloadHandler(
    filename = "mytable.tsv",
    content = function(file) {
      writeHeatmapToTable()
      file.copy("mytable.tsv", file)
    })

}
# ================================================================
# ------------- Run app ------------------------------------------
# ================================================================
shinyApp(ui = ui, server = server)
