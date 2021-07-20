#' Run our shiny app
#'
#' @export
#' @import bioDist
#' @import ConsensusClusterPlus
#' @import foreign
#' @import ggplot2
#' @import ggpubr
#' @import gplots
#' @import limma
#' @import nlme
#' @import preprocessCore
#' @import RColorBrewer
#' @import reshape2
#' @import Rtsne
#' @import scales
#' @import shiny
#' @import shinyjs
#' @import stringr
#' @import survival
#' @import survminer

pdacShiny <- function() {
  appDir <- system.file("shiny", package = "pdacR")
  runApp(appDir, display.mode = "normal")
}

