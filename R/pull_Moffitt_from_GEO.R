#' Download Moffitt array from GEO and parse
#'
#' @export
#' @import GEOquery
#' @import Biobase
#' @import org.Hs.eg.db
#' @import openxlsx
#' @import plyr

pull_Moffitt_from_GEO <- function() {

  # library(GEOquery)
  # library(Biobase)
  # library(org.Hs.eg.db)
  # library(openxlsx)
  # library(plyr)

  ## =============================
  # Pull Moffitt GEO dataset from GEO
  ## =============================
  gset <- getGEO(filename = system.file("extdata/Moffitt_GEO",
                                        "GSE71729_series_matrix.txt.gz",
                                        package = "pdac"),
                 GSEMatrix =TRUE,
                 AnnotGPL=FALSE,
                 getGPL=FALSE)

  ## =============================
  # Parse and organize sample metadata
  ## =============================
  samps <-gset@phenoData
  sampInfo <- samps@data
  sampInfo <- sampInfo[,c("title","geo_accession",
                          "characteristics_ch2",
                          "characteristics_ch2.1","characteristics_ch2.2",
                          "characteristics_ch2.3","characteristics_ch2.4",
                          "characteristics_ch2.5")]
  names(sampInfo)[names(sampInfo)=="title"] <- "submitted_sample_id"
  sampInfo$specimen_type <- "cell line"
  sampInfo$location <- as.character(NA)
  sampInfo$moffittTumor <-  as.character(NA)
  sampInfo$moffittStroma <-  as.character(NA)
  sampInfo$survival_months <- NA
  sampInfo$censor <-  as.character(NA)
  for(characteristics in c("characteristics_ch2","characteristics_ch2.1",
                           "characteristics_ch2.2","characteristics_ch2.3",
                           "characteristics_ch2.4","characteristics_ch2.5")){
    sampInfo$censor[sampInfo[,characteristics] == "death_event_1death_0censor: 0"] <- "censor"
    sampInfo$censor[sampInfo[,characteristics] == "death_event_1death_0censor: 1"] <- "death"
    sampInfo$moffittTumor[sampInfo[,characteristics] == "tumor_subtype_0na_1classical_2basal: 1"] <- "classical"
    sampInfo$moffittTumor[sampInfo[,characteristics] == "tumor_subtype_0na_1classical_2basal: 2"] <- "basal"
    sampInfo$moffittStroma[sampInfo[,characteristics] == "stroma_subtype_0na_1low_2normal_3activated: 1"] <- "low"
    sampInfo$moffittStroma[sampInfo[,characteristics] == "stroma_subtype_0na_1low_2normal_3activated: 2"] <- "normal"
    sampInfo$moffittStroma[sampInfo[,characteristics] == "stroma_subtype_0na_1low_2normal_3activated: 3"] <- "activated"
    for( i in grep(x = sampInfo[,characteristics],pattern = "survival_months")){
      tmp <- as.character(sampInfo[i,characteristics])
      tmp <- strsplit(x=tmp,split=": ")
      sampInfo$survival_months[i] <- as.numeric(tmp[[1]][2])
    }
    for( i in grep(x = sampInfo[,characteristics],pattern = "cell line/tissue")){
      tmp <- as.character(sampInfo[i,characteristics])
      tmp <- strsplit(x=tmp,split=": ")
      sampInfo$location[i] <- (tmp[[1]][2])
    }
    for( i in grep(x = sampInfo[,characteristics],pattern = "tissue type:")){
      tmp <- as.character(sampInfo[i,characteristics])
      tmp <- strsplit(x=tmp,split=": ")
      sampInfo$specimen_type[i] <- (tmp[[1]][2])
    }
  }
  sampInfo <- sampInfo[,-which(names(sampInfo) %in% c("characteristics_ch2","characteristics_ch2.1",
                           "characteristics_ch2.2","characteristics_ch2.3",
                           "characteristics_ch2.4","characteristics_ch2.5"))]
  sampInfo$specimen_type <- as.factor(sampInfo$specimen_type)
  sampInfo$location <- as.factor(sampInfo$location)
  sampInfo$moffittTumor <- as.factor(sampInfo$moffittTumor)
  sampInfo$moffittStroma <- as.factor(sampInfo$moffittStroma)
  sampInfo$censor <- as.factor(sampInfo$censor)


  ## =============================
  #  Parse expression dataframe
  ## =============================
  ex <- exprs(gset)

  ## =============================
  # Parse gene labels
  ## =============================
  featInfo <- data.frame(SYMBOL = rownames(ex))

  ## =============================
  # Include dataset metadata
  ## =============================
  metadata <- list(log.transformed = TRUE,
                   reference = "Moffitt RA et al, Nat Gen, 2015, PMID:26343385",
                   accession = "GEO: GSE71729",
                   description = "145 primary PDAC, 61 metastatic PDAC, 17 cell lines, 46 normal pancreas, 88 distant site normal tissue, microarray",
                   survivalA = "overall survival days",
                   survivalB = "None")


  ## =============================
  # Add sample info not on GEO
  ## =============================
  filename <- system.file("extdata/Moffitt_GEO",
                          "Moffitt_extended.xlsx",
                          package = "pdac")
  extra.data <- data.frame(lapply(read.xlsx(filename),factor))
  sampInfo <- join(sampInfo,extra.data)

  ## =============================
  # Correct data types
  ## =============================
  sampInfo$purity <- as.numeric(as.character(sampInfo$purity))
  summary(sampInfo)

  ## =============================
  # Preparation for classifier training
  ## =============================
  sampInfo$tumor.classifier.training <- FALSE
  sampInfo$stroma.classifier.training <- FALSE
  sampInfo[["cluster.MT"]] <- as.character(NA)
  sampInfo[["cluster.MS"]] <- as.character(NA)
  sampInfo$cluster.MT[sampInfo$moffittTumor %in% "basal"] <- "basal"
  sampInfo$cluster.MT[sampInfo$moffittTumor %in% "classical"] <- "classical"
  sampInfo$cluster.MS[sampInfo$moffittStroma %in% "activated"] <- "activated"
  sampInfo$cluster.MS[sampInfo$moffittStroma %in% "normal"] <- "normal"
  sampInfo$stroma.classifier.training[sampInfo$moffittStroma %in% c("activated","normal")] <- TRUE
  sampInfo$tumor.classifier.training[sampInfo$moffittTumor %in% c("basal","classical")] <- TRUE
  sampInfo$cluster.MT <- as.factor(sampInfo$cluster.MT)
  sampInfo$cluster.MS <- as.factor(sampInfo$cluster.MS)
  sampInfo$MoffittTumor <- sampInfo$cluster.MT
  sampInfo$MoffittStroma <- sampInfo$cluster.MS
  sampInfo$cluster.MT <- NULL
  sampInfo$cluster.MS <- NULL

  ## =============================
  # Wrangle survival data for shiny app
  ## =============================

  sampInfo$survivalA <- sampInfo$survival_months * (365/12)
  sampInfo$censorA.0yes.1no <- factor(sampInfo$censor,
                                      levels = c("censor","death", NaN),
                                      labels = c("0","1", "NA"))

  sampInfo$survival_months <- NULL
  sampInfo$censor <- NULL

  print(names(sampInfo))

  ## =============================
  # Set default parameters for app
  ## =============================
  default.tracks <- list(sample.tracks = c("specimen_type"),
                         gene.sets = c("Moffitt.Basal.25","Moffitt.Classical.25"),
                         filter.out.track = c("specimen_type"),
                         filter.out.levels = c("Normal"))

  ## =============================
  # Save parsed dataset
  ## =============================
  Moffitt_GEO_array <- list(sampInfo=sampInfo,featInfo=featInfo,ex=ex,metadata=metadata,default.tracks=default.tracks)
  save(list = c("Moffitt_GEO_array"),
       file = "./data/Moffitt_GEO_array.RData")

  return(NULL)
}
