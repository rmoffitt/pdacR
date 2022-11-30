#' Parse data from tables
#'
#' @export
#' @import org.Hs.eg.db
#' @import stringr
#' @import plyr
#' @import GEOquery


pull_Olive_from_GEO <- function() {
  # ------------------------------------------------------------------------------------
  # library(org.Hs.eg.db)
  # library(stringr)
  # library(plyr)
  # library(GEOquery)

  ## =============================
  ## Pull in sample information from GEO
  ## =============================
  gset <- getGEO(filename = system.file("extdata/Olive",
                                        "GSE93326-GPL11154_series_matrix.txt.gz",
                                        package = "pdacR"),
                 GSEMatrix =TRUE,
                 AnnotGPL=TRUE,
                 getGPL=FALSE)

  ## =============================
  ## Process sample information and add sample information from Olive 2019 supplementary paper in Gut
  ## =============================
  samps <- gset@phenoData
  sampInfo <- samps@data

  for (name in names(sampInfo)){
    if(is.character(sampInfo[[name]])){
      sampInfo[[name]] <- as.factor(sampInfo[[name]])
    }
  }

  summary(sampInfo)

  sampInfo <- sampInfo[,c("title",
                          "geo_accession",
                          "library_method:ch1",
                          "data_processing.3",
                          "patient_id:ch1",
                          "tumor:ch1")]

  names(sampInfo)[names(sampInfo)=="title"] <- "submitted_sample_id"

  print(summary(sampInfo))


  sampInfo$data_processing.3 <- NULL

  sampInfo$sampleID <- sapply(X = sampInfo$submitted_sample_id,
                              FUN = function(x){
                                y <- strsplit(x = as.character(x),
                                              split = " ")[[1]][1]
                                return(y)
                              })

  sampInfo$submitted_sample_id <- NULL

  sampInfo$patientID <- sapply(X = as.character(sampInfo$sampleID),
                               FUN = function(x){
                                 if(length(strsplit(x,split = "")[[1]]) == 7){
                                   y <- gsub(pattern = "_S|_E|_B",
                                             replacement = "_00",
                                             x)
                                   return(y)
                                 } else if(length(strsplit(x,split = "")[[1]]) == 8){
                                   y <- gsub(pattern = "_S|_E|_B",
                                             replacement = "_0",
                                             x)
                                   return(y)
                                 } else if(length(strsplit(x,split = "")[[1]]) == 9){
                                   y <- gsub(pattern = "_S|_E|_B",
                                             replacement = "_",
                                             x)
                                   return(y)
                                 }
                               })

  names(sampInfo) <- c(names(sampInfo)[1],
                       "LibraryPlatform",
                       "PatientID_original",
                       "TissueType", names(sampInfo)[5:6])

  sampInfo$sampleID <- as.factor(sampInfo$sampleID)
  sampInfo$patientID <- as.factor(sampInfo$patientID)



  supp_sampinfo <- system.file("extdata/Olive",
                               "gutjnl-2018-317706_supplementary_clinical_data.txt",
                               package = "pdacR")

  supp <- read.table(supp_sampinfo,
                     header = TRUE,
                     skip = 2,
                     sep = "\t")
  names(supp) <- c("patientID", names(supp)[2:length(names(supp))])

  print(head(supp))



  ## =============================
  ## Pull in and process expression data from GEO raw counts text file
  ## =============================
  expression_file <- system.file("extdata/Olive",
                                 "GSE93326_LCM-RNA-Seq-RawCounts.txt.gz",
                                 package = "pdacR")

  ex <- as.data.frame(read.table(expression_file,
                                 header = TRUE))

  ## =============================
  ## Process featInfo from Olive expression data
  ## =============================

  featinfo <- data.frame(SYMBOL = ex$Symbol)
  ex$Symbol <- NULL

  featinfo$ENSEMBL <- as.factor(AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                                      keys = as.character(featinfo$SYMBOL),
                                                      keytype = "SYMBOL",
                                                      column = "ENSEMBL",
                                                      multiVals = "first"))

  head(featinfo)

  ## =============================
  ## Organize expression and sample information
  ## =============================

  matched_samples <- intersect(names(ex), sampInfo$sampleID)

  sampInfo_v1 <- sampInfo[match(matched_samples, sampInfo$sampleID),]
  ex_v1 <- ex[, (match(matched_samples, names(ex)))]

  sampInfo_v2 <- join(sampInfo_v1,
                      supp,
                      by = "patientID")


  addtl_stroma <- names(ex)[which(!names(ex) %in% matched_samples)]
  addtl_stroma_clean <- gsub(pattern = "S", replacement = "", addtl_stroma)

  addtl_stroma_set <- data.frame(geo_accession = NA,
                                 LibraryPlatform = NA,
                                 PatientID_original = NA,
                                 TissueType = "stroma",
                                 sampleID = addtl_stroma,
                                 patientID = addtl_stroma_clean)

  addtl_stroma_set <- join(addtl_stroma_set,
                           supp,
                           by = "patientID")

  sampInfo_v3 <- rbind(sampInfo_v2, addtl_stroma_set)

  rm(sampInfo_v1, sampInfo_v2)

  sampInfo_v3 <- sampInfo_v3[, c("geo_accession",
                                 "LibraryPlatform",
                                 "patientID",
                                 "sampleID",
                                 "TissueType",
                                 "Age",
                                 "Race",
                                 "Gender",
                                 "Primary_Site",
                                 "Location",
                                 "AJCC_Stage",
                                 "Surgical_Margins",
                                 "Histological_diagnosis",
                                 "Tumor_grade",
                                 "Solid_pattern",
                                 "Mucinous_features",
                                 "Clear_cell_features",
                                 "Micropapillary_growth",
                                 "Features_Of_Possible_Concomittant_IPMN")]

  sampInfo_v3$TissueType[which(sampInfo_v3$TissueType %in% NA)] <- "stroma"


  addtl_stroma_ex <- ex[, match(addtl_stroma, names(ex))]

  ex_v2 <- cbind(ex_v1, addtl_stroma_ex)

  dim(sampInfo_v3)
  dim(ex_v2)

  summary(sampInfo_v3)
  print(data.frame(expression.name = names(ex_v2),
                   sampinfo_name = sampInfo_v3$patientID))


  metadata <- list(data.processing = "HTSeq for read counts",
                   log.transformed = FALSE,
                   reference = "Maurer C et al, Gut, 2019, PMID:30658994",
                   accession = "GEO: GSE93326",
                   description = "66 matched laser-capture-dissected stroma and epithelium samples, with a subset of triplicate bulk, stroma, epithelium matches and 57 additional unmatched stroma",
                   survivalA = "None",
                   survivalB = "None",
                   default_selections = list(filter_column = "TissueType",
                                             filter_levels = c("TissueType:stroma"),
                                             sampleTracks = c("Features_Of_Possible_Concomittant_IPMN","TissueType")))

  dataset <- list(ex_v2,
                  sampInfo_v3,
                  featinfo,
                  metadata)

  names(dataset) <- c("ex",
                      "sampInfo",
                      "featInfo",
                      "metadata")

  ## =============================
  # Consensus clustering to find tumor and stroma subtypes - DEPRECATED WITH PURIST
  ## =============================

  # # Tumor
  # # ----------------------
  # sampleset <- which(dataset$sampInfo$TissueType %in% "epithelium")
  # tmp.k <- 2
  # tmp.ncusts <- 2
  #
  # featureset <- which(dataset$featInfo$SYMBOL %in%
  #                       c(as.character(pdacR::gene_lists$Moffitt.Classical.25),
  #                         as.character(pdacR::gene_lists$Moffitt.Basal.25)))
  #
  # smallx <- t(scale(t(dataset$ex[featureset,sampleset])))
  #
  # sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx),
  #                                                          seed = 1234,
  #                                                          pFeature = 0.8,
  #                                                          pItem = 0.8,
  #                                                          maxK = 6,
  #                                                          reps=200,
  #                                                          distance="pearson",
  #                                                          clusterAlg="hc")[[tmp.k]]$consensusTree
  #
  # tmp.cluster <- c("basal","classical")[cutree(tree = sampletree, k = 2)]
  # dataset$sampInfo$MoffittTumor <- NA
  # dataset$sampInfo$MoffittTumor[sampleset] <- tmp.cluster
  #
  # ColSideColors <-  pdacR::getSideColors(sampInfo = dataset$sampInfo[sampleset,],
  #                                 sampleTracks = c("MoffittTumor"),
  #                                 colorlists = list(c("orange", "blue")),
  #                                 drop.levels = TRUE)
  #
  # RowSideColors <-  pdacR::getSideColors(sampInfo = data.frame(basal =dataset$featInfo$SYMBOL[featureset] %in%
  #                                                         pdacR::gene_lists$Moffitt.Basal.25,
  #                                                       classical =dataset$featInfo$SYMBOL[featureset] %in%
  #                                                         pdacR::gene_lists$Moffitt.Classical.25),
  #                                 sampleTracks = c("basal",
  #                                                  "classical"),
  #                                 colorlists = list(c=c("white","orange"),
  #                                                   b=c("white","blue")))
  # pdacR::heatmap.3(x = smallx,
  #           scale="row",
  #           labRow = dataset$featInfo$SYMBOL[featureset],
  #           col = colorRampPalette(c("blue", "white", "red"))(n = 299),
  #           Colv = as.dendrogram(sampletree),
  #           Rowv = TRUE,
  #           distfun = function(x) as.dist((1-cor(t(x)))/2),
  #           ColSideColors = ColSideColors$SideColors,
  #           ColSideColorsSize = 6,
  #           RowSideColorsSize = 6,
  #           RowSideColors = t(RowSideColors$SideColors),
  #           margins = c(5,20))
  # legend(xy.coords(x=.90,y=1),
  #        legend=c(ColSideColors$text),
  #        fill=c(ColSideColors$colors),
  #        border=FALSE, bty="n",
  #        y.intersp = 0.9, cex=0.5)
  #
  # # Stroma
  # # ----------------------
  # geneMeans <- rowMeans(dataset$ex)
  # genesToDelete <- which(geneMeans < .01)
  #
  # dataset$ex <- log2(1+dataset$ex[-genesToDelete,])
  # dataset$featInfo <- dataset$featInfo[-genesToDelete,]
  #
  # sampleset <- which(dataset$sampInfo$TissueType %in% "stroma")
  # tmp.k <- 2
  # tmp.ncusts <- 2
  #
  # featureset <- which(dataset$featInfo$SYMBOL %in%
  #                       c(as.character(pdacR::gene_lists$Moffitt.Normal.25),
  #                         as.character(pdacR::gene_lists$Moffitt.Activated.25)))
  #
  # smallx <- t(scale(t(dataset$ex[featureset,sampleset])))
  #
  # sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx),
  #                                                          seed = 1234,
  #                                                          pFeature = 0.8,
  #                                                          pItem = 0.8,
  #                                                          maxK = 6,
  #                                                          reps=200,
  #                                                          distance="pearson",
  #                                                          clusterAlg="hc")[[tmp.k]]$consensusTree
  #
  # tmp.cluster <- c("normal","activated")[cutree(tree = sampletree, k = 2)]
  # dataset$sampInfo$MoffittStroma <- NA
  # dataset$sampInfo$MoffittStroma[sampleset] <- tmp.cluster
  #
  # ColSideColors <-  pdacR::getSideColors(sampInfo = dataset$sampInfo[sampleset,],
  #                                 sampleTracks = c("MoffittStroma"),
  #                                 colorlists = list(c("brown", "lightblue")),
  #                                 drop.levels = TRUE)
  #
  # RowSideColors <-  pdacR::getSideColors(sampInfo = data.frame(normal =dataset$featInfo$SYMBOL[featureset] %in%
  #                                                         pdacR::gene_lists$Moffitt.Normal.25,
  #                                                       activated =dataset$featInfo$SYMBOL[featureset] %in%
  #                                                         pdacR::gene_lists$Moffitt.Activated.25),
  #                                 sampleTracks = c("normal",
  #                                                  "activated"),
  #                                 colorlists = list(c=c("white","lightblue"),
  #                                                   b=c("white","brown")))


  # pdacR::heatmap.3(x = smallx,
  #           scale="row",
  #           labRow = dataset$featInfo$SYMBOL[featureset],
  #           col = colorRampPalette(c("blue", "white", "red"))(n = 299),
  #           Colv = as.dendrogram(sampletree),
  #           Rowv = TRUE,
  #           distfun = function(x) as.dist((1-cor(t(x)))/2),
  #           ColSideColors = ColSideColors$SideColors,
  #           ColSideColorsSize = 6,
  #           RowSideColorsSize = 6,
  #           RowSideColors = t(RowSideColors$SideColors),
  #           margins = c(5,20))
  # legend(xy.coords(x=.90,y=1),
  #        legend=c(ColSideColors$text),
  #        fill=c(ColSideColors$colors),
  #        border=FALSE, bty="n",
  #        y.intersp = 0.9, cex=0.5)


  Olive_2019 <- dataset
  Olive_2019$sampInfo = Olive_2019$sampInfo[,order(colnames(Olive_2019$sampInfo))]
  Olive_2019$sampInfo = Olive_2019$sampInfo[,c(which(colnames(Olive_2019$sampInfo)=="patientID"),
                                               which(colnames(Olive_2019$sampInfo)!="patientID"))]
  ## =============================
  # Save dataset
  ## =============================
  # ------------------------------------------------------------------------------------
  saveRDS(Olive_2019,
          file = "./data/Olive_2019.rds",
          compress=T)
  return(NULL)
}
