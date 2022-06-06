#' Download from GEO
#'
#' @export
#' @import GEOquery
#' @import Biobase
#' @import org.Hs.eg.db
#' @import openxlsx
#' @import plyr

pull_Seino_from_GEO <- function() {

  # library(GEOquery)
  # library(Biobase)
  # library(org.Hs.eg.db)
  # library(openxlsx)
  # library(plyr)

  ## =============================
  # Read expression file and initial parse sample info
  ## =============================
  gset <- getGEO(filename = system.file("extdata/Seino",
                                        "GSE107610_series_matrix.txt",
                                        package = "pdacR"),
                 GSEMatrix =TRUE,
                 AnnotGPL=TRUE,
                 getGPL=FALSE)

  samps <-gset@phenoData
  sampInfo <- samps@data
  sampInfo <- sampInfo[,c("title",
                          "geo_accession",
                          "type",
                          "extract_protocol_ch1",
                          "disease state:ch1")]
  names(sampInfo)[names(sampInfo)=="title"] <- "submitted_sample_id"

  for (name in names(sampInfo)){
    sampInfo[[name]] <- as.factor(sampInfo[[name]])
  }

  ## =============================
  # Read expression file
  ## =============================
  ex <- exprs(gset)
  ex <- as.data.frame(ex)

  ## =============================
  # Organize feature info
  ## =============================
  featInfo <- data.frame(SYMBOL = rownames(ex),
                         filler = 1:nrow(ex))

  feature_file <- system.file("extdata/Seino",
                              "GSE107610_gene_annotations.xlsx",
                              package = "pdacR")

  feature_table <- read.xlsx(feature_file)

  feature_table$SYMBOL <- feature_table$ID

  tmp <- merge(featInfo,
               feature_table,
               by = "SYMBOL",
               all.x = TRUE)

  important_features <- tmp[, c("SYMBOL",
                                "Ensembl",
                                "Gene.Symbol",
                                "Entrez.Gene",
                                "SwissProt")]

  important_features$Seino_annotation <- important_features$SYMBOL
  important_features$SYMBOL <- important_features$Gene.Symbol
  important_features$ENTREZID <- important_features$Entrez.Gene
  important_features$ENSEMBLID <- important_features$Ensembl

  important_features$Gene.Symbol <- NULL
  important_features$Entrez.Gene <- NULL
  important_features$Ensembl <- NULL

  featInfo <- important_features

  ## =============================
  # Clean gene SYMBOLs
  ## =============================
  SYMBOL <- sapply(X = featInfo$SYMBOL,
                   FUN = function(x){
                     y <- strsplit(x = x,
                                   split = "///")
                     y <- y[[1]][1]
                     SYMBOL <- y
                   })
  featInfo$SYMBOL <- unname(SYMBOL)

  entrez <- sapply(X = featInfo$ENTREZID,
                   FUN = function(x){
                     y <- strsplit(x = x,
                                   split = "///")
                     y <- y[[1]][1]
                     entrez <- y
                   })
  featInfo$ENTREZID <- unname(entrez)

  ## =============================
  # Aggregate probe IDs to SYMBOLS by sum
  ## =============================
  tmp <- aggregate(x = ex,
                   by = featInfo[,c(1,4)],
                   FUN = function(x) sum(x))
  featInfo <- data.frame(SYMBOL = tmp[[1]],
                         ENTREZID = tmp[[2]])
  ex  <- tmp[-1,c(-1,-2)]
  tmp <- featInfo[-1,]
  featInfo <- tmp

  ## =============================
  # Add metadata
  ## =============================
  metadata <- list(log.transformed = FALSE,
                   reference = "Seino T et al, Cell Stem Cell, 2018, PMID:29337182",
                   accession = "GEO: GSE107610",
                   description = "A library of 2 normal, 39 PDAC, 6 normal-like PDAC, and 10 engineered organoids",
                   survivalA = "None",
                   survivalB = "None",
                   default_selections = list(filter_column = "cancer_normal",
                                             filter_levels = c("normal",    "normalLike", "sgRNA-treated")))

  ## =============================
  # Add sample information not in GEO
  ## =============================
  diff_meth_file <- system.file("extdata/Seino",
                                "GSE107610_methylation.xlsx",
                                package = "pdacR")

  exp_cases_file <- system.file("extdata/Seino",
                                "GSE107610_experiment_cases.xlsx",
                                package = "pdacR")

  wnt_subtypes_file <- system.file("extdata/Seino",
                                   "GSE107610_wnt_subtypes.xlsx",
                                   package = "pdacR")

  diff_meth <- read.xlsx(diff_meth_file)
  exp_cases <- read.xlsx(exp_cases_file)
  wnt_subtypes <- read.xlsx(wnt_subtypes_file)


  sampInfo[, c("No.",
               "Gene.expression.array",
               "Whole.exome.sequence",
               "Methylation.array",
               "NL.counterpart",
               "Blood",
               "Engineered")]  <- NA

  sampInfo[c(1:4,12:28,32:57),
           c("No.",
             "Gene.expression.array",
             "Whole.exome.sequence",
             "Methylation.array",
             "NL.counterpart",
             "Blood",
             "Engineered")] <- exp_cases[c(1:5,
                                           7:8,
                                           11:45,
                                           47:51),]

  wnt_subtypes$No. <- wnt_subtypes$sample

  sampInfo <- join(sampInfo,
                   wnt_subtypes,
                   by = "No.")

  featInfo_methylation <- data.frame(SYMBOL = diff_meth$Gene.Symbol,
                                     target_ID = diff_meth$TargetID,
                                     pval = diff_meth$p_value,
                                     FDR = diff_meth$FDR)

  diff_meth[, c("Gene.Symbol",
                "TargetID",
                "p_value",
                "FDR")] <- NULL

  ## =============================
  # Organize dataset
  ## =============================
  Seino_GEO_array <- list(sampInfo = sampInfo,
                          featInfo = featInfo,
                          featInfo_meth = featInfo_methylation,
                          ex = ex,
                          diff_meth = diff_meth,
                          metadata = metadata)

  Seino_GEO_array$sampInfo <- Seino_GEO_array$sampInfo[, -which(names(Seino_GEO_array$sampInfo) %in% c("type",
                                                                                                       "extract_protocol_ch1",
                                                                                                       "sample"))]

  names(Seino_GEO_array$sampInfo) <- c("submitted_sample_id",
                                       "geo_accession",
                                       "cancer.normal",
                                       "patient_number",
                                       "expression.array",
                                       "whole.exome.seq",
                                       "meth.array",
                                       "NL.counterpart",
                                       "blood",
                                       "engineered",
                                       "wnt_stubtype")

  ## =============================
  # Fix sample annotation for easy use in GUI
  ## =============================
  Seino_GEO_array$sampInfo$cancer.normal <- as.character(Seino_GEO_array$sampInfo$cancer.normal)
  engineered <- grep(Seino_GEO_array$sampInfo$submitted_sample_id,
                     pattern = "engineered")
  normal_like <- grep(Seino_GEO_array$sampInfo$submitted_sample_id,
                      pattern = "normal like")
  more_engineered <- grep(Seino_GEO_array$sampInfo$submitted_sample_id,
                          pattern = "RNA")

  Seino_GEO_array$sampInfo$cancer.normal[engineered] <- "engineered"
  Seino_GEO_array$sampInfo$cancer.normal[normal_like] <- "normalLike"
  Seino_GEO_array$sampInfo$cancer.normal[more_engineered] <- "sgRNA-treated"

  Seino_GEO_array$sampInfo$cancer_normal <- as.factor(Seino_GEO_array$sampInfo$cancer.normal)
  Seino_GEO_array$sampInfo$cancer.normal <- NULL


  ## =============================
  # Consensus clustering to find tumor and stroma subtypes
  ## =============================
  #
  #   # Tumor
  #   # ----------------------
  #   dataset <- Seino_GEO_array
  #
  #   sampleset <- which(dataset$sampInfo$cancer.normal %in% "cancer")
  #   tmp.k <- 2
  #   tmp.ncusts <- 2
  #
  #   featureset <- which(dataset$featInfo$SYMBOL %in%
  #                         c(as.character(pdacR::gene_lists$Moffitt.Classical.25),
  #                           as.character(pdacR::gene_lists$Moffitt.Basal.25)))
  #
  #   smallx <- t(scale(t(dataset$ex[featureset,sampleset])))
  #
  #   sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx),
  #                                                            seed = 1234,
  #                                                            pFeature = 0.8,
  #                                                            pItem = 0.8,
  #                                                            maxK = 6,
  #                                                            reps=200,
  #                                                            distance="pearson",
  #                                                            clusterAlg="hc")[[tmp.k]]$consensusTree
  #
  #   tmp.cluster <- c("basal","classical")[cutree(tree = sampletree, k = 2)]
  #   dataset$sampInfo$MoffittTumor <- NA
  #   dataset$sampInfo$MoffittTumor[sampleset] <- tmp.cluster
  #
  #   ColSideColors <-  pdacR::getSideColors(sampInfo = dataset$sampInfo[sampleset,],
  #                                          sampleTracks = c("MoffittTumor"),
  #                                          colorlists = list(c("orange", "blue")),
  #                                          drop.levels = TRUE)
  #
  #   RowSideColors <-  pdacR::getSideColors(sampInfo = data.frame(basal =dataset$featInfo$SYMBOL[featureset] %in%
  #                                                                  pdacR::gene_lists$Moffitt.Basal.25,
  #                                                                classical =dataset$featInfo$SYMBOL[featureset] %in%
  #                                                                  pdacR::gene_lists$Moffitt.Classical.25),
  #                                          sampleTracks = c("basal",
  #                                                           "classical"),
  #                                          colorlists = list(c=c("white","orange"),
  #                                                            b=c("white","blue")))
  #   pdacR::heatmap.3(x = smallx,
  #                    scale="row",
  #                    labRow = dataset$featInfo$SYMBOL[featureset],
  #                    col = colorRampPalette(c("blue", "white", "red"))(n = 299),
  #                    Colv = as.dendrogram(sampletree),
  #                    Rowv = TRUE,
  #                    distfun = function(x) as.dist((1-cor(t(x)))/2),
  #                    ColSideColors = ColSideColors$SideColors,
  #                    ColSideColorsSize = 6,
  #                    RowSideColorsSize = 6,
  #                    RowSideColors = t(RowSideColors$SideColors),
  #                    margins = c(5,20))
  #   legend(xy.coords(x=.90,y=1),
  #          legend=c(ColSideColors$text),
  #          fill=c(ColSideColors$colors),
  #          border=FALSE, bty="n",
  #          y.intersp = 0.9, cex=0.5)
  #
  #   # Stroma
  #   # ----------------------
  #   dataset$sampInfo$MoffittStroma <- NA

  ## =============================
  # Save dataset
  ## =============================

  Seino_GEO_array$sampInfo = Seino_GEO_array$sampInfo[,order(colnames(Seino_GEO_array$sampInfo))]
  Seino_GEO_array$sampInfo = Seino_GEO_array$sampInfo[,c(which(colnames(Seino_GEO_array$sampInfo)=="submitted_sample_id"),
                                                         which(colnames(Seino_GEO_array$sampInfo)!="submitted_sample_id"))]

  saveRDS(Seino_GEO_array,
          file = "./data/Seino_GEO_array.rds",
          compress = T)
  return(NULL)
}
