#' Parse Moffitt supplemental RNAseq corresponding to GEO array
#'
#' @export

parse_Moffitt_supplement <- function() {


  ## =============================
  # Read expression file and parse expression dataframe, along with gene names
  ## =============================
  expression_file <- system.file("extdata/Moffitt_supplement",
                                 "Moffitt_S2.txt",
                                 package = "pdacR")
  input.data <- read.table(file = expression_file,
                           header = FALSE,
                           skip = 7,
                           sep = "\t")

  input.head <- read.table(file = expression_file,
                           header = FALSE,
                           nrows = 6,
                           fill = TRUE,
                           sep = "\t")

  input.head <- data.frame(t(input.head[,-(1:3)]),
                           stringsAsFactors = FALSE)
  names(input.head) <- as.character(input.head[1,])
  input.head <- input.head[-1,]
  input.head <- input.head[-length(input.head[[1]]),]

  sampInfo <- input.head

  featInfo <- data.frame(SYMBOL=input.data[,1],
                         species=input.data[,2])

  ex <- (input.data[,c(-(1:4),-(length(input.data)))])

  ## =============================
  # Parse and organize sample information
  ## =============================
  sampInfo <- data.frame(sampInfo)
  sampInfo$sample_type <- ""
  sampInfo$sample_type[sampInfo$CAF == "1"] <-"CAF"
  sampInfo$sample_type[sampInfo$PDX == "1"] <-"PDX"
  sampInfo$sample_type[sampInfo$CellLine == "1"] <-"CellLine"
  sampInfo$sample_type[sampInfo$Primary == "1"] <-"Primary"
  sampInfo$sample_type <- factor(sampInfo$sample_type)
  sampInfo <- sampInfo[,names(sampInfo) %in% c("publicID","KRAS.mutation","sample_type")]
  sampInfo$KRAS.mutation[sampInfo$KRAS.mutation %in% ""] <- NA
  sampInfo$KRAS.mutation <- as.factor(sampInfo$KRAS.mutation)
  summary(sampInfo)

  ## =============================
  # Add dataset metadata for reference
  ## =============================
  metadata <- list(log.transformed = FALSE,
                   reference = "Moffitt RA et al, Nat Gen, 2015, PMID:26343385",
                   accession = "NA: supplementary file from reference",
                   description = "15 primary tumor, 37 patient derived xenograft, 3 cell lines, 6 cancer associated fibroblast lines, RNAseq",
                   survivalA = "None",
                   survivalB = "None",
                   exp.type = "RNAseq",
                   default_selections = list(filter_column = "sample_type",
                                             filter_levels = c("CAF"),
                                             sampleTracks = c("sample_type")))

  ## =============================
  # Add gene species
  ## =============================
  featInfo$species <- factor(featInfo$species,labels = c("Hs","Mm"))

  ## =============================
  # Aggregation by species and symbol
  ## =============================
  tmp <- aggregate(x = ex,
                   by = featInfo,
                   FUN = function(x) sum(x))
  featInfo <- data.frame(SYMBOL = tmp[[1]],
                         species = tmp[[2]])
  ex  <- tmp[,c(-1,-2)]

  ## =============================
  # Prepare dataset for classifier training
  ## =============================
  # sampInfo$tumor.classifier.training <- FALSE
  # sampInfo$stroma.classifier.training <- FALSE
  # sampInfo[["MoffittTumor"]] <- as.character(NA)
  # sampInfo[["MoffittStroma"]] <- as.character(NA)
  #
  # # Call consensus subtypes for dataset
  # #############
  # # vvvvvvv
  # sampleset <- which(sampInfo$sample_type %in% "PDX")
  # featureset <- which(featInfo$SYMBOL %in%
  #                       c(as.character(pdacR::gene_lists$Moffitt.Basal.25),
  #                         as.character(pdacR::gene_lists$Moffitt.Classical.25)))
  # # ^^^^^^^
  # smallx <- t(scale(scale=FALSE,
  #                   center=TRUE,
  #                   x=t(log2(1+ex[featureset,sampleset]))))
  # sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = smallx,
  #                                                          seed = 1234,
  #                                                          pFeature = 0.8,
  #                                                          pItem = 0.8,
  #                                                          maxK = 6,
  #                                                          reps=200,
  #                                                          distance="pearson",
  #                                                          clusterAlg="kmdist")[[2]]$consensusTree
  # tmp.cluster <- c("classical","basal")[cutree(tree = sampletree, k = 2)]
  # sampInfo$MoffittTumor[sampleset] <- tmp.cluster
  # ggplot(data = data.frame(classical = (colMeans(ex[featInfo$SYMBOL %in%
  #                                                pdacR::gene_lists$Moffitt.Classical.25,sampleset])),
  #                          basal = (colMeans(ex[featInfo$SYMBOL %in%
  #                                                pdacR::gene_lists$Moffitt.Basal.25,sampleset])),
  #                          MoffittTumor = sampInfo$MoffittTumor[sampleset]),
  #        aes(x = basal,y = classical, color = MoffittTumor)) +
  #   geom_point(size=2)
  #############

  sampInfo = sampInfo[,order(colnames(sampInfo))]
  sampInfo = sampInfo[,c(which(colnames(sampInfo) == "publicID"),
                         which(colnames(sampInfo) != "publicID"))]
  ## =============================
  # Save dataset
  ## =============================
  Moffitt_S2 <- list(sampInfo=sampInfo,
                     featInfo=featInfo,
                     ex=ex,
                     metadata=metadata)

  Moffitt_S2.Hs <- list(ex = ex[which(featInfo$species == "Hs"),],
             sampInfo = sampInfo,
             metadata = metadata,
             featInfo = featInfo[which(featInfo$species == "Hs"),])

  Moffitt_S2.Mm <- list(ex = ex[which(featInfo$species == "Mm"),],
             sampInfo = sampInfo,
             metadata = metadata,
             featInfo = featInfo[which(featInfo$species == "Mm"),])

  saveRDS(Moffitt_S2,
       file = "./data/Moffitt_S2.rds",
       compress = T)
  saveRDS(Moffitt_S2.Hs,
       file = "./data/Moffitt_S2.Hs.rds",
       compress = T)
  saveRDS(Moffitt_S2.Mm,
       file = "./data/Moffitt_S2.Mm.rds",
       compress = T)

  return(NULL)
}
