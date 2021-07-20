#' Download Chen dataset from GEO
#'
#' @export
#' @import GEOquery
#' @import Biobase

pull_Chen_from_GEO <- function() {

  # library(GEOquery)
  # library(Biobase)

  ## =============================
  # Pull dataset from GEO
  ## =============================
  gset <- getGEO(filename = system.file("extdata/Chen_array_2015",
                                        "GSE57495_series_matrix.txt.gz",
                                        package = "pdac"),
                 GSEMatrix =TRUE,
                 AnnotGPL=FALSE,
                 getGPL=FALSE)

  ## =============================
  # Parse and organize sample information
  ## =============================
  samps <-gset@phenoData
  sampInfo <- samps@data
  sampInfo <- sampInfo[,c("title","geo_accession","source_name_ch1",
                          "characteristics_ch1",
                          "characteristics_ch1.1","characteristics_ch1.2")]
  names(sampInfo)[names(sampInfo)=="title"] <- "submitted_sample_id"
  sampInfo$specimen_type <- "primaryPDAC"
  sampInfo$location <- "pancreas"
  sampInfo$stage <- NA
  sampInfo$survival_months <- NA
  sampInfo$censor <- NA

  sampInfo$censor[sampInfo[,"characteristics_ch1.1"] == "vital.status: ALIVE"] <- "censor"
  sampInfo$censor[sampInfo[,"characteristics_ch1.1"] == "vital.status: DEAD"] <- "death"
  for(i in 1:length(sampInfo[[1]])){
    tmp <- as.character(sampInfo[i,"characteristics_ch1"])
    tmp <- strsplit(x=tmp,split=": ")
    sampInfo$survival_months[i] <- as.numeric(tmp[[1]][2])
  }

  for(i in 1:length(sampInfo[[1]])){
    tmp <- as.character(sampInfo[i,"characteristics_ch1.2"])
    tmp <- strsplit(x=tmp,split=": ")
    sampInfo$stage[i] <- as.character(tmp[[1]][2])
  }

  sampInfo <- sampInfo[,-which(names(sampInfo) %in% c("source_name_ch1","characteristics_ch1","characteristics_ch1.1",
                           "characteristics_ch1.2"))]

  sampInfo$specimen_type <- as.factor(sampInfo$specimen_type)
  sampInfo$location <- as.factor(sampInfo$location)
  sampInfo$stage <- as.factor(sampInfo$stage)
  sampInfo$censor <- as.factor(sampInfo$censor)

  sampInfo$survivalA <- sampInfo$survival_months * (365/12)
  sampInfo$censorA.0yes.1no <- factor(sampInfo$censor,
                                      levels = c("censor","death"),
                                      labels = c("0","1"))

  sampInfo$survival_months <- NULL
  sampInfo$censor <- NULL

  ## =============================
  # Organize expression data
  ## =============================
  ex <- exprs(gset)
  featInfo <- data.frame(PROBES = row.names(ex))
  plot(x = ex[featInfo$PROBES %in% "merck-NM_002639_at",],
       y = ex[featInfo$PROBES %in% "merck-BQ217236_a_at",])
  which(featInfo$PROBES %in% "merck-NM_002639_at")
  which(featInfo$PROBES %in% "merck-BQ217236_a_at")

  ## =============================
  # Pull matching for featInfo
  ## =============================
  gpl <- read.table(file = system.file("extdata/Chen_array_2015",
                                       "GPL15048.soft",
                                       package = "pdac"),
                    na.strings = "",
                    sep="\t",
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    skip = 1049,
                    nrows = 60607)
  neworder <- match(table = gpl$ID,x = featInfo$PROBES)
    gpl <- gpl[neworder,]

  featInfo <- data.frame(SYMBOL = gpl$GeneSymbol,
                         ENTREZID = gpl$EntrezGeneID,
                         PROBES = gpl$ID)

  plot(x = ex[featInfo$PROBES %in% "merck-NM_002639_at",],
       y = ex[featInfo$PROBES %in% "merck-BQ217236_a_at",])
  which(featInfo$PROBES %in% "merck-NM_002639_at")
  which(featInfo$PROBES %in% "merck-BQ217236_a_at")

  ## =============================
  # Aggregation by entrez and symbol
  ## =============================
  print("Dimension before aggregation")
  print(dim(ex))
  tmp <- aggregate(x = ex,
                   by = featInfo[,1:2],
                   FUN = function(x) sum(x))
  featInfo <- data.frame(SYMBOL = tmp[[1]],
                         ENTREZID = tmp[[2]])
  ex  <- tmp[,c(-1,-2)]
  print(dim(ex))

  ## =============================
  # Include dataset metadata for reference
  ## =============================
  metadata <- list(log.transformed = TRUE,
                   reference = "Chen DT et al, PLos One, 2015, PMID:26247463",
                   accession = "GEO: GSE57495",
                   description = "63 fresh frozen macrodissected PDAC tumor samples from Moffitt Cancer Center, microarray",
                   survivalA = "overall survival days",
                   survivalB = "None")

  ## =============================
  # Perform consensus clustering for tumor subtyping
  ## =============================
  sampInfo$tumor.classifier.training <- FALSE
  sampInfo$stroma.classifier.training <- FALSE
  sampInfo[["MoffittTumor"]] <- as.character(NA)
  sampInfo[["MoffittStroma"]] <- as.character(NA)

  #############

  sampleset <- which(sampInfo$specimen_type %in% "primaryPDAC")
  featureset <- which(featInfo$SYMBOL %in%
                        c(as.character(pdac::gene_lists$Moffitt.Basal.25),
                          as.character(pdac::gene_lists$Moffitt.Classical.25)))

  smallx <- t(scale(scale=FALSE,
                    center=TRUE,
                    x=t((ex[featureset,sampleset]))))
  sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = smallx,
                                                           seed = 1234,
                                                           pFeature = 0.8,
                                                           pItem = 0.8,
                                                           maxK = 6,
                                                           reps=200,
                                                           distance="pearson",
                                                           clusterAlg="kmdist")[[2]]$consensusTree
  tmp.cluster <- c("basal","classical")[cutree(tree = sampletree, k = 2)]
  sampInfo$MoffittTumor[sampleset] <- tmp.cluster
  ggplot(data = data.frame(classical = (colMeans(ex[featInfo$SYMBOL %in%
                                                      pdac::gene_lists$Moffitt.Classical.25,sampleset])),
                           basal = (colMeans(ex[featInfo$SYMBOL %in%
                                                  pdac::gene_lists$Moffitt.Basal.25,sampleset])),
                           MoffittTumor = sampInfo$MoffittTumor[sampleset]),
         aes(x = basal,y = classical, color = MoffittTumor)) +
    geom_point(size=2)




  featureset <- which(featInfo$SYMBOL %in%
                        c(as.character(pdac::gene_lists$Moffitt.Normal.25),
                          as.character(pdac::gene_lists$Moffitt.Activated.25)))

  smallx <- t(scale(scale=FALSE,
                    center=TRUE,
                    x=t((ex[featureset,sampleset]))))
  sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = smallx,
                                                           seed = 1234,
                                                           pFeature = 0.8,
                                                           pItem = 0.8,
                                                           maxK = 6,
                                                           reps=200,
                                                           distance="euclidean",
                                                           clusterAlg="kmdist")[[2]]$consensusTree
  tmp.cluster <- c("activated","normal")[cutree(tree = sampletree, k = 2)]
  sampInfo$MoffittStroma[sampleset] <- tmp.cluster
  ggplot(data = data.frame(normal = (colMeans(ex[featInfo$SYMBOL %in%
                                                      pdac::gene_lists$Moffitt.Normal.25,sampleset])),
                           activated = (colMeans(ex[featInfo$SYMBOL %in%
                                                  pdac::gene_lists$Moffitt.Activated.25,sampleset])),
                           MoffittStroma = sampInfo$MoffittStroma[sampleset]),
         aes(x = normal,y = activated, color = MoffittStroma)) +
    geom_point(size=2)
  #############

  ## =============================
  # Alter variable types for GUI use
  ## =============================

  cols_to_change <- c("geo_accession",
                      "MoffittTumor",
                      "MoffittStroma")

  for(col in cols_to_change){
    sampInfo[[col]] <- as.factor(sampInfo[[col]])
  }

  ## =============================
  # Save dataset
  ## =============================
  Chen_GEO_array <- list(sampInfo=sampInfo,
                         featInfo=featInfo,
                         ex=ex,
                         metadata=metadata)

  save(list = c("Chen_GEO_array"),
       file = "./data/Chen_GEO_array.RData")

  return(NULL)
}
