#' Parse data from tables for TCGA_PAAD
#'
#' @export
#' @import org.Hs.eg.db
#' @import stringr
#' @import plyr
#' @import openxlsx

parse_TCGA_PAAD <- function() {

  # library(org.Hs.eg.db)
  # library(stringr)
  # library(plyr)
  # library(openxlsx)

  ## =============================
  # Parse expression data
  ## =============================
  expression_file <- system.file("extdata/TCGA_PAAD",
                                 "PAAD.183.FullSet.fixed_20141109.txt",
                                 package = "pdacR")
  ex <- read.table(file = expression_file,header = TRUE)
  ex <- ex[-1,]
  rownames(ex) <- ex[,1]
  ex <- ex[,-1]

  featInfo <- data.frame(SYMBOL = sapply(rownames(ex),FUN = function(x) {strsplit(x,"|",fixed=TRUE)[[1]][1]}),
                         ENTREZID = sapply(rownames(ex),FUN = function(x) {strsplit(x,"|",fixed=TRUE)[[1]][2]}))

  ## =============================
  # Pull sample information
  ## =============================
  sample_table_file <- system.file("extdata/TCGA_PAAD",
                                   "Supplemental_Table_1_PAAD_sample_table_for_TCGA.xlsx",
                                   package = "pdacR")
  #---------------
  sample_table_1 <- read.xlsx(xlsxFile = sample_table_file,
                              sheet = 1)
  sample_table_1$Patient <- sapply(X = sample_table_1$Tumor.Sample.ID,
                                   FUN = function(x) {
                                     gsub(pattern = "-01A",
                                          replacement = "",
                                          x = x)})
  sample_table_1$Decision <- "whitelist"
  names(sample_table_1)
  #---------------
  sample_table_2 <- read.xlsx(xlsxFile = sample_table_file,
                              sheet = 2)
  sample_table_2$Tumor.Sample.ID <- sapply(X = sample_table_2$Patient,
                                           FUN = function(x) {
                                             paste(x,"01A",sep="-")})
  sample_table_2 <- sample_table_2[,(names(sample_table_2) %in% c("Patient","Decision","Batch","Pathology.Exclusion","Tumor.Sample.ID"))]
  names(sample_table_2) <- gsub(x=names(sample_table_2),pattern="Pathology.Exclusion",replacement = "Histological.type.by.RHH")
  summary(sample_table_2)
  #---------------
  sample_table_3 <- read.xlsx(xlsxFile = sample_table_file,
                              sheet = 3)
  sample_table_3$Tumor.Sample.ID <- sapply(X = sample_table_3$Patient,
                                           FUN = function(x) {
                                             paste(x,"01A",sep="-")})
  sample_table_3 <- sample_table_3[,(names(sample_table_3) %in% c("Patient","Decision","Batch","Notes","Tumor.Sample.ID"))]
  names(sample_table_3) <- gsub(x=names(sample_table_3),pattern="Notes",replacement = "Histological.type.by.RHH")

  summary(sample_table_3)
  #---------------
  sample_table_4 <- read.xlsx(xlsxFile = sample_table_file,
                              sheet = 4)
  sample_table_4$Tumor.Sample.ID <- sapply(X = sample_table_4$Sample.ID,
                                           FUN = function(x) {
                                             paste(x,"A",sep="")})
  sample_table_4 <- sample_table_4[,(names(sample_table_4) %in% c("Tumor.Sample.ID","Other.Notes","Pathology.Notes"))]
  sample_table_4$Decision <- "normal"
  sample_table_4$Patient <- gsub(sample_table_4$Tumor.Sample.ID,pattern="-11",replacement = "")
  names(sample_table_4) <- gsub(x=names(sample_table_4),pattern="Pathology.Notes",replacement = "Histological.type.by.RHH")
  summary(sample_table_4)

  ## =============================
  # Join sample information tables
  ## =============================

  sample_table <- join(sample_table_1,
                       sample_table_2,
                       type="full",
                       by = "Tumor.Sample.ID")

  sample_table <- join(sample_table,
                       sample_table_3,
                       type="full",
                       by = "Tumor.Sample.ID")
  sample_table <- join(sample_table,
                       sample_table_4,
                       type="full",
                       by = "Tumor.Sample.ID")
  sample_table <- join(sample_table,
                       data.frame(Tumor.Sample.ID = "TCGA-HZ-A9TJ-06A",
                                  Decision = "metastasis" ),
                       type="full",
                       by = "Tumor.Sample.ID")

  ## =============================
  # Match expression and sample data order and remove cases without metadata
  ## =============================
  data_cols_fixed <- sapply(names(ex),FUN=function(x) {substr(x = gsub(x = x,pattern = ".",replacement = "-",fixed = TRUE),
                                                              start = 1,stop = 16 )})
  col_order <- match(x = data_cols_fixed,
                     table = sample_table$Tumor.Sample.ID)
  sample_table <- sample_table[col_order,]

  no_metadata = which(is.na(sample_table$Tumor.Sample.ID))
  sample_table <- sample_table[-no_metadata,]
  ex <- ex[,-no_metadata]

  sampInfo <- sample_table

  ## =============================
  # Clean up sample information strings and adjust sample information classes
  ## =============================
  names(sampInfo) <-   gsub(pattern = "\\(|\\)|=",
                            replacement = "",
                            x = names(sampInfo),
                            fixed = FALSE)

  sampInfo <- (data.frame(lapply(X = sampInfo,FUN = function(x) {
    if(is.character(x)){
      return(as.factor(x))
    }
    return(x)
  })))

  # to numeric
  sampInfo$Follow.up.days <-
    as.numeric(as.character(sampInfo$Follow.up.days))
  sampInfo$Initial.Slide.Tumor.Cellularity <-
    as.numeric(as.character(sampInfo$Initial.Slide.Tumor.Cellularity))

  # to factor
  # sampInfo$mRNA.Bailey.Clusters.All.150.Samples.1squamous.2immunogenic.3progenitor.4ADEX <-
  #   as.factor(sampInfo$mRNA.Bailey.Clusters.All.150.Samples.1squamous.2immunogenic.3progenitor.4ADEX)
  sampInfo$Bailey.Clusters = factor(sampInfo$mRNA.Bailey.Clusters.All.150.Samples.1squamous.2immunogenic.3progenitor.4ADEX,
                                    levels = c(1,2,3,4),
                                    labels = c("Squamous","Immunogenic","Progenitor","ADEX"))
  sampInfo$mRNA.Bailey.Clusters.All.150.Samples.1squamous.2immunogenic.3progenitor.4ADEX = NULL
  # sampInfo$mRNA.Collisson.clusters.All.150.Samples.1classical.2exocrine.3QM <-
  #   as.factor(sampInfo$mRNA.Collisson.clusters.All.150.Samples.1classical.2exocrine.3QM)
  sampInfo$Collisson.Clusters = factor(sampInfo$mRNA.Collisson.clusters.All.150.Samples.1classical.2exocrine.3QM,
                                    levels = c(1,2,3),
                                    labels = c("Classical","Exocrine","QM"))
  sampInfo$mRNA.Collisson.clusters.All.150.Samples.1classical.2exocrine.3QM <- NULL
  # sampInfo$mRNA.Moffitt.clusters.All.150.Samples.1basal.2classical <-
  #   as.factor(sampInfo$mRNA.Moffitt.clusters.All.150.Samples.1basal.2classical)
  sampInfo$Moffitt.Clusters = factor(sampInfo$mRNA.Moffitt.clusters.All.150.Samples.1basal.2classical,
                                       levels = c(1,2),
                                       labels = c("Basal-like","Classical"))
  sampInfo$mRNA.Moffitt.clusters.All.150.Samples.1basal.2classical <- NULL

  sampInfo$Moffitt.HighPurity.Clusters = factor(sampInfo$mRNA.Moffitt.clusters.76.High.Purity.Samples.Only.1basal.2classical,
                                     levels = c(1,2),
                                     labels = c("Basal-like","Classical"))
  sampInfo$mRNA.Moffitt.clusters.76.High.Purity.Samples.Only.1basal.2classical <- NULL

  sampInfo$RPPA.Clusters.All.150.Samples <-
    as.factor(sampInfo$RPPA.Clusters.All.150.Samples)
  # sampInfo$KRAS.Mutated.1.or.0 <-
  #   as.factor(sampInfo$KRAS.Mutated.1.or.0)
  sampInfo$KRAS.Mutation.Status <- factor(sampInfo$KRAS.Mutated.1.or.0,
                                          levels = c(0,1),
                                          labels = c("No Mutation","Mutated"))
  sampInfo$KRAS.Mutated.1.or.0 <- NULL
  sampInfo$lncRNA.Clusters.All.150.Samples <-
    as.factor(sampInfo$lncRNA.Clusters.All.150.Samples)
  sampInfo$miRNA.Clusters.All.150.Samples <-
    as.factor(sampInfo$miRNA.Clusters.All.150.Samples)
  sampInfo$Methylation.Clusters.All.150.Samples <-
    as.factor(sampInfo$Methylation.Clusters.All.150.Samples)

  summary(sampInfo)

  ## =============================
  # Coerce expression data to numeric
  ## =============================
  ex <- apply(X=ex,MARGIN=2,FUN=as.numeric)
  dimnames(ex)[[2]] <- sampInfo$Tumor.Sample.ID

  # Prepare dataset for consensus clustering for tumor and stroma subtyping - DEPRECATED
  # ================================================
  # sampInfo$tumor.classifier.training <- FALSE
  # sampInfo$stroma.classifier.training <- FALSE
  # sampInfo[["cluster.MT"]] <- as.character(NA)
  # sampInfo[["cluster.MS"]] <- as.character(NA)
  # sampInfo$cluster.MT <- c('basal','classical')[sampInfo$mRNA.Moffitt.clusters.All.150.Samples.1basal.2classical]
  # sampleset <- which(sampInfo$Decision %in% "whitelist")
  # sampInfo$tumor.classifier.training[sampleset] <- TRUE
  # sampleset <- c(sampleset,which(sampInfo$Decision %in% "pseudonormal"))
  # sampInfo$stroma.classifier.training[sampleset] <- TRUE
  #
  # ## =============================
  # # Make stroma cluster calls
  # ## =============================
  # featureset <- which(featInfo$SYMBOL %in%
  #                       c(as.character(gene_lists$Moffitt.Activated.25),
  #                         as.character(gene_lists$Moffitt.Normal.25)))
  # sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = t(scale(t(as.matrix(log2(1+ex[featureset,sampleset]))))),
  #                                                          seed = 1234,
  #                                                          pFeature = 0.8,
  #                                                          pItem = 0.8,
  #                                                          maxK = 6,
  #                                                          reps=200,
  #                                                          distance="pearson",
  #                                                          clusterAlg="hc")[[2]]$consensusTree
  # tmp.cluster <- c("activated","normal")[cutree(tree = sampletree, k = 2)]
  # sampInfo$cluster.MS[sampleset] <- tmp.cluster

  # Convert cluster calls to factor
  # =====================================================================
  # sampInfo$cluster.MT <- as.factor(sampInfo$cluster.MT)
  # sampInfo$cluster.MS <- as.factor(sampInfo$cluster.MS)
  # sampInfo$MoffittTumor <- sampInfo$cluster.MT
  # sampInfo$MoffittStroma <- sampInfo$cluster.MS
  # sampInfo$cluster.MT <- NULL
  # sampInfo$cluster.MS <- NULL

  ## =============================
  # Wrangle survival data for shiny app
  ## =============================
  sampInfo$survivalA <- sampInfo$Follow.up.days
  sampInfo$Days.to.death <- NULL
  sampInfo$Follow.up.days <- NULL
  sampInfo$censorA.0yes.1no <- factor(sampInfo$Censored.1yes.0no,
                                      levels = c(0,1,NaN),
                                      labels = c("1","0","NA"))
  sampInfo$Censored.1yes.0no <- NULL

  ## =============================
  # Organize dataset
  ## =============================
  sampInfo$Filter_by = sampInfo$Decision
  sampInfo = sampInfo[,-which(colnames(sampInfo) == "Decision")]
  sampInfo = sampInfo[,order(colnames(sampInfo))]

  TCGA_PAAD <- list()
  TCGA_PAAD$ex <- ex
  TCGA_PAAD$sampInfo <- sampInfo
  TCGA_PAAD$metadata <- list(log.transformed = FALSE,
                             reference = "Raphael BJ, et al, Cancer Cell, 2017, PMID:28810144",
                             accession = "DbGAaP: phs000178, BroadGDAC: PAAD",
                             description = "181 macrodissected primary tumors, 150 of which are whitelisted by histologic identification of neoplastic cellularity",
                             survivalA = "overall survival days",
                             survivalB = "None",
                             exp.type = "RNAseq",
                             default_selections = list(filter_column = "Filter_by",
                                                       filter_levels = c("Filter_by:exclude",
                                                                         "Filter_by:metastasis",
                                                                         "Filter_by:normal",
                                                                         "Filter_by:pseudonormal"),
                                                       geneSets = c("Moffitt.Basal.25","Moffitt.Classical.25"),
                                                       sampleTracks = c("ABSOLUTE Purity")))

  TCGA_PAAD$featInfo <- featInfo
  TCGA_PAAD$sampInfo = TCGA_PAAD$sampInfo[,c(which(colnames(TCGA_PAAD$sampInfo) == "Tumor.Sample.ID"),
                                             which(!(colnames(TCGA_PAAD$sampInfo) == "Tumor.Sample.ID")))]


  ## =============================
  # Save dataset
  ## =============================
  saveRDS(TCGA_PAAD,
       file = "./data/TCGA_PAAD.rds",
       compress = T)
  return(NULL)
}
