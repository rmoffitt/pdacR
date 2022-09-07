#' Parse data from tables
#'
#' @export
#' @import reshape2
#' @import org.Hs.eg.db
#' @import illuminaHumanv4.db
#' @import plyr


parse_ICGC <- function() {

  # ------------------------------------------------------------------------------------
  # supplemental table #1 from Bailey et al, Nature, 2016
  Bailey_sup1_fname <- system.file("extdata/ICGC",
                                   "Bailey_sup_table_1.txt",
                                   package = "pdacR")
  Bailey_sup1 <- read.table(file = Bailey_sup1_fname,
                            sep = "\t",
                            header = TRUE,
                            skip = 1)
  # by donor  : Bailey_sup1$icgc_id ........................... (ICGC_0245)
  rm("Bailey_sup1_fname")
  Bailey_sup1$qpure_score <- as.numeric(as.character(Bailey_sup1$qpure_score))
  Bailey_sup1$Gender <- as.character(Bailey_sup1$Gender)
  Bailey_sup1$Gender[Bailey_sup1$Gender %in% c("N/A","Not documented")] <- NA
  Bailey_sup1$Gender[Bailey_sup1$Gender %in% c("Female")] <- "female"
  Bailey_sup1$Gender[Bailey_sup1$Gender %in% c("Male")] <- "male"
  Bailey_sup1$Gender <- factor(Bailey_sup1$Gender)
  Bailey_sup1$Age.at.Diagnosis.in.Years <- as.numeric(as.character(Bailey_sup1$Age.at.Diagnosis.in.Years))
  Bailey_sup1$survival_months <- as.numeric(as.character(Bailey_sup1$Length.of.Follow.Up..months.))
  Bailey_sup1$censored <- as.character(Bailey_sup1$Status)
  Bailey_sup1$censored[Bailey_sup1$censored %in% c("Unknown",
                                                   "Not documented")] <- NA
  Bailey_sup1$censored[Bailey_sup1$censored %in% c("Alive - Without Disease",
                                                   "Alive - With Disease",
                                                   "Deceased - Of Other Cause",
                                                   "Lost to Follow Up",
                                                   "Alive - Disease Status Unknown")] <- "censor"
  Bailey_sup1$censored[Bailey_sup1$censored %in% c("Deceased - Of Disease",
                                                   "Deceased - Of Unknown Cause")] <- "death"
  Bailey_sup1$censored <- as.factor(Bailey_sup1$censored)
  # ------------------------------------------------------------------------------------
  # supplemental table #2 from Bailey et al, Nature, 2016
  Bailey_sup2_fname <- system.file("extdata/ICGC",
                                   "Bailey_sup_table_2.csv",
                                   package = "pdacR")
  Bailey_sup2 <- read.table(file = Bailey_sup2_fname,
                            sep = ",",
                            header = TRUE,
                            skip = 0)
  # by donor  : Bailey_sup2$icgc_id ........................... (ICGC_0006)
  rm("Bailey_sup2_fname")

  # ------------------------------------------------------------------------------------
  # PACA-CA sample list downloaded from ICGC data portal
  PACA_CA_samplesheet_fname <- system.file("extdata/ICGC",
                                           "specimen.tsv",
                                           package = "pdacR")
  PACA_CA_samplesheet <- read.table(file = PACA_CA_samplesheet_fname,
                                    sep = "\t",
                                    header = TRUE)
  # PACA_CA_samplesheet_seq <- droplevels(PACA_CA_samplesheet[PACA_CA_samplesheet$sequencing_strategy %in% "RNA-Seq",
  #                                                           -which(names(PACA_CA_samplesheet) %in% "raw_data_accession")])
  PACA_CA_samplesheet_seq <- PACA_CA_samplesheet[!duplicated(PACA_CA_samplesheet),]
  # by sample : PACA_CA_samplesheet_seq$submitted_sample_id ... (MPCC_0001_Pa_C)
  # by sample : PACA_CA_samplesheet_seq$icgc_sample_id ........ (SA412478)
  # by donor  : PACA_CA_samplesheet_seq$submitted_donor_id .... (PCSI_0022)

  # PACA_CA_samplesheet_array <- droplevels(PACA_CA_samplesheet[PACA_CA_samplesheet$sequencing_strategy %in% "non-NGS",
  #                                                           -which(names(PACA_CA_samplesheet) %in% "raw_data_accession")])
  # PACA_CA_samplesheet_array <- PACA_CA_samplesheet_array[!duplicated(PACA_CA_samplesheet_array),]
  # by sample : PACA_CA_samplesheet_array$submitted_sample_id . (PCSI_0001_Pa_P)
  # by sample : PACA_CA_samplesheet_array$icgc_sample_id ...... (SA412485)
  # by donor  : PACA_CA_samplesheet_array$submitted_donor_id .. (PCSI_0002)

  rm("PACA_CA_samplesheet_fname")

  # ------------------------------------------------------------------------------------
  # PACA-CA exp_seq downloaded from ICGC data portal
  ## File is too large for github without LFS. Parse code here for transparency but file is not present.
  ## If you want to replicate parsing, download the relevant files from https://dcc.icgc.org/search?filters=%7B%22donor%22:%7B%22availableDataTypes%22:%7B%22is%22:%5B%22exp_seq%22%5D%7D,%22projectId%22:%7B%22is%22:%5B%22PACA-CA%22%5D%7D%7D%7D
  ## Sequencing Based Gene Expression ^
  PACA_CA_exp_seq_fname <- system.file("extdata/ICGC",
                                       "PACA_CA_Seq_Updated.tsv.gz",
                                       package = "pdacR")
  PACA_CA_exp_seq <- read.table(gzfile(PACA_CA_exp_seq_fname ),
                                sep = "\t",
                                header = TRUE)
  PACA_CA_seq <- list()
  PACA_CA_seq$ex <- reshape2::dcast(data = PACA_CA_exp_seq,
                                    formula = gene_id ~ icgc_specimen_id,
                                    value.var = "normalized_read_count",
                                    fun.aggregate = sum)
  rownames(PACA_CA_seq$ex ) <- PACA_CA_seq$ex$gene_id
  PACA_CA_seq$featInfo <- data.frame(ENSEMBL = PACA_CA_seq$ex$gene_id)
  PACA_CA_seq$ex  <- PACA_CA_seq$ex [,-which(names(PACA_CA_seq$ex ) %in% "gene_id")]

  PACA_CA_seq$featInfo$ENTREZID <- as.factor(AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                                                   keys = as.character(PACA_CA_seq$featInfo$ENSEMBL),
                                                                   keytype = "ENSEMBL",
                                                                   column = "ENTREZID",
                                                                   multiVals = "first"))
  PACA_CA_seq$featInfo$SYMBOL <- as.factor(AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                                                 keys = as.character(PACA_CA_seq$featInfo$ENSEMBL),
                                                                 keytype = "ENSEMBL",
                                                                 column = "SYMBOL",
                                                                 multiVals = "first"))

  PACA_CA_seq$sampInfo <- data.frame(icgc_specimen_id = names(PACA_CA_seq$ex))
  PACA_CA_seq$sampInfo <- cbind(PACA_CA_seq$sampInfo,
                                PACA_CA_samplesheet_seq[match(PACA_CA_seq$sampInfo$icgc_specimen_id,
                                                              PACA_CA_samplesheet_seq$icgc_specimen_id),
                                                        c("submitted_specimen_id",
                                                          "icgc_donor_id",
                                                          "submitted_donor_id",
                                                          "specimen_type",
                                                          "specimen_type_other")])

  PACA_CA_seq$sampInfo$dataset <- "PACA_CA_seq"

  rm("PACA_CA_samplesheet_seq")
  rm("PACA_CA_exp_seq")
  rm("PACA_CA_exp_seq_fname")

  ## Join in clinical information
  clinical_info_fname = system.file("extdata/ICGC",
                                    "donor.tsv",
                                    package = "pdacR")
  clinical_info <- read.table(clinical_info_fname,
                              sep = "\t", header = T)

  PACA_CA_seq$sampInfo = dplyr::full_join(PACA_CA_seq$sampInfo, clinical_info, by = "icgc_donor_id")
  PACA_CA_seq$sampInfo$survivalA = as.integer(PACA_CA_seq$sampInfo$donor_survival_time)
  PACA_CA_seq$sampInfo$censorA.0yes.1no = factor(ifelse(PACA_CA_seq$sampInfo$donor_vital_status == "alive",0,1))
  PACA_CA_seq$sampInfo$censorA.0yes.1no[which(PACA_CA_seq$sampInfo$donor_vital_status == "")] = NA

  PACA_CA_seq$sampInfo$survivalB = as.integer(PACA_CA_seq$sampInfo$donor_relapse_interval)
  PACA_CA_seq$sampInfo$censorB.0yes.1no = factor(ifelse(PACA_CA_seq$sampInfo$donor_relapse_type == "",0,1))
  PACA_CA_seq$sampInfo$censorB.0yes.1no[which(PACA_CA_seq$sampInfo$disease_status_last_followup == "")] = NA

  ## Identify primary sample per patient for survival filtering
  dupes = PACA_CA_seq$sampInfo$icgc_donor_id[which(duplicated(PACA_CA_seq$sampInfo$icgc_donor_id))]
  PACA_CA_seq$sampInfo$survival_patient_list = factor(ifelse(PACA_CA_seq$sampInfo$icgc_donor_id %in% dupes &
                                                               PACA_CA_seq$sampInfo$specimen_type != "Primary tumour - other", "drop","keep"))

  ### Deprecated with new release
  # # ------------------------------------------------------------------------------------
  # # PACA-CA from ICGC data portal has incorrect specimen codes
  # head(PACA_CA_seq$sampInfo)
  # PACA_CA_seq$sampInfo$actual_type <- factor(
  #   x = sapply(X = PACA_CA_seq$sampInfo$submitted_sample_id,
  #              FUN = function(x){
  #                x <- as.character(x)
  #                y <- strsplit(x = x, split = "_")[[1]]
  #                z <- y[length(y)]
  #                return(z)
  #              }),
  #   levels = c("C","X","P"),
  #   labels = c("Cell line","Xenograft","Primary tumour")
  # )
  # summary(PACA_CA_seq$sampInfo)
  #
  #
  # # ------------------------------------------------------------------------------------
  # # PACA-CA exp_array downloaded from ICGC data portal returns no results!
  # rm("PACA_CA_samplesheet_array")


  # ------------------------------------------------------------------------------------
  # PACA-AU sample list downloaded from ICGC data portal
  PACA_AU_samplesheet_fname <- system.file("extdata/ICGC",
                                           "sample.PACA-AU.1478117594130.tsv",
                                           package = "pdacR")
  PACA_AU_samplesheet <- read.table(file = PACA_AU_samplesheet_fname,
                                    sep = "\t",
                                    header = TRUE)
  PACA_AU_samplesheet_seq <- droplevels(PACA_AU_samplesheet[PACA_AU_samplesheet$sequencing_strategy %in% "RNA-Seq",
                                                            -which(names(PACA_AU_samplesheet) %in% "raw_data_accession")])
  PACA_AU_samplesheet_seq <- PACA_AU_samplesheet_seq[!duplicated(PACA_AU_samplesheet_seq),]
  # by sample : PACA_AU_samplesheet_seq$submitted_sample_id ... (8013915)
  # by sample : PACA_AU_samplesheet_seq$icgc_sample_id ........ (SA412478)
  # by donor  : PACA_AU_samplesheet_seq$submitted_donor_id .... (ICGC_0004)

  PACA_AU_samplesheet_array <- droplevels(PACA_AU_samplesheet[PACA_AU_samplesheet$sequencing_strategy %in% "non-NGS",
                                                              -which(names(PACA_AU_samplesheet) %in% "raw_data_accession")])
  PACA_AU_samplesheet_array <- PACA_AU_samplesheet_array[!duplicated(PACA_AU_samplesheet_array),]
  # by sample : PACA_AU_samplesheet_array$submitted_sample_id . (8031121)
  # by sample : PACA_AU_samplesheet_array$icgc_sample_id ...... (SA407710)
  # by donor  : PACA_AU_samplesheet_array$submitted_donor_id .. (ICGC_0327)

  rm("PACA_AU_samplesheet_fname")

  # ------------------------------------------------------------------------------------
  # PACA-AU exp_seq downloaded from ICGC data portal
  PACA_AU_exp_seq_fname <- system.file("extdata/ICGC",
                                       "PACA-AU.exp_seq.tsv.gz",
                                       package = "pdacR")
  PACA_AU_exp_seq <- read.table(gzfile(PACA_AU_exp_seq_fname ),
                                sep = "\t",
                                header = TRUE)

  PACA_AU_seq <- list()
  PACA_AU_seq$ex <- reshape2::dcast(data = PACA_AU_exp_seq,
                                    formula = gene_id ~ icgc_sample_id,
                                    value.var = "normalized_read_count",
                                    fun.aggregate = sum)
  rownames(PACA_AU_seq$ex ) <- PACA_AU_seq$ex$gene_id
  PACA_AU_seq$featInfo <- data.frame(ENSEMBL = PACA_AU_seq$ex$gene_id)
  PACA_AU_seq$ex  <- PACA_AU_seq$ex [,-which(names(PACA_AU_seq$ex ) %in% "gene_id")]

  PACA_AU_seq$featInfo$ENTREZID <- as.factor(AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                                                   keys = as.character(PACA_AU_seq$featInfo$ENSEMBL),
                                                                   keytype = "ENSEMBL",
                                                                   column = "ENTREZID",
                                                                   multiVals = "first"))
  PACA_AU_seq$featInfo$SYMBOL <- as.factor(AnnotationDbi::mapIds(x = org.Hs.eg.db,
                                                                 keys = as.character(PACA_AU_seq$featInfo$ENSEMBL),
                                                                 keytype = "ENSEMBL",
                                                                 column = "SYMBOL",
                                                                 multiVals = "first"))

  PACA_AU_seq$sampInfo <- data.frame(icgc_sample_id = names(PACA_AU_seq$ex))
  PACA_AU_seq$sampInfo <- cbind(PACA_AU_seq$sampInfo,
                                PACA_AU_samplesheet_seq[match(PACA_AU_seq$sampInfo$icgc_sample_id,
                                                              PACA_AU_samplesheet_seq$icgc_sample_id),
                                                        c("submitted_sample_id",
                                                          "icgc_specimen_id",
                                                          "submitted_specimen_id",
                                                          "icgc_donor_id",
                                                          "submitted_donor_id",
                                                          "specimen_type")])

  PACA_AU_seq$sampInfo$dataset <- "PACA_AU_seq"

  rm("PACA_AU_samplesheet_seq")
  rm("PACA_AU_exp_seq")
  rm("PACA_AU_exp_seq_fname")

  # ------------------------------------------------------------------------------------
  # PACA-AU exp_array downloaded from ICGC data portal
  PACA_AU_exp_array_fname <- system.file("extdata/ICGC",
                                         "PACA-AU.exp_array.tsv.gz",
                                         package = "pdacR")
  PACA_AU_exp_array <- read.table(gzfile(PACA_AU_exp_array_fname ),
                                  sep = "\t",
                                  header = TRUE)

  PACA_AU_array <- list()
  PACA_AU_array$ex <- reshape2::dcast(data = PACA_AU_exp_array,
                                      formula = gene_id ~ icgc_sample_id,
                                      value.var = "normalized_expression_value",
                                      fun.aggregate = sum)
  rownames(PACA_AU_array$ex ) <- PACA_AU_array$ex$gene_id
  PACA_AU_array$featInfo <- data.frame(PROBEID = PACA_AU_array$ex$gene_id)
  PACA_AU_array$ex  <- PACA_AU_array$ex [,-which(names(PACA_AU_array$ex ) %in% "gene_id")]

  PACA_AU_array$featInfo$ENTREZID <- as.factor(AnnotationDbi::mapIds(x = illuminaHumanv4.db,
                                                                     keys = as.character(PACA_AU_array$featInfo$PROBEID),
                                                                     keytype = "PROBEID",
                                                                     column = "ENTREZID",
                                                                     multiVals = "first"))
  PACA_AU_array$featInfo$ENSEMBL <- as.factor(AnnotationDbi::mapIds(x = illuminaHumanv4.db,
                                                                    keys = as.character(PACA_AU_array$featInfo$PROBEID),
                                                                    keytype = "PROBEID",
                                                                    column = "ENSEMBL",
                                                                    multiVals = "first"))
  PACA_AU_array$featInfo$SYMBOL <- as.factor(AnnotationDbi::mapIds(x = illuminaHumanv4.db,
                                                                   keys = as.character(PACA_AU_array$featInfo$PROBEID),
                                                                   keytype = "PROBEID",
                                                                   column = "SYMBOL",
                                                                   multiVals = "first"))

  PACA_AU_array$sampInfo <- data.frame(icgc_sample_id = names(PACA_AU_array$ex))
  PACA_AU_array$sampInfo <- cbind(PACA_AU_array$sampInfo,
                                  PACA_AU_samplesheet[match(PACA_AU_array$sampInfo$icgc_sample_id,
                                                            PACA_AU_samplesheet$icgc_sample_id),
                                                      c("submitted_sample_id",
                                                        "icgc_specimen_id",
                                                        "submitted_specimen_id",
                                                        "icgc_donor_id",
                                                        "submitted_donor_id",
                                                        "specimen_type")])

  PACA_AU_array$sampInfo$dataset <- "PACA_AU_array"

  rm("PACA_AU_samplesheet_array")
  rm("PACA_AU_exp_array")
  rm("PACA_AU_exp_array_fname")
  # ------------------------------------------------------------------------------------

  names(Bailey_sup1)[names(Bailey_sup1)=="icgc_id"] <- "submitted_donor_id"
  names(Bailey_sup2)[names(Bailey_sup2)=="icgc_id"] <- "submitted_donor_id"

  PACA_AU_array$sampInfo <- join(x = PACA_AU_array$sampInfo,
                                 by = "submitted_donor_id",
                                 y = Bailey_sup1[,!(names(Bailey_sup1) %in%
                                                      c("EXOME","WGS","SNP","RNASeq","HT.12.Expression.array",
                                                        "X450K.Methylation.arrays","SV.analysis"))])
  PACA_AU_seq$sampInfo <- join(x = PACA_AU_seq$sampInfo,
                               by = "submitted_donor_id",
                               y = Bailey_sup1[,!(names(Bailey_sup1) %in%
                                                    c("EXOME","WGS","SNP","RNASeq","HT.12.Expression.array",
                                                      "X450K.Methylation.arrays","SV.analysis"))])

  PACA_AU_array$sampInfo <- join(x = PACA_AU_array$sampInfo,
                                 by = "submitted_donor_id",
                                 y = Bailey_sup2)
  PACA_AU_seq$sampInfo <- join(x = PACA_AU_seq$sampInfo,
                               by = "submitted_donor_id",
                               y = Bailey_sup2)

  # Adjust survival data for GUI
  ## =============================

  # PACA AU array
  # -------------------------------
  PACA_AU_array$sampInfo$survivalA <- "NA"
  PACA_AU_array$sampInfo$survivalA[which(!PACA_AU_array$sampInfo$Length.of.Follow.Up..months. %in% c("Not documented", "NA"))] <-
    as.character(PACA_AU_array$sampInfo$Length.of.Follow.Up..months.[which(!PACA_AU_array$sampInfo$Length.of.Follow.Up..months. %in% c("Not documented", "NA"))])
  PACA_AU_array$sampInfo$survivalA <- as.numeric(PACA_AU_array$sampInfo$survivalA) * (365/12)

  PACA_AU_array$sampInfo$censorA.0yes.1no <- as.character(PACA_AU_array$sampInfo$censored)
  PACA_AU_array$sampInfo$censorA.0yes.1no[which(PACA_AU_array$sampInfo$Length.of.Follow.Up..months. %in% c("Not documented", "NA"))] <- "NA"
  PACA_AU_array$sampInfo$censorA.0yes.1no[which(PACA_AU_array$sampInfo$censorA.0yes.1no %in% "censor")] <- 0
  PACA_AU_array$sampInfo$censorA.0yes.1no[which(PACA_AU_array$sampInfo$censorA.0yes.1no %in% "death")] <- 1
  PACA_AU_array$sampInfo$censorA.0yes.1no <- as.numeric(PACA_AU_array$sampInfo$censorA.0yes.1no)

  PACA_AU_array$sampInfo$survival_months <- NULL
  PACA_AU_array$sampInfo$censored <- NULL
  PACA_AU_array$sampInfo$Length.of.Follow.Up..months. <- NULL

  # PACA AU seq
  # -------------------------------
  PACA_AU_seq$sampInfo$survivalA <- "NA"
  PACA_AU_seq$sampInfo$survivalA[which(!PACA_AU_seq$sampInfo$Length.of.Follow.Up..months. %in% c("Not documented", "NA"))] <-
    as.character(PACA_AU_seq$sampInfo$Length.of.Follow.Up..months.[which(!PACA_AU_seq$sampInfo$Length.of.Follow.Up..months. %in% c("Not documented", "NA"))])
  PACA_AU_seq$sampInfo$survivalA <- as.numeric(PACA_AU_seq$sampInfo$survivalA) * (365/12)

  # change this so that 1 = no (not censored) and 0 = yes (censored)
  PACA_AU_seq$sampInfo$censorA.0yes.1no <- as.character(PACA_AU_seq$sampInfo$censored)
  PACA_AU_seq$sampInfo$censorA.0yes.1no[which(PACA_AU_seq$sampInfo$Length.of.Follow.Up..months. %in% c("Not documented", "NA"))] <- "NA"
  PACA_AU_seq$sampInfo$censorA.0yes.1no[which(PACA_AU_seq$sampInfo$censorA.0yes.1no %in% "censor")] <- 0
  PACA_AU_seq$sampInfo$censorA.0yes.1no[which(PACA_AU_seq$sampInfo$censorA.0yes.1no %in% "death")] <- 1
  PACA_AU_seq$sampInfo$censorA.0yes.1no <- as.numeric(PACA_AU_seq$sampInfo$censorA.0yes.1no)

  PACA_AU_seq$sampInfo$survival_months <- NULL
  PACA_AU_seq$sampInfo$censored <- NULL
  PACA_AU_seq$sampInfo$Length.of.Follow.Up..months. <- NULL


  ## =============================
  PACA_AU_array$metadata <- list(log.transformed = TRUE,
                                 reference = "Bailey P et al, Nature, 2016, PMID:26909576 ",
                                 accession = c("EGAS00001000154",
                                               "GSE36924",
                                               "GSE49149",
                                               "PACA-AU"),
                                 description = "131 primary pancreatic tumors: 101 PDAC, 30 other",
                                 survivalA = "overall survival days",
                                 survivalB = "None",
                                 exp.type = "Array",
                                 default_selections = list(filter_column = "HistoSubtype",
                                                           filter_levels = levels(as.factor(PACA_AU_array$HistoSubtype))[-5]))

  PACA_AU_seq$metadata <- list(log.transformed = TRUE,
                               reference = "Bailey P et al, Nature, 2016, PMID:26909576 ",
                               accession = c("EGAS00001000154",
                                             "GSE36924",
                                             "GSE49149",
                                             "PACA-AU"),
                               description = "82 primary pancreatic tumors, 8 cell lines, and 2 metastatic tumors",
                               survivalA = "overall survival days",
                               survivalB = "None",
                               exp.type = "RNAseq",
                               default_selections = list(filter_column = "HistoSubtype",
                                                         filter_levels = levels(as.factor(PACA_AU_seq$HistoSubtype))[-5]))

  PACA_CA_seq$metadata <- list(log.transformed = FALSE,
                               reference = "ICGC Data Portal",
                               accession = "PACA-CA",
                               description = "40 cell lines, 11 xenografts, 16 mets, and 195 pancreatic primary tumor",
                               survivalA = "overall survival days FILTER BY SURVIVAL_PATIENT_LIST",
                               survivalB = "progression or relapse interval",
                               exp.type = "RNAseq",
                               default_selections = list(sampleTracks = "specimen_type"))


  ## =============================
  # Consensus clustering to find tumor and stroma subtypes
  ## =============================

  # PACA_AU_array
  # ----------------------------------------------------
  # Tumor
  # ----------------------
  # dataset <- PACA_AU_array
  # sampleset <- 1:nrow(dataset$sampInfo)
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
  # tmp.cluster <- c("classical","basal")[cutree(tree = sampletree, k = 2)]
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

  # =======================
  #PACA_AU_array <- dataset
  # =======================

  # PACA_AU_seq
  # ----------------------------------------------------
  # Tumor
  # ----------------------
  # dataset <- PACA_AU_seq
  # sampleset <- which((dataset$sampInfo$HistoSubtype %in% "Pancreatic Ductal Adenocarcinoma")
  #                    & (dataset$sampInfo$Sample.type %in% "Primary tumour"))
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
  # tmp.cluster <- c("classical","basal")[cutree(tree = sampletree, k = 2)]
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
  # tmp.cluster <- c("activated","normal")[cutree(tree = sampletree, k = 2)]
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
  #
  # # =======================
  # PACA_AU_seq <- dataset
  # # =======================

  # PACA_CA_seq
  # ----------------------------------------------------
  # Tumor
  # ----------------------
  # dataset <- PACA_CA_seq
  # sampleset <- which(dataset$sampInfo$actual_type %in% "Cell line")
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
  # tmp.cluster <- c("activated","normal")[cutree(tree = sampletree, k = 2)]
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

  # =======================
  # PACA_CA_seq <- dataset
  # =======================

  ### Reorganize sampInfo for ease of interepretation
  PACA_CA_seq$sampInfo = PACA_CA_seq$sampInfo[,order(colnames(PACA_CA_seq$sampInfo))]
  PACA_CA_seq$sampInfo = PACA_CA_seq$sampInfo[,c(which(colnames(PACA_CA_seq$sampInfo) == "icgc_specimen_id"),
                                                 which(!(colnames(PACA_CA_seq$sampInfo) == "icgc_specimen_id")))]

  PACA_AU_seq$sampInfo = PACA_AU_seq$sampInfo[,order(colnames(PACA_AU_seq$sampInfo))]
  PACA_AU_seq$sampInfo = PACA_AU_seq$sampInfo[,c(which(colnames(PACA_AU_seq$sampInfo) == "icgc_sample_id"),
                                                 which(!(colnames(PACA_AU_seq$sampInfo) == "icgc_sample_id")))]

  PACA_AU_array$sampInfo = PACA_AU_array$sampInfo[,order(colnames(PACA_AU_array$sampInfo))]
  PACA_AU_array$sampInfo = PACA_AU_array$sampInfo[,c(which(colnames(PACA_AU_array$sampInfo) == "icgc_sample_id"),
                                                     which(!(colnames(PACA_AU_array$sampInfo) == "icgc_sample_id")))]

  # ------------------------------------------------------------------------------------
  # save(list = c("Bailey_sup1","Bailey_sup2","PACA_AU_samplesheet","PACA_CA_samplesheet"),
  #      file = "./data/ICGC_metadata.RData",
  #      compress = T)

  saveRDS(PACA_AU_array,
          file = "./data/PACA_AU_array.rds",
          compress = T)
  saveRDS(PACA_AU_seq,
          file = "./data/PACA_AU_seq.rds",
          compress = T)
  saveRDS(PACA_CA_seq,
          file = "./data/PACA_CA_seq.rds",
          compress = T)
  return(NULL)
}
