#' Parse data from tables for Puleo et al, 2018
#'
#' @export
#' @import org.Hs.eg.db
#' @import stringr
#' @import plyr
#' @import openxlsx
#' @import hgu219.db

parse_Puleo_Array <- function() {

  # library(org.Hs.eg.db)
  # library(stringr)
  # library(plyr)
  # library(openxlsx)
  # library(hgu219.db)

  ## =============================
  # Read expression file and merge expression to master data frame
  ## =============================

  expression_file1 <- system.file("extdata/Puleo",
                                 "ProcessedExpression1.txt",
                                 package = "pdac")

  expression_file2 <- system.file("extdata/Puleo",
                                 "ProcessedExpression2.txt",
                                 package = "pdac")

  expression_file3 <- system.file("extdata/Puleo",
                                 "ProcessedExpression3.txt",
                                 package = "pdac")

  expression_file4 <- system.file("extdata/Puleo",
                                 "ProcessedExpression4.txt",
                                 package = "pdac")

  ex1 <- read.table(expression_file1,
                   sep = "\t",
                   header = TRUE)

  ex2 <- read.table(expression_file2,
                   sep = "\t",
                   header = TRUE)

  ex3 <- read.table(expression_file3,
                   sep = "\t",
                   header = TRUE)

  ex4 <- read.table(expression_file4,
                   sep = "\t",
                   header = TRUE)

  ex <- cbind(ex1,
              ex2,
              ex3,
              ex4)

  s_at <- grep(pattern = "s",
               row.names(ex))

  x_at <- grep(pattern = "x",
               row.names(ex))

  ex_clean <- ex[-c(s_at,x_at),]

  print(dim(ex_clean))

  ## =============================
  # Choose best case PROBEID-SYMBOL matches
  ## =============================

  df<- select(hgu219.db,
              keys = row.names(ex_clean),
              column = c("PROBEID","SYMBOL", "ENTREZID"))

  df <- df[-which(duplicated(df$PROBEID)),]

  tmp <- aggregate(x = ex_clean,
                   by = df[,2:3],
                   FUN = function(x) sum(x))


  featInfo <- data.frame(SYMBOL = tmp[[1]],
                         ENTREZID = tmp[[2]])

  ex <- tmp[,-c(1,2)]

  print(dim(ex))

  ## =============================
  # Read sample information and ensure matching with expression file
  ## =============================

  sample_file <- system.file("extdata/Puleo",
                             "E-MTAB-6134.sdrf.txt",
                             package = "pdac")


  sample_info <- read.table(file = sample_file,
                            sep = "\t",
                            header = TRUE)



  print(summary(sample_info))

  ## =============================
  # Remove unwanted sample info factors and rename factors
  ## =============================
  print(dim(sample_info))

  dropfactors <- which(names(sample_info) %in% c("Characteristics.organism.",
                                                 "Characteristics.developmental.stage.",
                                                 "Characteristics.organism.part.",
                                                 "Characteristics.disease.",
                                                 "Characteristics.sampling.site.",
                                                 "Term.Source.REF",
                                                 "Term.Accession.Number",
                                                 "Term.Source.REF.1",
                                                 "Term.Accession.Number.1",
                                                 "Term.Source.REF.2",
                                                 "Term.Accession.Number.2",
                                                 "Material.Type",
                                                 "Description",
                                                 "Protocol.REF",
                                                 "Protocol.REF.1",
                                                 "Protocol.REF.2",
                                                 "Labeled.Extract.Name",
                                                 "Label",
                                                 "Protocol.REF.3",
                                                 "Protocol.REF.4",
                                                 "Assay.Name",
                                                 "Technology.Type",
                                                 "Array.Design.REF",
                                                 "Term.Source.REF.3",
                                                 "Comment..ArrayExpress.FTP.file.",
                                                 "Protocol.REF.5",
                                                 "Derived.Array.Data.File",
                                                 "Comment..Derived.ArrayExpress.FTP.file.",
                                                 "Extract.Name",
                                                 "Factor.Value.wholetumclassif."))

  print(length(dropfactors))

  sampInfo <- sample_info[, -dropfactors]

  print(dim(sampInfo))

  names(sampInfo) <- c("Sample.name",
                       "Sex",
                       "Tumor.grade",
                       "TNM.tumor.grade",
                       "Resection.margin",
                       "Clinical.center",
                       "FFPEblock.age.years",
                       "FFPEtime",
                       "OS.delay.months",
                       "ostime",
                       "OS.censor.0yes.1no",
                       "DFS.delay.months",
                       "DFStime",
                       "DFS.censor.0yes.1no",
                       "HighTumorCellClass",
                       "WholeTumorClass",
                       "Average.VAF",
                       "KRAS.mut.1yes.0no",
                       "p53.mut.1yes.0no",
                       "CDKN2a.mut.1yes.0no",
                       "Array.file")

  sampInfo$Tumor.grade <- factor(sampInfo$Tumor.grade,
                                 levels = c("not available",
                                            "tumour grading G1",
                                            "tumour grading G2",
                                            "tumour grading G3"),
                                 labels = c("NA", "G1", "G2", "G3"))

  sampInfo$Resection.margin <- factor(sampInfo$Resection.margin,
                                      levels = c("not available",
                                                 "resection margin R0",
                                                 "resection margin R1"),
                                      labels = c("NA", "R0", "R1"))

  sampInfo$OS.delay.months <- as.character(sampInfo$OS.delay.months)
  sampInfo$OS.delay.months[which(sampInfo$OS.delay.months %in% "not available")] <- "NA"

  sampInfo$OS.censor.0yes.1no <- factor(sampInfo$OS.censor.0yes.1no,
                                      levels = c("not available",
                                                 "1",
                                                 "0"),
                                      labels = c("NA", "0", "1"))

  sampInfo$DFS.delay.months <- as.character(sampInfo$DFS.delay.months)
  sampInfo$DFS.delay.months[which(sampInfo$DFS.delay.months %in% "not available")] <- "NA"

  sampInfo$DFS.censor.0yes.1no <- factor(sampInfo$DFS.censor.0yes.1no,
                                        levels = c("not available",
                                                   "1",
                                                   "0"),
                                        labels = c("NA", "0", "1"))

  sampInfo$HighTumorCellClass <- as.character(sampInfo$HighTumorCellClass)
  sampInfo$HighTumorCellClass[which(sampInfo$HighTumorCellClass %in% "not available")] <- "NA"
  sampInfo$HighTumorCellClass <- as.factor(sampInfo$HighTumorCellClass)

  sampInfo$Average.VAF <- as.character(sampInfo$Average.VAF)
  sampInfo$Average.VAF[which(sampInfo$Average.VAF %in% "not available")] <- "NA"

  sampInfo$KRAS.mut.1yes.0no <- factor(sampInfo$KRAS.mut.1yes.0no,
                                       levels = c("mutation in KRas",
                                                  "no mutation in KRas",
                                                  "not available"),
                                       labels = c("1", "0", "NA"))

  sampInfo$p53.mut.1yes.0no <- factor(sampInfo$p53.mut.1yes.0no,
                                      levels = c("mutation in TP53",
                                                 "no mutation in TP53",
                                                 "not available"),
                                      labels = c("1", "0", "NA"))

  sampInfo$CDKN2a.mut.1yes.0no <- factor(sampInfo$CDKN2a.mut.1yes.0no,
                                         levels = c("mutation in CDKN2a",
                                                    "no mutation in CDKN2a",
                                                    "not available"),
                                         labels = c("1", "0", "NA"))

  print(summary(sampInfo))

  movetonumeric <- numeric()
  for(i in 1:ncol(sampInfo)){
    if(i == which(names(sampInfo) %in% c("WholeTumorClass",
                                         "HighTumorCellClass"))){
     print(NULL)
    }
    else if(is.character(sampInfo[,i])){
      movetonumeric <- c(movetonumeric, i)
    }
  }

  for(column in movetonumeric){
    sampInfo[,column] <- as.numeric(sampInfo[,column])
  }

  sampInfo <- sampInfo[, c("Sample.name",
                           "Array.file",
                           "Clinical.center",
                           names(sampInfo)[((which(names(sampInfo) %in% "Sample.name")+1):
                              (which(names(sampInfo) %in% "Clinical.center")-1))],
                           names(sampInfo)[((which(names(sampInfo) %in% "Clinical.center")+1):
                              (which(names(sampInfo) %in% "Array.file")-1))])]

  print(summary(sampInfo))

  ## =============================
  # Wrangle survival data to consistent naming and days format
  ## =============================

  sampInfo$survivalA <- sampInfo$OS.delay.months * (365/12)
  sampInfo$censorA.0yes.1no <- sampInfo$OS.censor.0yes.1no

  sampInfo$survivalB <- sampInfo$DFS.delay.months * (365/12)
  sampInfo$censorB.0yes.1no <- sampInfo$DFS.censor.0yes.1no

  remove_cols <- which(names(sampInfo) %in% c("OS.delay.months",
                                             "ostime",
                                             "OS.censor.0yes.1no",
                                             "DFS.delay.months",
                                             "DFStime",
                                             "DFS.censor.0yes.1no"))

  sampInfo <- sampInfo[,-remove_cols]

  print(dim(sampInfo))

  ## =============================
  # Add metadata
  ## =============================

  metadata <- list(log.transformed = "FALSE",
                   reference = "Puleo F et al, Gastroenterology, 2018, PMID:30165049",
                   accession = "E-MTAB-6134",
                   description = "309 primary PDAC FFPE block tumor samples processed by Affymetrix array",
                   survivalA = "overall survival days",
                   survivalB = "disease-free survival days",
                   exp.type = "Array")

  ## =============================
  # Compile dataset and organize samples
  ## =============================

  dataset <- list(ex, sampInfo, featInfo, metadata)
  names(dataset) <- c("ex", "sampInfo", "featInfo", "metadata")

  matched <- match(names(dataset$ex), dataset$sampInfo$Sample.name)
  dataset$sampInfo <- dataset$sampInfo[matched,]

  Puleo_array <- dataset

  ## =============================
  # Consensus clustering to find tumor and stroma subtypes
  ## =============================

  # Tumor
  # ----------------------
  dataset <- Puleo_array
  sampleset <- 1:nrow(dataset$sampInfo)
  tmp.k <- 2
  tmp.ncusts <- 2

  featureset <- which(dataset$featInfo$SYMBOL %in%
                        c(as.character(pdac::gene_lists$Moffitt.Classical.25),
                          as.character(pdac::gene_lists$Moffitt.Basal.25)))

  smallx <- t(scale(t(dataset$ex[featureset,sampleset])))

  sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx),
                                                           seed = 1234,
                                                           pFeature = 0.8,
                                                           pItem = 0.8,
                                                           maxK = 6,
                                                           reps=200,
                                                           distance="pearson",
                                                           clusterAlg="kmdist")[[tmp.k]]$consensusTree

  tmp.cluster <- c("basal","classical")[cutree(tree = sampletree, k = 2)]
  dataset$sampInfo$MoffittTumor <- NA
  dataset$sampInfo$MoffittTumor[sampleset] <- tmp.cluster

  ColSideColors <-  getSideColors(sampInfo = dataset$sampInfo[sampleset,],
                                  sampleTracks = c("MoffittTumor"),
                                  colorlists = list(c("orange", "blue")),
                                  drop.levels = TRUE)

  RowSideColors <-  getSideColors(sampInfo = data.frame(basal =dataset$featInfo$SYMBOL[featureset] %in%
                                                          pdac::gene_lists$Moffitt.Basal.25,
                                                        classical =dataset$featInfo$SYMBOL[featureset] %in%
                                                          pdac::gene_lists$Moffitt.Classical.25),
                                  sampleTracks = c("basal",
                                                   "classical"),
                                  colorlists = list(c=c("white","orange"),
                                                    b=c("white","blue")))
  heatmap.3(x = smallx,
            scale="row",
            labRow = dataset$featInfo$SYMBOL[featureset],
            col = colorRampPalette(c("blue", "white", "red"))(n = 299),
            Colv = as.dendrogram(sampletree),
            Rowv = TRUE,
            distfun = function(x) as.dist((1-cor(t(x)))/2),
            ColSideColors = ColSideColors$SideColors,
            ColSideColorsSize = 6,
            RowSideColorsSize = 6,
            RowSideColors = t(RowSideColors$SideColors),
            margins = c(5,20))
  legend(xy.coords(x=.90,y=1),
         legend=c(ColSideColors$text),
         fill=c(ColSideColors$colors),
         border=FALSE, bty="n",
         y.intersp = 0.9, cex=0.5)

  # Stroma
  # ----------------------
  tmp.k <- 2
  tmp.ncusts <- 2

  featureset <- which(dataset$featInfo$SYMBOL %in%
                        c(as.character(pdac::gene_lists$Moffitt.Normal.25),
                          as.character(pdac::gene_lists$Moffitt.Activated.25)))

  smallx <- t(scale(t(dataset$ex[featureset,sampleset])))

  sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx),
                                                           seed = 1234,
                                                           pFeature = 0.8,
                                                           pItem = 0.8,
                                                           maxK = 6,
                                                           reps=200,
                                                           distance="euclidean",
                                                           clusterAlg="kmdist")[[tmp.k]]$consensusTree

  tmp.cluster <- c("activated","normal")[cutree(tree = sampletree, k = 2)]
  dataset$sampInfo$MoffittStroma <- NA
  dataset$sampInfo$MoffittStroma[sampleset] <- tmp.cluster

  ColSideColors <-  getSideColors(sampInfo = dataset$sampInfo[sampleset,],
                                  sampleTracks = c("MoffittStroma"),
                                  colorlists = list(c("brown", "lightblue")),
                                  drop.levels = TRUE)

  RowSideColors <-  getSideColors(sampInfo = data.frame(normal =dataset$featInfo$SYMBOL[featureset] %in%
                                                          pdac::gene_lists$Moffitt.Normal.25,
                                                        activated =dataset$featInfo$SYMBOL[featureset] %in%
                                                          pdac::gene_lists$Moffitt.Activated.25),
                                  sampleTracks = c("normal",
                                                   "activated"),
                                  colorlists = list(c=c("white","lightblue"),
                                                    b=c("white","brown")))
  heatmap.3(x = smallx,
            scale="row",
            labRow = dataset$featInfo$SYMBOL[featureset],
            col = colorRampPalette(c("blue", "white", "red"))(n = 299),
            Colv = as.dendrogram(sampletree),
            Rowv = TRUE,
            distfun = function(x) as.dist((1-cor(t(x)))/2),
            ColSideColors = ColSideColors$SideColors,
            ColSideColorsSize = 6,
            RowSideColorsSize = 6,
            RowSideColors = t(RowSideColors$SideColors),
            margins = c(5,20))
  legend(xy.coords(x=.90,y=1),
         legend=c(ColSideColors$text),
         fill=c(ColSideColors$colors),
         border=FALSE, bty="n",
         y.intersp = 0.9, cex=0.5)

  Puleo_array <- dataset

  ## =============================
  # Save dataset
  ## =============================
  save(list = c("Puleo_array"),
       file = "./data/Puleo_array.RData",
       compress = T)

  return(NULL)
}

