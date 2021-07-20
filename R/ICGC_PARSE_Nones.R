#' Download Chen dataset from GEO
#'
#' @export
#' @import GEOquery
#' @import Biobase
#' @import readr
#' @import rlang

ICGC_PARSE_Nones <- function() {

  ## =============================
survival <- read.xlsx(system.file("extdata/Nones", "ICGC_GEO_Map_05112015.xlsx", package = "pdac"))

x <- getGEO(GEO = "GSE50827",
            GSEMatrix =TRUE,
            AnnotGPL=FALSE,
            getGPL=TRUE)
## =============================
x-> gset
gset <- gset$GSE50827_series_matrix.txt.gz
samps <-gset@phenoData
sampInfo <- samps@data
sampInfo <- sampInfo[,c("title","geo_accession","source_name_ch1",
                        "relation", "disease:ch1", "tissue:ch1")]
names(sampInfo)[names(sampInfo)=="title"] <- "submitted_sample_id"
sampInfo$specimen_type <- "primaryPDAC"
sampInfo$location <- "pancreas"
sampInfo$stage <- NA
sampInfo$survival_months <- NA
sampInfo$censor <- NA
sampInfo$specimen_type <- as.factor(sampInfo$specimen_type)
sampInfo$location <- as.factor(sampInfo$location)
sampInfo$stage <- as.factor(sampInfo$stage)
sampInfo$censor <- survival$Moffitt_censor_call
sampInfo$survival_days <- survival$Moffitt_survival_days
sampInfo$survivalA <- survival$Moffitt_survival_days
sampInfo$censorA.1yes.0no <- survival$Moffitt_censor_call

## =============================
ex <- Biobase::exprs(gset)
featInfo <- data.frame(PROBES = row.names(ex),
                       ENTREZID = row.names(ex),
                       SYMBOL = row.names(ex))

## =============================
gpl <- read_delim(system.file("extdata/Nones", "GPL10558.soft", package = "pdac"),
                  "\t", escape_double = FALSE, trim_ws = TRUE,
                  skip = 165479)

gpl <- as.data.frame(lapply(gpl,as.factor))

neworder <- match(table = gpl$ID,x = featInfo$PROBES)
gpl <- gpl[neworder,]

featInfo <- data.frame(PROBES = gpl$ID,
                       ENTREZID = gpl$Entrez_Gene_ID,
                       SYMBOL = gpl$Symbol)

#CHECK
plot(x = ex[featInfo$SYMBOL %in% "S100A2",],
     y = ex[featInfo$SYMBOL %in% "KRT17",])

which(featInfo$SYMBOL %in% "S100A2")
#CHECK
## =============================
# aggregation by entrez and symbol

print("Dimension before aggregation")
print(dim(ex))

tmp <- aggregate(x = ex,
                 by = featInfo[,1:3],
                 FUN = function(x) sum(x))

featInfo <- data.frame(SYMBOL = tmp[[3]],
                       ENTREZID = tmp[[2]])
ex  <- tmp[,c(-1,-2,-3)] #taking out PROBES, SYMBOLS, ENTREZID from tms so its only gsms

print(dim(ex))
## =============================
#metadata <- list(log.transformed = TRUE)
metadata <- list(log.transformed = TRUE,
                 reference = "Nones et al, GUT BMJ 2015, PMID: 25336113 ",
                 accession = "GEO: GSE50827",
                 description = "SOX9 shown to stimulate ductal gene ex and acceleration of premalignant lesions preceding pdac",
                 survivalA = "overall survival days",
                 survivalB = "None")




## =============================
sampInfo$tumor.classifier.training <- FALSE
sampInfo$stroma.classifier.training <- FALSE
sampInfo[["cluster.MT"]] <- as.character(NA)
sampInfo[["cluster.MS"]] <- as.character(NA)

## =============================
sampleset <- which(sampInfo$specimen_type %in% "primaryPDAC")
featureset <- which(featInfo$SYMBOL %in%
                      c(as.character(pdac::gene_lists$Moffitt.Basal.25),
                        as.character(pdac::gene_lists$Moffitt.Classical.25)))
## =============================
small_x <- t(scale(scale=FALSE,
                  center=TRUE,
                  x=t((ex[featureset,sampleset]))))

sampletree <- ConsensusClusterPlus::ConsensusClusterPlus(d = small_x,
                                                         seed = 1234,
                                                         pFeature = 0.8,
                                                         pItem = 0.8,
                                                         maxK = 6,
                                                         reps=200,
                                                         distance="pearson",
                                                         clusterAlg="kmdist")[[2]]$consensusTree
tmp.cluster <- c("basal","classical")[cutree(tree = sampletree, k = 2)]
sampInfo$cluster.MT[sampleset] <- tmp.cluster
sampInfo$tumor.classifier.training[sampleset] <- TRUE
ggplot(data = data.frame(classical = (colMeans(ex[featInfo$SYMBOL %in%
                                                    pdac::gene_lists$Moffitt.Classical.25,sampleset])),
                         basal = (colMeans(ex[featInfo$SYMBOL %in%
                                                pdac::gene_lists$Moffitt.Basal.25,sampleset])),
                         cluster.MT = sampInfo$cluster.MT[sampleset]),
       aes(x = basal,y = classical, color = cluster.MT)) +
  geom_point(size=2)
#############
## =============================
Nones_GEO_array <- list(sampInfo=sampInfo,featInfo=featInfo,ex=ex,metadata=metadata)
save(list = c("Nones_GEO_array"),
     file = "./data/Nones_GEO_array.RData")

return(NULL)
}
