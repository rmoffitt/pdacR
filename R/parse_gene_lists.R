#' Parse data from tables
#'
#' @export
#' @import org.Hs.eg.db
#' @import openxlsx

parse_gene_lists <- function() {
  library(org.Hs.eg.db)
  library(openxlsx)
  gene_lists <- list()

  ## =============================
  # from Moffitt et al, Nature Genetics, 2015
  ## =============================
  tmpFname <- system.file("extdata/gene_lists", "Moffitt-basal25.txt", package = "pdacR")
  MoffittBasal25 <- read.table(file = tmpFname,
                               sep = "\t",
                               header = FALSE,
                               skip = 0)
  gene_lists$Moffitt.Basal.25 <- MoffittBasal25[[1]]
  # --------------------------
  tmpFname <- system.file("extdata/gene_lists", "Moffitt-classical25.txt", package = "pdacR")
  MoffittClassical25 <- read.table(file = tmpFname,
                                   sep = "\t",
                                   header = FALSE,
                                   skip = 0)
  gene_lists$Moffitt.Classical.25 <- MoffittClassical25[[1]]
  # --------------------------
  tmpFname <- system.file("extdata/gene_lists", "Moffitt-normal25.txt", package = "pdacR")
  MoffittNormal25 <- read.table(file = tmpFname,
                                sep = "\t",
                                header = FALSE,
                                skip = 0)
  gene_lists$Moffitt.Normal.25 <- MoffittNormal25[[1]]
  # --------------------------
  tmpFname <- system.file("extdata/gene_lists", "Moffitt-activated25.txt", package = "pdacR")
  MoffittActivated25 <- read.table(file = tmpFname,
                                   sep = "\t",
                                   header = FALSE,
                                   skip = 0)
  gene_lists$Moffitt.Activated.25 <- MoffittActivated25[[1]]

  gene_lists$Moffitt.Tumor <- rbind(data.frame(symbols = MoffittBasal25[[1]], type = "basal"),
                                    data.frame(symbols = MoffittClassical25[[1]], type = "classical"))

  ## =============================
  # from Puleo et al, Gastroenterology, 2018
  ## =============================
  tmpFname <- system.file("extdata/Puleo", "1-s2.0-S0016508518349199-mmc3.xlsx", package = "pdacR")
  Puleo_centroids <-  read.xlsx(tmpFname)

  gene_lists$Puleo.Centroids <- data.frame(
    symbols = Puleo_centroids[,1],
    type = apply(X = Puleo_centroids[,-1],
                 MARGIN = 1,
                 FUN = function(x){ names(x)[which.max(x)] })
  )
  write.table(x = gene_lists$Puleo.Centroids,
              file = "./data/Puleo_readable_list.tsv",
              sep ="\t",
              row.names = FALSE,
              col.names = FALSE)


  ## =============================
  # from Collisson et al, Nature Medicine, 2011
  ## =============================
  tmpFname <- system.file("extdata/gene_lists", "Collisson-classical.csv", package = "pdacR")
  CollissonClassical <- read.table(file = tmpFname,
                                   sep = ",",
                                   header = FALSE,
                                   skip = 0)
  gene_lists$Collisson.Classical <- CollissonClassical[[1]]
  # --------------------------
  tmpFname <- system.file("extdata/gene_lists", "Collisson-exocrine.csv", package = "pdacR")
  CollissonExocrine<- read.table(file = tmpFname,
                                 sep = ",",
                                 header = FALSE,
                                 skip = 0)
  gene_lists$Collisson.Exocrine <- CollissonExocrine[[1]]
  # --------------------------
  tmpFname <- system.file("extdata/gene_lists", "Collisson-QM.csv", package = "pdacR")
  CollissonQM <- read.table(file = tmpFname,
                            sep = ",",
                            header = FALSE,
                            skip = 0)
  gene_lists$Collisson.QM <- CollissonQM[[1]]

  ## =============================
  # from Bailey et al, Nature, 2016
  ## =============================
  tmpFname <- system.file("extdata/gene_lists", "ICGC_ADEX_vs_rest.txt", package = "pdacR")
  ICGC_ADEX <- read.table(file = tmpFname,
                          sep = "\t",
                          header = TRUE,
                          skip = 0)
  gene_lists$ICGC.ADEX.Up <- ICGC_ADEX$Symbol[ICGC_ADEX$t>0]
  gene_lists$ICGC.ADEX.Down <- ICGC_ADEX$Symbol[ICGC_ADEX$t<0]
  # --------------------------
  tmpFname <- system.file("extdata/gene_lists", "ICGC_Immunogenic_vs_rest.txt", package = "pdacR")
  ICGC_Immunogenic <- read.table(file = tmpFname,
                                 sep = "\t",
                                 header = TRUE,
                                 skip = 0)
  gene_lists$ICGC.Immunogenic.Up <- ICGC_Immunogenic$Symbol[ICGC_Immunogenic$t>0]
  gene_lists$ICGC.Immunogenic.Down <- ICGC_Immunogenic$Symbol[ICGC_Immunogenic$t<0]
  # --------------------------
  tmpFname <- system.file("extdata/gene_lists", "ICGC_Progenitor_vs_rest.txt", package = "pdacR")
  ICGC_Progenitor <- read.table(file = tmpFname,
                                sep = "\t",
                                header = TRUE,
                                skip = 0)
  gene_lists$ICGC.Progenitor.Up <- ICGC_Progenitor$Symbol[ICGC_Progenitor$t>0]
  gene_lists$ICGC.Progenitor.Down <- ICGC_Progenitor$Symbol[ICGC_Progenitor$t<0]
  # --------------------------
  tmpFname <- system.file("extdata/gene_lists", "ICGC_Squamous_vs_rest.txt", package = "pdacR")
  ICGC_Squamous <- read.table(file = tmpFname,
                              sep = "\t",
                              header = TRUE,
                              skip = 0)
  gene_lists$ICGC.Squamous.Up <- ICGC_Squamous$Symbol[ICGC_Squamous$t>0]
  gene_lists$ICGC.Squamous.Down <- ICGC_Squamous$Symbol[ICGC_Squamous$t<0]
  # --------------------------
  tmpFname <- system.file("extdata/gene_lists", "ICGC_SAM.txt", package = "pdacR")
  ICGC_SAM <- read.table(file = tmpFname,
                         sep = "\t",
                         header = FALSE,
                         skip = 0)
  checkICGC <- function(x,A,I,S,P){
    if((x %in% A) & !(x %in% I) & !(x %in% S) & !(x %in% P) ){return("ADEX")}
    if(!(x %in% A) & (x %in% I) & !(x %in% S) & !(x %in% P) ){return("Immunogenic")}
    if(!(x %in% A) & !(x %in% I) & (x %in% S) & !(x %in% P) ){return("Squamous")}
    if(!(x %in% A) & !(x %in% I) & !(x %in% S) & (x %in% P) ){return("Progenitor")}
    return("not unique")
  }
  gene_lists$ICGC.SAM <- data.frame(symbols = ICGC_SAM[,2],
                                    type = sapply(X = ICGC_SAM[,2],
                                                  FUN = function(x){checkICGC(x,
                                                                              gene_lists$ICGC.ADEX.Up,
                                                                              gene_lists$ICGC.Immunogenic.Up ,
                                                                              gene_lists$ICGC.Squamous.Up ,
                                                                              gene_lists$ICGC.Progenitor.Up )})
  )

  # --------------------------
  for(i in c("ADEX","Immunogenic","Squamous","Progenitor")){
    shortlist = gene_lists$ICGC.SAM$symbols[gene_lists$ICGC.SAM$type == i]
    gene_lists[[paste("ICGC",i,"Up.unique",sep = ".")]] <- intersect(gene_lists[[paste("ICGC",i,"Up",sep = ".")]],shortlist) ## Parse down to Unique Up
    gene_lists[[paste("ICGC",i,"Down.unique",sep = ".")]] <- intersect(gene_lists[[paste("ICGC",i,"Down",sep = ".")]],shortlist) ## Parse down to Unique Down
  }
  write.table(x = gene_lists$ICGC.SAM,
              file = "./data/Bailey_readable_list.tsv",
              sep ="\t",
              row.names = FALSE,
              col.names = FALSE)

  ## =============================
  # from Scarpa et al, Nature, 2017
  ## =============================
  tmpFname <- system.file("extdata/gene_lists", "Scarpa_PanNET_genes.txt", package = "pdacR")
  Scarpa_PanNET <- read.table(file = tmpFname,
                              sep = "\t",
                              header = FALSE,
                              skip = 0)
  gene_lists$Scarpa_PanNET <- Scarpa_PanNET[[1]]


  ## =============================
  # from Moffitt et al, Nature Genetics, 2015
  ## =============================
  tmpFname <- system.file("extdata/gene_lists", "ng.3398-S4.xlsx", package = "pdacR")
  MoffittFullTable <- read.xlsx(xlsxFile =  tmpFname)
  gene_lists$Moffitt.Top5s <- data.frame(SYMBOL = c(),
                                         class = c())
  gene_lists$Moffitt.Top25s <- data.frame(SYMBOL = c(),
                                          class = c())
  for(i in 1:14){
    tmp.score <-
      MoffittFullTable[,i+2] -
      apply(X = MoffittFullTable[,c(-1,-2,-(i+2))],
            MARGIN = 1,
            FUN=max)
    rank.order <- order(-tmp.score)
    gene_lists$Moffitt.Top5s <-
      rbind(gene_lists$Moffitt.Top5s ,
            data.frame(SYMBOL = (MoffittFullTable$symbol[rank.order[1:5]]),
                       class = (names(MoffittFullTable)[i+2])))
    gene_lists$Moffitt.Top25s <-
      rbind(gene_lists$Moffitt.Top25s ,
            data.frame(SYMBOL = MoffittFullTable$symbol[rank.order[1:25]],
                       class = names(MoffittFullTable)[i+2]) )
    gene_lists[[paste(c("Moffitt",
                        names(MoffittFullTable)[i+2],
                        "top25"),
                      collapse=".")]] <-
      MoffittFullTable$symbol[rank.order[1:25]]
    gene_lists[[paste(c("Moffitt",
                        names(MoffittFullTable)[i+2],
                        "top250"),
                      collapse=".")]] <-
      MoffittFullTable$symbol[rank.order[1:250]]
    gene_lists[[paste(c("Moffitt",
                        names(MoffittFullTable)[i+2],
                        "top100"),
                      collapse=".")]] <-
      MoffittFullTable$symbol[rank.order[1:100]]
  }

  ## =============================
  # CIBERSORT from Newman et al, Nature Methods, 2015
  ## =============================
  tmpFname <- system.file("extdata/gene_lists", "cibersort.tsv", package = "pdacR")
  cibersort <- read.table(file = tmpFname,
                          sep = "\t",
                          header = TRUE,
                          skip = 0)
  for(i in 2:length(cibersort)){
    print(names(cibersort)[i])
    gene_lists[[paste(c("CIBERSORT",names(cibersort)[i]),collapse=".")]] <-
      as.character(cibersort[cibersort[,i]==1,1])
  }

  ## =============================
  # from Chan_Seng_Yue et al, Nature Genetics, 2019
  ## =============================
  tmpFname <- system.file("extdata/gene_lists", "Notta_supplement.xlsx", package = "pdacR")
  Chan_Seng_Yue <- read.xlsx(tmpFname,
                             sheet = "Supplementary Table 4",
                             startRow = 2)
  Chan_Seng_Yue = Chan_Seng_Yue[,c("Sig..2.genes","Sig..10.genes","Sig..1.genes","Sig..6.genes")]
  names(Chan_Seng_Yue) = c("BasalA","BasalB","ClassicalA","ClassicalB")

  for(i in 1:ncol(Chan_Seng_Yue)){
    gene_lists[[paste0("Chan_Seng_Yue.",names(Chan_Seng_Yue)[i])]] <- as.character(na.omit(Chan_Seng_Yue[,i]))
  }

  ## =============================
  # Save gene lists
  ## =============================
  print(summary(gene_lists))
  saveRDS(gene_lists,
          file = "./data/gene_lists.rds")
  return(NULL)
}
