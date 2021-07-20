#' Parse data from tables
#'
#' @export
#' @import SingleCellExperiment
parse_GSE125588_mice <- function() {
   # This will take quite a while in current format
   library(SingleCellExperiment)
   # ------------------------------------------------------------------------------------
   norm_panc <- system.file("extdata/Mouse_scRNA", "GSM3577882_parsed.rds", package = "pdac")
   early_KIC <- system.file("extdata/Mouse_scRNA", "GSM3577883_parsed.rds", package = "pdac")
   late_KIC <- system.file("extdata/Mouse_scRNA", "GSM3577884_parsed.rds", package = "pdac")
   late_KPfC <- system.file("extdata/Mouse_scRNA", "GSM3577885_parsed.rds", package = "pdac")
   late_KPC <- system.file("extdata/Mouse_scRNA", "GSM3577886_parsed.rds", package = "pdac")
   # Norm --------------------------------------------------------
   Normal_Pancreas <- list()
   Normal_Pancreas$ex <- as.matrix(readRDS(norm_panc))
   Normal_Pancreas$sampInfo <- data.frame(barcodes = colnames(readRDS(norm_panc)),
                                    disease = "normal pancreas")
   Normal_Pancreas$metadata = list(log.transformed = F,
                                   exp.type = "scRNA",
                                    reference = "Hosein & Brekken, MD Anderson, 2019",
                                    accession = "GSE125588")
   Normal_Pancreas$featInfo = data.frame(SYMBOL = rownames(readRDS(norm_panc)))

   a <- which(duplicated(Normal_Pancreas$featInfo))
   b <- as.vector(Normal_Pancreas$featInfo$SYMBOL[a])
   b <- b[!duplicated(b)]
   count = 0

   for (i in b){
      count = count + 1
      print(i)
      g <- which(Normal_Pancreas$featInfo$SYMBOL == i)
      print(g)
      Normal_Pancreas$ex[g[1],] <- colSums(Normal_Pancreas$ex[g,])
      Normal_Pancreas$ex <- Normal_Pancreas$ex[-(g[-1]),]
      Normal_Pancreas$featInfo <- data.frame(SYMBOL = Normal_Pancreas$featInfo$SYMBOL[-(g[-1])])
      print(paste(count, "duplicates trimmed out of", length(b)))
      print("---------------------------------")
   }


   # early KIC --------------------------------------------------------
   Early_KIC <- list()
   Early_KIC$ex <- as.matrix(readRDS(early_KIC))
   Early_KIC$sampInfo <- data.frame(barcodes = colnames(readRDS(early_KIC)),
                                    disease = "early p48/KRAS/CDKN2A")
   Early_KIC$metadata = list(log.transformed = F,
                             exp.type = "scRNA",
                              reference = "Hosein & Brekken, MD Anderson, 2019",
                              accession = "GSE125588")
   Early_KIC$featInfo = data.frame(SYMBOL = rownames(readRDS(early_KIC)))

   a <- which(duplicated(Early_KIC$featInfo))
   b <- as.vector(Early_KIC$featInfo$SYMBOL[a])
   b <- b[!duplicated(b)]
   count = 0

   for (i in b){
      count = count + 1
      print(i)
      g <- which(Early_KIC$featInfo$SYMBOL == i)
      print(g)
      Early_KIC$ex[g[1],] <- colSums(Early_KIC$ex[g,])
      Early_KIC$ex <- Early_KIC$ex[-(g[-1]),]
      Early_KIC$featInfo <- data.frame(SYMBOL = Early_KIC$featInfo$SYMBOL[-(g[-1])])
      print(paste(count, "duplicates trimmed out of", length(b)))
      print("---------------------------------")
   }

   # late KIC --------------------------------------------------------
   Late_KIC <- list()
   Late_KIC$ex <- as.matrix(readRDS(late_KIC))
   Late_KIC$sampInfo <- data.frame(barcodes = colnames(readRDS(late_KIC)),
                                    disease = "late p48/KRAS/CDKN2A")
   Late_KIC$metadata = list(log.transformed = F,
                            exp.type = "scRNA",
                             reference = "Hosein & Brekken, MD Anderson, 2019",
                             accession = "GSE125588")
   Late_KIC$featInfo = data.frame(SYMBOL = rownames(readRDS(late_KIC)))

   a <- which(duplicated(Late_KIC$featInfo))
   b <- as.vector(Late_KIC$featInfo$SYMBOL[a])
   b <- b[!duplicated(b)]
   count = 0

   for (i in b){
      count = count + 1
      print(i)
      g <- which(Late_KIC$featInfo$SYMBOL == i)
      print(g)
      Late_KIC$ex[g[1],] <- colSums(Late_KIC$ex[g,])
      Late_KIC$ex <- Late_KIC$ex[-(g[-1]),]
      Late_KIC$featInfo <- data.frame(SYMBOL = Late_KIC$featInfo$SYMBOL[-(g[-1])])
      print(paste(count, "duplicates trimmed out of", length(b)))
      print("---------------------------------")
   }

   # Late KPfC --------------------------------------------------------
   Late_KPfC <- list()
   Late_KPfC$ex <- as.matrix(readRDS(late_KPfC))
   Late_KPfC$sampInfo <- data.frame(barcodes = colnames(readRDS(late_KPfC)),
                              disease = "p48/KRAS/p53")
   Late_KPfC$metadata = list(log.transformed = F,
                              reference = "Hosein & Brekken, MD Anderson, 2019",
                              accession = "GSE125588")
   Late_KPfC$featInfo = data.frame(SYMBOL = rownames(readRDS(late_KPfC)))

   a <- which(duplicated(Late_KPfC$featInfo))
   b <- as.vector(Late_KPfC$featInfo$SYMBOL[a])
   b <- b[!duplicated(b)]
   count = 0

   for (i in b){
      count = count + 1
      print(i)
      g <- which(Late_KPfC$featInfo$SYMBOL == i)
      print(g)
      Late_KPfC$ex[g[1],] <- colSums(Late_KPfC$ex[g,])
      Late_KPfC$ex <- Late_KPfC$ex[-(g[-1]),]
      Late_KPfC$featInfo <- data.frame(SYMBOL = Late_KPfC$featInfo$SYMBOL[-(g[-1])])
      print(paste(count, "duplicates trimmed out of", length(b)))
      print("---------------------------------")
   }

   # Late KPC --------------------------------------------------------
   Late_KPC <- list()
   Late_KPC$ex <- as.matrix(readRDS(late_KPC))
   Late_KPC$sampInfo <- data.frame(barcodes = colnames(readRDS(late_KPC)),
                             disease = "PDX1/KRAS/p53")
   Late_KPC$metadata = list(log.transformed = F,
                            exp.type = "scRNA",
                             reference = "Hosein & Brekken, MD Anderson, 2019",
                             accession = "GSE125588")
   Late_KPC$featInfo = data.frame(SYMBOL = rownames(readRDS(late_KPC)))

   a <- which(duplicated(Late_KPC$featInfo))
   b <- as.vector(Late_KPC$featInfo$SYMBOL[a])
   b <- b[!duplicated(b)]
   count = 0

   for (i in b){
      count = count + 1
      print(i)
      g <- which(Late_KPC$featInfo$SYMBOL == i)
      print(g)
      Late_KPC$ex[g[1],] <- colSums(Late_KPC$ex[g,])
      Late_KPC$ex <- Late_KPC$ex[-(g[-1]),]
      Late_KPC$featInfo <- data.frame(SYMBOL = Late_KPC$featInfo$SYMBOL[-(g[-1])])
      print(paste(count, "duplicates trimmed out of", length(b)))
      print("---------------------------------")
   }

   # Combined GSE125588 ---------------------------------------------
   GSE_125588 <- list()
   GSE_125588$ex <- cbind(Normal_Pancreas$ex, Early_KIC$ex, Late_KIC$ex, Late_KPfC$ex, Late_KPC$ex)
   GSE_125588$sampInfo <- rbind(Normal_Pancreas$sampInfo, Early_KIC$sampInfo, Late_KIC$sampInfo, Late_KPfC$sampInfo, Late_KPC$sampInfo)
   GSE_125588$metadata <- list(log.transformed = F,
                                    reference = "Hosein & Brekken, MD Anderson, 2019",
                                    accession = "GSE125588")
   GSE_125588$featInfo = data.frame(SYMBOL = Late_KPC$featInfo$SYMBOL)

   ## =============================
   # Save dataset
   ## =============================
   save(list = c("Normal_Pancreas"),
        file = "./data/Normal_Pancreas.RData")

   save(list = c("Early_KIC"),
        file = "./data/Early_KIC.RData")

   save(list = c("Late_KIC"),
        file = "./data/Late_KIC.RData")

   save(list = c("Late_KPfC"),
        file = "./data/Late_KPfC.RData")

   save(list = c("Late_KPC"),
        file = "./data/Late_KPC.RData")

    save(list = c("GSE_125588"),
         file = "./data/Combined_scRNA.RData")
   return(NULL)
}

