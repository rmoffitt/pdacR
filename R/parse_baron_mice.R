#' Parse data for baron
#'
#' @export
#' @import SingleCellExperiment
parse_baron_mice <- function() {
   library(SingleCellExperiment)
   baron <- readRDS("~/pdac/inst/extdata/Mouse_scRNA/baron-mouse.rds")

   Baron_Dataset <- list()
   Baron_Dataset$ex <- counts(baron)
   Baron_Dataset$sampInfo <- data.frame(barcodes = colnames(baron),
                                     disease = "PDX1/KRAS/p53")
   samp <- data.frame(colData(baron)@listData)[,1:3]
   Baron_Dataset$sampInfo <- cbind(Baron_Dataset$sampInfo, samp)
   Baron_Dataset$sampInfo$mouse <- factor(Baron_Dataset$sampInfo$mouse)
   Baron_Dataset$metadata = list(log.transformed = F,
                                 exp.type = "scRNA",
                              reference = "Baron M, et al, Cell Systems, 2016",
                              accession = "doi: 10.1016/j.cels.2016.08.011")
   Baron_Dataset$featInfo = data.frame(SYMBOL = rownames(baron))

   save(list = c("Baron_Dataset"),
      file = "./data/Baron_Dataset.RData")
   return(NULL)
}
