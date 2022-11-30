#' This code is here for transparency, but for data storage reasons the full combined rds file is not available.
#'
#' @import Seurat
#' @import readxl
#' @export
#'

parse_scOh_atlas <- function(){
  ## Won't be available as is, here for transparency and posterity
  atlas <- readRDS("/premDisk/khoh/Github/misc/data/pdac_full_combined.rds")
  atlas@meta.data = atlas@meta.data[,which(colnames(atlas@meta.data) %in%
                                             c("Patient",
                                               "Dataset")
  )]
  atlas@meta.data$Patient[which(is.na(atlas$Patient))] = "Powers1"
  tmp_converted = Convert_GUI_data(atlas, from = "Seurat", to ="Moffitt",
                                   pseudobulk = T, bulk_by = "Patient")
  tmp_converted$metadata$log.transformed = FALSE
  rownames(tmp_converted$sampInfo) = NULL

  SupplementaryTable4Metadata <- read_excel("./inst/extdata/scAtlas/SupplementaryTable4Metadata.xlsx")
  SupplementaryTable4Metadata  = SupplementaryTable4Metadata[,c("Patient","Condition","Note")]

  tmp_converted$sampInfo = dplyr::left_join(tmp_converted$sampInfo, SupplementaryTable4Metadata, by = "Patient")
  tmp_converted$sampInfo = tmp_converted$sampInfo[,order(colnames(tmp_converted$sampInfo))]
  tmp_converted$sampInfo = tmp_converted$sampInfo[,c(which(colnames(tmp_converted$sampInfo) == "Patient"),
                                                     which(colnames(tmp_converted$sampInfo) != "Patient"))]
  tmp_converted$metadata$default_selections = list(filter_column = "Condition",
                                                   filter_levels = "Condition:Normal",
                                                   sampleTracks = "Dataset")
  saveRDS(tmp_converted,
          file = "./data/scAtlas.pseudobulked.rds",
          compress = T)

}
