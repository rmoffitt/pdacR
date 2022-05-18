#' This code is here for transparency, but for data storage reasons the full combined rds file is not available.
#'
#' @import Seurat
#' @export
#'

parse_scOh_atlas <- function(){
  ## Won't be available as is, here for transparency and posterity
  atlas <- readRDS("/premDisk/khoh/Github/scTME/data/pdac_full_combined.rds")
  atlas@meta.data = atlas@meta.data[,which(colnames(atlas@meta.data) %in%
                                         c("Patient",
                                           "Dataset")
  )]
  atlas@meta.data$Patient[which(is.na(atlas$Patient))] = "Powers1"
  tmp_converted = Convert_GUI_data(atlas, from = "Seurat", to ="Moffitt",
                                   pseudobulk = T, bulk_by = "Patient")
  tmp_converted$metadata$log.transformed = FALSE
  rownames(tmp_converted$sampInfo) = NULL

  tmp$sampInfo = dplyr::left_join(tmp$sampInfo, SupplementaryTable4Metadata, by = "Patient")
  tmp$sampInfo = tmp$sampInfo[,order(colnames(tmp$sampInfo))]
  saveRDS(tmp_converted,
          file = "./data/scAtlas.pseudobulked.rds",
          compress = T)

  ## If we want to save dataset level objects
  # for(i in levels(as.factor(atlas$Dataset))){
  #   cells = which(atlas$Dataset == i)
  #   tmp = subset(atlas, cells = cells)
  #   tmp@meta.data = tmp@meta.data[,which(colnames(tmp@meta.data) %in%
  #                                         c("Patient",
  #                                           "Dataset")
  #                                       )]
  #   tmp_converted = Convert_GUI_data(atlas, from = "Seurat", to ="Moffitt",
  #                                    pseudobulk = T, bulk_by = "Patient")
  #   tmp_converted$metadata$log.transformed = FALSE
  #   rownames(tmp_converted$sampInfo) = NULL
  #   saveRDS(tmp_converted,
  #           file = paste0("/datadrive/shared/pdacR_datasets/",i,".pseudobulked.rds"),
  #           compress = T)
  # }
}
