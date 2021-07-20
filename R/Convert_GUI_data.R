#' Fully flexible conversion of .RData
#'
#' \code{Convert_GUI_data} will convert between a standard list format experiment file, SummarizedExperiment class, and Seurat Object
#' @param object An object of either class 'list','SummarizedExperiment', or 'Seurat'
#' @param from Character string of length one - the class of the object you're putting in
#' @param to Character string of length one - the class of the object you're getting
#' @return Returns an object of classyou requested that contains all info from inital experiment
#'
#' @export
#' @import SummarizedExperiment
#' @import Seurat

Convert_GUI_data <- function(object, from, to) {
  # From SE to either GUI List or Seurat
  if(from %in% c('Summ','summ','Summarized','SummarizedExperiment')){
    if(class(object) != "SummarizedExperiment"){
      warning("Object not a Summarized Experiment, check again")
    }
    else{
      if(to %in% c("GUI","gui","Moff","moff","moffitt","Moffitt","list")){
        new_obj <- list(ex = assay(object),
                        sampInfo = colData(object),
                        featInfo = rowData(object),
                        metaData = list(log.transformed = F))
      }
      if(to %in% c("Seurat","seurat")){
        new_obj <- Seurat::CreateSeuratObject(counts = assay(object),
                                              meta.data = colData(object))
      }
    }
  }
  # From GUI List to either Seurat or SE
  if(from %in% c("GUI","gui","Moff","moff","moffitt","Moffitt","list")){
    if(class(object) != "list"){
      warning("Object not a GUI list, check again")
    }
    else{
      if(to %in% c('Summ','summ','Summarized','SummarizedExperiment')){
        new_obj <- SummarizedExperiment(assays = object$ex,
                                        rowData = object$featInfo,
                                        colData = object$sampInfo,
                                        metadata = object$metadata)
      }
      if(to %in% c("Seurat","seurat")){
        warning('This is only helpful for single cell data, is your experiment single cell?')
        if(is.null(rownames(object$ex))){
          rownames(object$ex) = object$featInfo[,1]
        }
        if(is.null(rownames(object$sampInfo))){
          rownames(object$sampInfo) = colnames(object$ex)
        }
        new_obj <- Seurat::CreateSeuratObject(counts = object$ex,
                                              meta.data = object$sampInfo)
      }
    }
  }
  ## From Seurat to SE or GUI List
  if(from %in% c("Seurat","seurat")){
    if(class(object) != "Seurat"){
      warning("Object not a Seurat Object, check again")
    }
    else{
      if(to %in% c('Summ','summ','Summarized','SummarizedExperiment')){
        new_obj <- SummarizedExperiment(assays = Seurat::GetAssayData(object),
                                        rowData = rownames(object),
                                        colData = object@meta.data,
                                        metadata = object@Misc)
      }
      if(to %in% c("GUI","gui","Moff","moff","moffitt","Moffitt","list")){
        new_obj <- list(ex = as.matrix(Seurat::GetAssayData(object)),
                        sampInfo = object@meta.data,
                        featInfo = rownames(object),
                        metadata = list(object@Misc))
      }
    }
  }
  return(new_obj)
}
