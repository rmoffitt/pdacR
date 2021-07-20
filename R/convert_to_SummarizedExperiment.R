#' Create SummarizedExperiment object from pdac list style object
#'
#' @export
#' @import BiocManager

convert_to_SummarizedExperiment <- function(dataset) {
  # SummarizedExperiment is a matrix-like container where rows represent
  # features of interest (e.g. genes, transcripts, exons, etc.) and columns
  # represent samples. The objects contain one or more assays, each
  # represented by a matrix-like object of numeric or other mode. The rows
  # of a SummarizedExperiment object represent features of interest.
  # Information about these features is stored in a DataFrame object,
  # accessible using the function rowData(). Each row of the  DataFrame
  # provides information on the feature in the corresponding row of the
  # SummarizedExperiment object. Columns of the DataFrame represent different
  # attributes of the features of interest, e.g., gene or transcript IDs, etc.

  # Analagous parts :
  # rowData(se)   <- dataset$featInfo
  # colData(se)   <- dataset$sampInfo
  # assays(se)    <- dataset$ex
  # metadata(se)  <- dataset$metadata

  nrows <- dim(dataset$ex)[1]
  ncols <- dim(dataset$ex)[2]

  expression <- matrix(dataset$ex, nrows)
  colData    <- DataFrame(dataset$sampInfo)
  rowData    <- DataFrame(dataset$featInfo)

  se <- SummarizedExperiment(assays  = list(expression=expression),
                       rowData  = rowData,
                       colData  = colData,
                       metadata = dataset$metadata)

  print(se)
  return(se)
}
