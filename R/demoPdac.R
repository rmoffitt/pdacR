#' Demo and testing script
#'
#' @export

demoPdac <- function() {
#library(pdac)

  data_set_list <- pdac::data_set_list
  for(selectedvariable in data_set_list$variablenames){
      print(as.character(selectedvariable))
      dataset <- get(as.character(selectedvariable))
      print(str(dataset$sampInfo))
  }
}
