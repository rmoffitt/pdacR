#' generate PAMG for given dataset
#'
#' @export
#' @import pdacmolgrad

implement_PAMG = function(data){
    temp = pdacmolgrad::projectMolGrad(log2(1+data$ex), geneSymbols = data$featInfo$SYMBOL)
    names(temp) <- paste0("molgrad_",names(temp))
    #data$sampInfo = cbind(data$sampInfo, temp)
    return(as.list(temp))
}
