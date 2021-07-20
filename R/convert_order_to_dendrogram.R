#' create a flat dendrogram based on a ordered list of integers
#'
#' @export

convert_order_to_dendrogram <- function(sampleorder){
  a <- list()  
  for(i in 2:length(sampleorder)){
    if(i == 2){
      a$merge = matrix(-sampleorder[1:2],nrow=1,ncol=2)
    }else {
      a$merge <- rbind(a$merge,c(i-2,-sampleorder[i]))
    }
  }
  a$height <- rep(0,length(sampleorder)-1) 
  a$order  <- sampleorder       
  a$labels <- sampleorder
  class(a) <- "hclust"
  d<- as.dendrogram(a)
  return(d)
}