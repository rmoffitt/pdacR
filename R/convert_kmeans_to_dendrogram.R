#' converts a list of integer labels to a pseudo-tree for plotting above or along side a heat map
#' Authors:  Richard Moffitt, Steven Botez, Peter Ly
#'
#' @export

convert_kmeans_to_dendrogram <- function(sample_order) {
  a <- list()
  a$merge <- matrix(nrow = 0, ncol = 2)
  what_row <- 1
  k_index <- list()
  k <- 1
  num_clusters <- length(unique(sample_order))

  # First loop
  while (k <= num_clusters) {
    clusterSampIndices <- which(k == sample_order)
    i <- 1
    while (i <= length(clusterSampIndices)) {
      if (length(clusterSampIndices) == 1) {
        # Only one sample in cluster
        k_index <- c(k_index, -clusterSampIndices[i])
        break
      } else if (i == 1) {
        # Establish first tree of cluster
        a$merge <-
          rbind(a$merge,
                -c(clusterSampIndices[i], clusterSampIndices[i + 1]))
        a$height[what_row] <- i/length(clusterSampIndices)
        i <- i + 1
        if (is.na(clusterSampIndices[i + 1])) {
          # Last sample of cluster K has been processed
          k_index <- c(k_index, what_row)
        }
      } else {
        a$merge <- rbind(a$merge, c(what_row - 1, -clusterSampIndices[i]))
        a$height[what_row] <- i/length(clusterSampIndices)
        if (is.na(clusterSampIndices[i + 1])) {
          # Last sample of cluster K has been processed
          k_index <- c(k_index, what_row)
        }
      }
      i <- i + 1
      what_row <- what_row + 1
    }
    k <- k + 1
  }

  # Second loop
  for (i in 2:length(k_index)) {
    if (i == 2) {
      a$merge <- rbind(a$merge, c(k_index[[1]], k_index[[2]]))
    } else {
      a$merge <- rbind(a$merge, c(what_row -1, k_index[[i]]))
    }
    a$height[what_row] <- 2
    what_row <- what_row + 1
  }

  a$order  <- order(sample_order)
  a$labels <- sample_order
  class(a) <- "hclust"
  d <- as.dendrogram(a)
  return(d)
}
