#' classify a data set with a trained TSP classifier
#'
#' @export

create.classif = function(dat, classifier, dec= NULL, drop = F, labels = NULL, fit){
  indmat = matrix(-1, ncol(dat), nrow(classifier$TSPs))
  for(i in 1:nrow(classifier$TSPs)){
    p1 = which(rownames(dat) == classifier$TSPs[i,1])
    p2 = which(rownames(dat) == classifier$TSPs[i,2])
    if(length(p1) == 0 | length(p2) == 0){
      indmat[i,] = 0
      print(paste("skipping",classifier$TSPs[i,]))
    }else{
      indmat[,i] = (dat[p1,] > dat[p2,])^2
    }
  }
  colnames(indmat) = paste("indmat", 1:ncol(indmat), sep = "")
  X=cbind(rep(1, nrow(indmat)), indmat)
  p2 = exp(X%*%c(fit$a0,as.numeric(fit$beta)))/(1+exp(X%*%c(fit$a0,as.numeric(fit$beta))))
  data = indmat
  
  tp = (p2 > .5)^2
  names(tp) = names(p2) = colnames(dat)
  if(ncol(dat) != nrow(data)){
    print(dim(dat))
    print(dim(data))
    print(dim(indmat))
    print(dim(train_sub3))
  }
  indmat <- data.frame(t(indmat))
  names(indmat) <- names(dat)
  rownames(indmat) <- paste(classifier$TSPs[,1],classifier$TSPs[,2],sep =".")
  t = list(class = tp, predprob = p2, dat = dat,classifier = classifier, dec = dec, labels = labels, indmat = indmat)
  return(t)
}
