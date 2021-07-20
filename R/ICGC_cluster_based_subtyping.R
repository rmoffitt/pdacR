#' @title Cluster Based Subtyping
#' @description Perform clustering, label samples, and re-save data sets
#' @export

ICGC_cluster_based_subtyping <- function(data_set_list = c("summary")){
  pdf("ICGC_cluster_based_subtyping.pdf")
  #######################################################################
    dataset <- pdac::PACA_AU_seq
    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,9,"PACA_AU_seq",adj = 0,cex = 3)
    p <- p;text(1,5,paste("Total Samples",length(dataset$sampInfo[[1]])),adj = 0)

    # -------------------------------------------------------------
    # Boilerplate
    plot.hist(dataset)
    normalized.expression <- preprocessCore::normalize.quantiles(as.matrix(log2(1+dataset$ex)))
    dataset$sampInfo$tumor.classifier.training <- FALSE
    dataset$sampInfo$tumor.classifier.outlier <- FALSE
    dataset$sampInfo[["cluster.MT"]] <- as.character(NA)
    dataset$sampInfo[["cluster.MT.scaled"]] <- as.character(NA)

    # -------------------------------------------------------------
    # Select cases to cluster based on metadata

    sampleset <- which(dataset$sampInfo$HistoSubtype %in% "Pancreatic Ductal Adenocarcinoma")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- TRUE
    sampleset <- which(dataset$sampInfo$HistoSubtype %in% "PDA - Adenosquamous carcinoma")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- TRUE

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("HistoSubtype is \n Pancreatic Ductal Adenocarcinoma or \n PDA - Adenosquamous carcinoma",
                          sum(dataset$sampInfo$tumor.classifier.training)),adj = 0)


    sampleset <- which(dataset$sampInfo$Sample.type %in% "Cell line ")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- FALSE

    p <- p;text(1,2,paste("Sample.type is not Cell line",
                          sum(dataset$sampInfo$tumor.classifier.training)),adj = 0)

    sampleset <- which(dataset$sampInfo$HistoSubtype %in% "PDA -  signet ring")
    dataset$sampInfo$tumor.classifier.outlier[sampleset] <- TRUE
    sampleset <- which(dataset$sampInfo$HistoSubtype %in% "Acinar Cell Carcinoma")
    dataset$sampInfo$tumor.classifier.outlier[sampleset] <- TRUE


    # -------------------------------------------------------------
    # Boilerplate
    sampleset <- which(dataset$sampInfo$tumor.classifier.training)
    featureset <- which(dataset$featInfo$SYMBOL %in%
                          c(as.character(pdac::gene_lists$Moffitt.Basal.25),
                            as.character(pdac::gene_lists$Moffitt.Classical.25)))
    smallx <- (normalized.expression[featureset,sampleset])
    smallx.scaled <- t(scale(t(normalized.expression[featureset,sampleset])))

    # -------------------------------------------------------------
    # unscaled clustering
    cluster.result.c <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx),
                                                             seed = 1234,
                                                             pFeature = 0.8,
                                                             pItem = 0.8,
                                                             maxK = 3,
                                                             reps=200,
                                                             distance="pearson",
                                                             clusterAlg="kmdist")[[2]]$consensusTree
    cluster.cut <- data.frame(cuts = (c('classical','basal',NA,NA)[cutree(cluster.result.c,k=6)]),stringsAsFactors = FALSE)
    dataset$sampInfo$cluster.MT[sampleset] <- cluster.cut[[1]]
    visualize.it(dataset,smallx,cluster.result.c,featureset,cluster.cut)

    # -------------------------------------------------------------
    # scaled clustering
    cluster.result.c <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx.scaled),
                                                                   seed = 1234,
                                                                   pFeature = 0.8,
                                                                   pItem = 0.8,
                                                                   maxK = 3,
                                                                   reps=200,
                                                                   distance="pearson",
                                                                   clusterAlg="kmdist")[[2]]$consensusTree
    cluster.cut <- data.frame(cuts = (c('classical','basal',NA)[cutree(cluster.result.c,k=7)]),stringsAsFactors = FALSE)
    dataset$sampInfo$cluster.MT.scaled[sampleset] <- cluster.cut[[1]]
    visualize.it(dataset,smallx,cluster.result.c,featureset,cluster.cut)

    # -------------------------------------------------------------
    # adjustments to training set
    dataset$sampInfo$tumor.classifier.outlier[is.na(dataset$sampInfo$cluster.MT) &
                                                is.na(dataset$sampInfo$cluster.MT.scaled) &
                                                dataset$sampInfo$tumor.classifier.training  ] <- TRUE
    dataset$sampInfo$tumor.classifier.training[is.na(dataset$sampInfo$cluster.MT)] <- FALSE

    # -------------------------------------------------------------
    # Wrapup
    dataset$sampInfo$cluster.MT <- factor(dataset$sampInfo$cluster.MT, levels = c("basal","classical"))
    dataset$sampInfo$cluster.MT.scaled <- factor(dataset$sampInfo$cluster.MT.scaled, levels = c("basal","classical"))
    dataset$sampInfo$tumor.classifier.training <- factor(dataset$sampInfo$tumor.classifier.training,
                                                         levels = c(FALSE,TRUE))
    dataset$sampInfo$tumor.classifier.outlier <- factor(dataset$sampInfo$tumor.classifier.outlier,
                                                        levels = c(FALSE,TRUE))
    PACA_AU_seq_plus <- dataset
    save(file = "data/PACA_AU_seq_plus.RData",list = "PACA_AU_seq_plus")

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("after cluster trimming",sum(!is.na(dataset$sampInfo$cluster.MT))),adj = 0)
  #######################################################################
    dataset <- pdac::PACA_CA_seq
    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,9,"PACA_CA_seq",adj = 0,cex = 3)
    p <- p;text(1,5,paste("Total Samples",length(dataset$sampInfo[[1]])),adj = 0)

    # -------------------------------------------------------------
    # Boilerplate
    plot.hist(dataset)
    normalized.expression <- preprocessCore::normalize.quantiles(as.matrix(log2(1+dataset$ex)))
    dataset$sampInfo$tumor.classifier.training <- FALSE
    dataset$sampInfo$tumor.classifier.outlier <- FALSE
    dataset$sampInfo[["cluster.MT"]] <- as.character(NA)
    dataset$sampInfo[["cluster.MT.scaled"]] <- as.character(NA)

    # -------------------------------------------------------------
    # Select cases to cluster based on metadata

    sampleset <- which(dataset$sampInfo$actual_type %in% "Cell line")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- TRUE

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("actual_type is Cell line",
                          sum(dataset$sampInfo$tumor.classifier.training)),adj = 0)
    # -------------------------------------------------------------
    # Boilerplate
    sampleset <- which(dataset$sampInfo$tumor.classifier.training)
    featureset <- which(dataset$featInfo$SYMBOL %in%
                          c(as.character(pdac::gene_lists$Moffitt.Basal.25),
                            as.character(pdac::gene_lists$Moffitt.Classical.25)))
    smallx <- (normalized.expression[featureset,sampleset])
    smallx.scaled <- t(scale(t(normalized.expression[featureset,sampleset])))

    # -------------------------------------------------------------
    # unscaled clustering
    cluster.result.c <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx),
                                                                   seed = 1234,
                                                                   pFeature = 0.8,
                                                                   pItem = 0.8,
                                                                   maxK = 5,
                                                                   reps=200,
                                                                   distance="pearson",
                                                                   clusterAlg="kmdist")[[2]]$consensusTree
    cluster.cut <- data.frame(cuts = (c(NA)[cutree(cluster.result.c,k=2)]),stringsAsFactors = FALSE)
    dataset$sampInfo$cluster.MT[sampleset] <- cluster.cut[[1]]
    visualize.it(dataset,smallx,cluster.result.c,featureset,cluster.cut)

    # -------------------------------------------------------------
    # scaled clustering
    cluster.result.c <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx.scaled),
                                                                   seed = 1234,
                                                                   pFeature = 0.8,
                                                                   pItem = 0.8,
                                                                   maxK = 3,
                                                                   reps=200,
                                                                   distance="pearson",
                                                                   clusterAlg="kmdist")[[2]]$consensusTree
    cluster.cut <- data.frame(cuts = (c(NA)[cutree(cluster.result.c,k=2)]),stringsAsFactors = FALSE)
    dataset$sampInfo$cluster.MT.scaled[sampleset] <- cluster.cut[[1]]
    visualize.it(dataset,smallx,cluster.result.c,featureset,cluster.cut)

    # -------------------------------------------------------------
    # adjustments to training set

    dataset$sampInfo$tumor.classifier.training[is.na(dataset$sampInfo$cluster.MT)] <- FALSE

    # -------------------------------------------------------------
    # Wrapup
    dataset$sampInfo$cluster.MT <- factor(dataset$sampInfo$cluster.MT, levels = c("basal","classical"))
    dataset$sampInfo$cluster.MT.scaled <- factor(dataset$sampInfo$cluster.MT.scaled, levels = c("basal","classical"))
    dataset$sampInfo$tumor.classifier.training <- factor(dataset$sampInfo$tumor.classifier.training,
                                                         levels = c(FALSE,TRUE))
    dataset$sampInfo$tumor.classifier.outlier <- factor(dataset$sampInfo$tumor.classifier.outlier,
                                                        levels = c(FALSE,TRUE))
    PACA_CA_seq_plus <- dataset
    save(file = "data/PACA_CA_seq_plus.RData",list = "PACA_CA_seq_plus")

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("after cluster trimming",sum(!is.na(dataset$sampInfo$cluster.MT))),adj = 0)
  #######################################################################
    dataset <- pdac::PACA_AU_array
    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,9,"PACA_AU_array",adj = 0,cex = 3)
    p <- p;text(1,5,paste("Total Samples",length(dataset$sampInfo[[1]])),adj = 0)

    # -------------------------------------------------------------
    # Boilerplate
    plot.hist(dataset)
    normalized.expression <- preprocessCore::normalize.quantiles(as.matrix((dataset$ex)))
    dataset$sampInfo$tumor.classifier.training <- FALSE
    dataset$sampInfo$tumor.classifier.outlier <- FALSE
    dataset$sampInfo[["cluster.MT"]] <- as.character(NA)
    dataset$sampInfo[["cluster.MT.scaled"]] <- as.character(NA)

    # -------------------------------------------------------------
    # Select cases to cluster based on metadata

    sampleset <- which(dataset$sampInfo$HistoSubtype %in% "Pancreatic Ductal Adenocarcinoma")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- TRUE
    sampleset <- which(dataset$sampInfo$HistoSubtype %in% "PDA - Adenosquamous carcinoma")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- TRUE

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("HistoSubtype is \n Pancreatic Ductal Adenocarcinoma or \n PDA - Adenosquamous carcinoma",
                          sum(dataset$sampInfo$tumor.classifier.training)),adj = 0)

    sampleset <- which(dataset$sampInfo$Sample.type %in% "Cell line ")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- FALSE

    p <- p;text(1,2,paste("Sample.type is not Cell line",
                          sum(dataset$sampInfo$tumor.classifier.training)),adj = 0)

    sampleset <- which(dataset$sampInfo$HistoSubtype %in% "PDA -  signet ring")
    dataset$sampInfo$tumor.classifier.outlier[sampleset] <- TRUE
    sampleset <- which(dataset$sampInfo$HistoSubtype %in% "Acinar Cell Carcinoma")
    dataset$sampInfo$tumor.classifier.outlier[sampleset] <- TRUE

    # -------------------------------------------------------------
    # Boilerplate
    sampleset <- which(dataset$sampInfo$tumor.classifier.training)
    featureset <- which(dataset$featInfo$SYMBOL %in%
                          c(as.character(pdac::gene_lists$Moffitt.Basal.25),
                            as.character(pdac::gene_lists$Moffitt.Classical.25)))
    smallx <- (normalized.expression[featureset,sampleset])
    smallx.scaled <- t(scale(t(normalized.expression[featureset,sampleset])))

    # -------------------------------------------------------------
    # unscaled clustering
    cluster.result.c <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx),
                                                                   seed = 1234,
                                                                   pFeature = 0.8,
                                                                   pItem = 0.8,
                                                                   maxK = 3,
                                                                   reps=200,
                                                                   distance="pearson",
                                                                   clusterAlg="kmdist")[[2]]$consensusTree
    cluster.cut <- data.frame(cuts = (c('classical','basal',NA,NA)[cutree(cluster.result.c,k=4)]),stringsAsFactors = FALSE)
    dataset$sampInfo$cluster.MT[sampleset] <- cluster.cut[[1]]
    visualize.it(dataset,smallx,cluster.result.c,featureset,cluster.cut)

    # -------------------------------------------------------------
    # scaled clustering
    cluster.result.c <- ConsensusClusterPlus::ConsensusClusterPlus(d = as.matrix(smallx.scaled),
                                                                   seed = 1234,
                                                                   pFeature = 0.8,
                                                                   pItem = 0.8,
                                                                   maxK = 3,
                                                                   reps=200,
                                                                   distance="pearson",
                                                                   clusterAlg="kmdist")[[2]]$consensusTree
    cluster.cut <- data.frame(cuts = (c('classical','basal',NA)[cutree(cluster.result.c,k=9)]),stringsAsFactors = FALSE)
    dataset$sampInfo$cluster.MT.scaled[sampleset] <- cluster.cut[[1]]
    visualize.it(dataset,smallx,cluster.result.c,featureset,cluster.cut)

    # -------------------------------------------------------------
    # adjustments to training set
    dataset$sampInfo$tumor.classifier.outlier[is.na(dataset$sampInfo$cluster.MT) &
                                                is.na(dataset$sampInfo$cluster.MT.scaled) &
                                                dataset$sampInfo$tumor.classifier.training  ] <- TRUE
    dataset$sampInfo$tumor.classifier.training[is.na(dataset$sampInfo$cluster.MT)] <- FALSE

    # -------------------------------------------------------------
    # Wrapup
    dataset$sampInfo$cluster.MT <- factor(dataset$sampInfo$cluster.MT, levels = c("basal","classical"))
    dataset$sampInfo$cluster.MT.scaled <- factor(dataset$sampInfo$cluster.MT.scaled, levels = c("basal","classical"))
    dataset$sampInfo$tumor.classifier.training <- factor(dataset$sampInfo$tumor.classifier.training,
                                                         levels = c(FALSE,TRUE))
    dataset$sampInfo$tumor.classifier.outlier <- factor(dataset$sampInfo$tumor.classifier.outlier,
                                                        levels = c(FALSE,TRUE))

    PACA_AU_array_plus <- dataset
    save(file = "data/PACA_AU_array_plus.RData",list = "PACA_AU_array_plus")
    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("after cluster trimming",sum(!is.na(dataset$sampInfo$cluster.MT))),adj = 0)
  #######################################################################
    dev.off()
}

visualize.it <- function(dataset,smallx,cluster.result.c,featureset,cluster.cut){
  cluster.result.r <- hclust(d = bioDist::cor.dist(x = (smallx)))
  heatmap.3(x = smallx, scale="row",labRow = dataset$featInfo$SYMBOL[featureset],
            col = colorRampPalette(c("blue", "white", "red"))(n = 299),
            Colv = as.dendrogram(cluster.result.c),
            Rowv = as.dendrogram(cluster.result.r),
            ColSideColors =  getSideColors(sampInfo = data.frame(cuts =
                                                                   factor(x = cluster.cut$cuts,
                                                                          levels = c("basal","classical","activated","normal"))),
                                           sampleTracks = "cuts",
                                           colorlists = list(c("orange","blue","brown","skyblue")),
                                           drop.levels = FALSE)$SideColors,
            RowSideColors =  t(getSideColors(sampInfo = data.frame(basal =dataset$featInfo$SYMBOL[featureset] %in%
                                                                     pdac::gene_lists$Moffitt.Basal.25,
                                                                   classical =dataset$featInfo$SYMBOL[featureset] %in%
                                                                     pdac::gene_lists$Moffitt.Classical.25,
                                                                   normal =dataset$featInfo$SYMBOL[featureset] %in%
                                                                     pdac::gene_lists$Moffitt.Normal.25,
                                                                   activated = dataset$featInfo$SYMBOL[featureset] %in%
                                                                     pdac::gene_lists$Moffitt.Activated.25),
                                             sampleTracks = c("basal",
                                                              "classical",
                                                              "normal",
                                                              "activated"),
                                             colorlists = list(b=c("white","orange"),
                                                               c=c("white","blue"),
                                                               n=c("white","skyblue"),
                                                               a=c("white","brown")))$SideColors))
}

plot.hist <- function(dataset){
  par(mfrow=c(2,1))
  h <- hist(log2(1+as.matrix(dataset$ex)),breaks=100,plot=FALSE)
  plot(y=sqrt(h$count),x=h$mids,type='h',main = "Logged")
  h <- hist(as.matrix(dataset$ex),breaks=100,plot=FALSE)
  plot(y=sqrt(h$count),x=h$mids,type='h',main = "As Is")
  }
