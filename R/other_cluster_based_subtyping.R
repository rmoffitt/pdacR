#' @title Cluster Based Subtyping
#'
#' @description Perform clustering, label samples, and re-save data sets
#' @export

other_cluster_based_subtyping <- function(data_set_list = c("summary")){
  pdf("other_cluster_based_subtyping.pdf")
  
    #######################################################################

    dataset <- pdac::Chen_GEO_array
    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,9,"Chen_GEO_array",adj = 0,cex = 3)
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
    summary(dataset$sampInfo)
    sampleset <- which(dataset$sampInfo$specimen_type %in% "primaryPDAC"  )
    dataset$sampInfo$tumor.classifier.training[sampleset] <- TRUE

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("specimen_type is primaryPDAC",
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
                                                                   maxK = 3,
                                                                   reps=200,
                                                                   distance="pearson",
                                                                   clusterAlg="kmdist")[[3]]$consensusTree
    cluster.cut <- data.frame(cuts = (c("basal","classical",NA)[cutree(cluster.result.c,k=4)]),stringsAsFactors = FALSE)
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
    cluster.cut <- data.frame(cuts = (c('basal','classical')[cutree(cluster.result.c,k=7)]),stringsAsFactors = FALSE)
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
    dataset$sampInfo$tumor.classifier.training <- factor(dataset$sampInfo$tumor.classifier.training,
                                                         levels = c(FALSE,TRUE))
    dataset$sampInfo$tumor.classifier.outlier <- factor(dataset$sampInfo$tumor.classifier.outlier,
                                                        levels = c(FALSE,TRUE))
    dataset$sampInfo$cluster.MT <- factor(dataset$sampInfo$cluster.MT, levels = c("basal","classical"))
    dataset$sampInfo$cluster.MT.scaled <- factor(dataset$sampInfo$cluster.MT.scaled, levels = c("basal","classical"))
    Chen_GEO_array_plus <- dataset
    save(file = "data/Chen_GEO_array_plus.RData",list = "Chen_GEO_array_plus")

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("after cluster trimming",sum(!is.na(dataset$sampInfo$cluster.MT))),adj = 0)
    #######################################################################

    dataset <- pdac::Moffitt_S2
    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,9,"Moffitt_S2",adj = 0,cex = 3)
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
    summary(dataset$sampInfo)

    sampleset <- which(dataset$sampInfo$sample_type %in% "PDX")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- TRUE

    sampleset <- which(dataset$sampInfo$sample_type %in% "CAF")
    dataset$sampInfo$tumor.classifier.outlier[sampleset] <- TRUE

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("sample_type is PDX",
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
                                                                   maxK = 3,
                                                                   reps=200,
                                                                   distance="pearson",
                                                                   clusterAlg="kmdist")[[2]]$consensusTree
    cluster.cut <- data.frame(cuts = (c('classical',NA,'basal','basal')[cutree(cluster.result.c,k=4)]),stringsAsFactors = FALSE)
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
    cluster.cut <- data.frame(cuts = (c('classical','basal')[cutree(cluster.result.c,k=3)]),stringsAsFactors = FALSE)
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
    dataset$sampInfo$tumor.classifier.training <- factor(dataset$sampInfo$tumor.classifier.training,
                                                         levels = c(FALSE,TRUE))
    dataset$sampInfo$tumor.classifier.outlier <- factor(dataset$sampInfo$tumor.classifier.outlier,
                                                        levels = c(FALSE,TRUE))
    dataset$sampInfo$cluster.MT <- factor(dataset$sampInfo$cluster.MT, levels = c("basal","classical"))
    dataset$sampInfo$cluster.MT.scaled <- factor(dataset$sampInfo$cluster.MT.scaled, levels = c("basal","classical"))
    Moffitt_S2_plus <- dataset
    save(file = "data/Moffitt_S2_plus.RData",list = "Moffitt_S2_plus")

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("after cluster trimming",sum(!is.na(dataset$sampInfo$cluster.MT))),adj = 0)
    #######################################################################

    dataset <- pdac::Moffitt_GEO_array
    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,9,"Moffitt_GEO_array",adj = 0,cex = 3)
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
    summary(dataset$sampInfo)

    sampleset <- which(dataset$sampInfo$specimen_type %in% "Primary")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- TRUE

    sampleset <- which(dataset$sampInfo$specimen_type %in% "Normal")
    dataset$sampInfo$tumor.classifier.outlier[sampleset] <- TRUE

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("specimen_type is Primary",
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
                                                                   maxK = 3,
                                                                   reps=200,
                                                                   distance="pearson",
                                                                   clusterAlg="kmdist")[[2]]$consensusTree
    cluster.cut <- data.frame(cuts = (c('classical','basal',NA)[cutree(cluster.result.c,k=7)]),stringsAsFactors = FALSE)
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
    cluster.cut <- data.frame(cuts = (c('basal','classical')[cutree(cluster.result.c,k=8)]),stringsAsFactors = FALSE)
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
    dataset$sampInfo$tumor.classifier.training <- factor(dataset$sampInfo$tumor.classifier.training,
                                                         levels = c(FALSE,TRUE))
    dataset$sampInfo$tumor.classifier.outlier <- factor(dataset$sampInfo$tumor.classifier.outlier,
                                                        levels = c(FALSE,TRUE))
    dataset$sampInfo$cluster.MT <- factor(dataset$sampInfo$cluster.MT, levels = c("basal","classical"))
    dataset$sampInfo$cluster.MT.scaled <- factor(dataset$sampInfo$cluster.MT.scaled, levels = c("basal","classical"))
    Moffitt_GEO_array_plus <- dataset
    save(file = "data/Moffitt_GEO_array_plus.RData",list = "Moffitt_GEO_array_plus")

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("after cluster trimming",sum(!is.na(dataset$sampInfo$cluster.MT))),adj = 0)

    #######################################################################

    dataset <- pdac::TCGA_PAAD
    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,9,"TCGA_PAAD",adj = 0,cex = 3)
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
    summary(dataset$sampInfo)

    sampleset <- which(dataset$sampInfo$Decision %in% "whitelist")
    dataset$sampInfo$tumor.classifier.training[sampleset] <- TRUE
    dataset$sampInfo$tumor.classifier.outlier[-sampleset] <- TRUE

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("Decision is whitelist",
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
                                                                   maxK = 3,
                                                                   reps=200,
                                                                   distance="pearson",
                                                                   clusterAlg="kmdist")[[2]]$consensusTree
    cluster.cut <- data.frame(cuts = (c(NA,'classical','basal')[cutree(cluster.result.c,k=7)]),stringsAsFactors = FALSE)
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
    cluster.cut <- data.frame(cuts = (c('basal','classical')[cutree(cluster.result.c,k=10)]),stringsAsFactors = FALSE)
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
    dataset$sampInfo$tumor.classifier.training <- factor(dataset$sampInfo$tumor.classifier.training,
                                                         levels = c(FALSE,TRUE))
    dataset$sampInfo$tumor.classifier.outlier <- factor(dataset$sampInfo$tumor.classifier.outlier,
                                                        levels = c(FALSE,TRUE))
    dataset$sampInfo$cluster.MT <- factor(dataset$sampInfo$cluster.MT, levels = c("basal","classical"))
    dataset$sampInfo$cluster.MT.scaled <- factor(dataset$sampInfo$cluster.MT.scaled, levels = c("basal","classical"))
    TCGA_PAAD_plus <- dataset
    save(file = "data/TCGA_PAAD_plus.RData",list = "TCGA_PAAD_plus")

    p <- plot(x=1:10,y=1:10,ann=FALSE,type="n",xaxt="n",yaxt="n")
    p <- p;text(1,8,paste("after cluster trimming",sum(!is.na(dataset$sampInfo$cluster.MT))),adj = 0)
    ############################################################################

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
