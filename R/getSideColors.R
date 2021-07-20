#' create a sidecolor object from a sampInfo and a list of tracks
#'
#' @export
#' @import RColorBrewer
#' @import scales

getSideColors <- function(sampInfo,
                          sampleTracks = names(sampInfo),
                          colorlists = c("gray94","blue","green","yellow","orange","red","black"),
                          displaynames = as.list(sampleTracks),
                          drop.levels = FALSE){

  if(is.list(colorlists)==FALSE){
    # one or no colormaps were provided, so use the same for all tracks
    colorlists <- rep(list(colorlists),length(sampleTracks))
  }

  legend.text = c(" ")
  legend.color = c("white")

  SideColors <- NULL
  for(thisTrack in sampleTracks){
    class.of.track <- class(sampInfo[,thisTrack])

    if(class.of.track %in% "matrix"){
      class.of.track <- class(1234)
    }
    switch(class.of.track,
           #----------------------------------------------------
           character = {
             class <- factor(sampInfo[,thisTrack])
           },
           #----------------------------------------------------
           logical = {
             class <- factor(sampInfo[,thisTrack],levels = c(FALSE,TRUE))
           },
           #----------------------------------------------------
           numeric   = {
             tmp <- sampInfo[,thisTrack]
             class <- factor(x = round(rescale(x = tmp,
                                               to = c(1,101))),
                             levels = 1:101)
             tmp <- range(tmp, na.rm = TRUE, finite = TRUE)
           },
           #----------------------------------------------------
           factor    = {
             class <- sampInfo[,thisTrack]
             if(drop.levels==TRUE){
               class <- droplevels(class)
             }
           }
           #----------------------------------------------------
    )

    classPalette <- colorRampPalette(colorlists[[which(sampleTracks %in% thisTrack)]])(n=length(levels(class)))
    classPalette <- classPalette[1:length(levels(class))]
    SideColors <- as.matrix(cbind(SideColors,classPalette[class]))
    colnames(SideColors)[dim(SideColors)[2]] <- displaynames[which(sampleTracks %in% thisTrack)]

    used_colors <- which(levels(class) %in% class)

    switch(class.of.track,
           #----------------------------------------------------
           logical = {
             legend.text = c(legend.text,displaynames[[which(sampleTracks %in% thisTrack)]],levels(class)[]," ")
             legend.color = c(legend.color,"white",classPalette[],"white")
           },
           #----------------------------------------------------
           character = {
             legend.text = c(legend.text,displaynames[[which(sampleTracks %in% thisTrack)]],levels(class)[used_colors]," ")
             legend.color = c(legend.color,"white",classPalette[used_colors],"white")
           },
           #----------------------------------------------------
           numeric   = {
             legend.text = c(legend.text,displaynames[[which(sampleTracks %in% thisTrack)]],
                             signif(c(tmp[1],
                                      tmp[1]*0.75+tmp[2]*0.25,
                                      tmp[1]*0.5+tmp[2]*0.5,
                                      tmp[1]*0.25+tmp[2]*0.75,
                                      tmp[2]),
                                    digits = 3),
                             " ")
             legend.color = c(legend.color,"white",
                              classPalette[c(1,26,51,76,101)],
                              "white")
           },
           #----------------------------------------------------
           factor    = {
             legend.text = c(legend.text,displaynames[[which(sampleTracks %in% thisTrack)]],levels(class)[used_colors]," ")
             legend.color = c(legend.color,"white",classPalette[used_colors],"white")
           }
           #----------------------------------------------------
    )
  }
  if(is.null(SideColors)){
    SideColors <- as.matrix(rep("white",nrow(sampInfo) )) # length(sampInfo[[1]])
  }
  return(list(SideColors=SideColors,text=legend.text,colors=legend.color))
}
