message("Ensuring appropriate install of all dependencies from both CRAN and BioConductor")
if (!require("plyr")) {
  install.packages("plyr",repos="http://cran.rstudio.com/")
}
if (!require("tidyverse")) {
  install.packages("tidyverse",repos="http://cran.rstudio.com/")
}
if (!require("devtools")) {
  install.packages("devtools",repos="http://cran.rstudio.com/")
}

if (!require("testthat")) {
  install.packages("testthat",repos="http://cran.rstudio.com/")
}
if (!require("reshape2")) {
  install.packages("reshape2",repos="http://cran.rstudio.com/")
}
if (!require("gplots")) {
  install.packages("gplots",repos="http://cran.rstudio.com/")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer",repos="http://cran.rstudio.com/")
}
if (!require("shiny")) {
  install.packages("shiny",repos="http://cran.rstudio.com/")
}
if (!require("scales")) {
  install.packages("scales",repos="http://cran.rstudio.com/")
}
if (!require("ggplot2")) {
  install.packages("ggplot2",repos="http://cran.rstudio.com/")
}
if (!require("openxlsx")) {
  install.packages("openxlsx",repos="http://cran.rstudio.com/")
}
if (!require("Rtsne")) {
  install.packages("Rtsne",repos="http://cran.rstudio.com/")
}
if (!require("GGally")) {
  install.packages("GGally",repos="http://cran.rstudio.com/")
}
if (!require("ggpubr")) {
  install.packages("ggpubr",repos="http://cran.rstudio.com/")
}
if (!require("ggrepel")) {
  install.packages("ggrepel",repos="http://cran.rstudio.com/")
}
if (!require("pROC")) {
  install.packages("pROC",repos="http://cran.rstudio.com/")
}
if (!require("survminer")) {
  install.packages("survminer",repos="http://cran.rstudio.com/")
}
if (!require("nlme")) {
  install.packages("nlme",repos="http://cran.rstudio.com/")
}
if (!require("foreign")) {
  install.packages("foreign",repos="http://cran.rstudio.com/")
}
if (!require("shinyjs")) {
  install.packages("shinyjs",repos="http://cran.rstudio.com/")
}
if (!require("decoderr")) {
  devtools::install_github("laurapeng/decoderr")
}
if (!require("rlang")) {
  install.packages("rlang",repos="http://cran.rstudio.com/")
}
if (!require("stringi")) {
  install.packages("stringi",repos="http://cran.rstudio.com/")
}
if (!require("shinythemes")) {
  install.packages("shinythemes",repos="http://cran.rstudio.com/")
}
if (!require("pdacmolgrad")) {
  devtools::install_github("RemyNicolle/pdacmolgrad")
}
if (!require("Seurat")) {
  install.packages("Seurat")
}
if (!require("mosaic")) {
  install.packages("mosaic")
}
if (!require("Cairo")) {
  install.packages("Cairo")
}

myversion <- as.numeric(R.Version()$major)+0.1*as.numeric(R.Version()$minor)
if( myversion < 3.5){
  source("https://bioconductor.org/biocLite.R")
  biocLite()

  if (!require("Biobase")) {
    biocLite("Biobase")
  }

  if (!require("AnnotationDbi")) {
    biocLite("AnnotationDbi")
  }

  if (!require("org.Hs.eg.db")) {
    biocLite("org.Hs.eg.db")
  }
  if (!require("illuminaHumanv4.db")) {
    biocLite("illuminaHumanv4.db")
  }

  if (!require("ConsensusClusterPlus")) {
    biocLite("ConsensusClusterPlus")
  }
  if (!require("GEOquery")) {
    biocLite("GEOquery")
  }
  if (!require("bioDist")) {
    biocLite("bioDist")
  }
  if (!require("preprocessCore")) {
    biocLite("preprocessCore")
  }
  if (!require("sva")) {
    biocLite("sva")
  }
  if (!require("limma")) {
    biocLite("limma")
  }
} else {
  if (!require("BiocManager")) {
    install.packages("BiocManager",repos="http://cran.rstudio.com/")
  }

  if (!require("Biobase")) {
    BiocManager::install("Biobase")
  }

  if (!require("AnnotationDbi")) {
    BiocManager::install("AnnotationDbi")
  }

  if (!require("org.Hs.eg.db")) {
    BiocManager::install("org.Hs.eg.db")
  }
  if (!require("illuminaHumanv4.db")) {
    BiocManager::install("illuminaHumanv4.db")
  }

  if (!require("ConsensusClusterPlus")) {
    BiocManager::install("ConsensusClusterPlus")
  }
  if (!require("GEOquery")) {
    BiocManager::install("GEOquery")
  }
  if (!require("bioDist")) {
    BiocManager::install("bioDist")
  }
  if (!require("preprocessCore")) {
    BiocManager::install("preprocessCore")
  }
  if (!require("sva")) {
    BiocManager::install("sva")
  }
  if (!require("yaml")) {
    BiocManager::install("yaml")
  }
  if (!require("hgu219.db")) {
    BiocManager::install("hgu219.db")
  }
  if (!require("oligo")) {
    BiocManager::install("oligo")
  }
  if (!require("SingleCellExperiment")) {
    BiocManager::install("SingleCellExperiment")
  }
  if (!require("ggpubr")) {
    BiocManager::install("ggpubr")
  }
  if (!require("DESeq2")) {
    BiocManager::install("DESeq2")
  }
  if (!require("limma")) {
    BiocManager::install("limma")
  }
}
message("pdacR dependency install complete")
