#' Parse moffitt tumor classifier from 2015 virtual microdissection paper
#'
#' @import openxlsx
#' @export

parse_moffitt_classifier <- function(){

  Moffitt_classifier_2015 <- pdacR::MoffittTumorClassifier_Tue_Jun_9_2015_cv.ktsp.500.100

  Moffitt_classifier_2019 <- classifs$oct25_equivalent_freeze

  # ------------------------------------------

  save(file = "./data/Moffitt_classifier_2015.RData",
       list = "Moffitt_classifier_2015")

  save(file = "./data/Moffitt_classifier_2019.RData",
       list = "Moffitt_classifier_2019")

  # ------------------------------------------
  classifier_list <- data.frame(labels = c("Moffitt_Classifier_2015",
                                           "PurIST_2019"),
                                variablenames =c("Moffitt_classifier_2015",
                                                 "Moffitt_classifier_2019"))
  classifier_list$colors = list(c("blue","white","orange"),
                                c("blue","white","orange"))

  save(file = "./data/classifier_list.RData",
       list = "classifier_list")

  return(NULL)
}
