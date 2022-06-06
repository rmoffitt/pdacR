#' generate_data_set_list
#'
#' @export

generate_data_set_list <- function() {
  data_set_list <- data.frame(labels=c(
    "Chen array, 2015",
    "CPTAC3-PDA, 2022",
    "ICGC ILLUMINA arrays 2013",
    "ICGC PACA-AU array, 2016",
    "ICGC PACA-AU seq, 2016",
    "ICGC PACA-CA seq, 2016",
    "Moffitt array, 2015",
    "Moffitt seq, 2015",
    "Olive seq, 2019",
    "Puleo array, 2018",
    "scAtlas Pseudobulked",
    "Seino array, 2018",
    "TCGA PAAD, 2017"),
    variablenames =c("Chen_GEO_array",
                     "CPTAC_exp",
                     "Nones_GEO_array",
                     "PACA_AU_array",
                     "PACA_AU_seq",
                     "PACA_CA_seq",
                     "Moffitt_GEO_array",
                     "Moffitt_S2.Hs",
                     "Olive_2019",
                     "Puleo_array",
                     "scAtlas.pseudobulked",
                     "Seino_GEO_array",
                     "TCGA_PAAD"))

  save(list = c("data_set_list"),
       file = "./data/data_set_list.RData")

  # for(i in data_set_list$variablenames){
  #   writeLines('----------------------------')
  #   writeLines(i)
  #   dataset <- get(as.character(i))
  #   thesecols <- grep(names(dataset$sampInfo),pattern = "cluster.MT|tumor.classifier")
  #   print(summary(dataset$sampInfo[,thesecols]))
  #   if("MoffittTumor" %in% names(dataset$sampInfo)){
  #     thesecols <- grep(names(dataset$sampInfo),pattern = "MoffittTumor")
  #     print(table(dataset$sampInfo[,thesecols]))
  #   }
  #   else{
  #     thesecols <- grep(names(dataset$sampInfo),pattern = "cluster.MT")
  #     print(table(dataset$sampInfo[,thesecols]))}
  #   thesecols <- grep(names(dataset$sampInfo),pattern = "tumor.classifier")
  #   print(summary(dataset$sampInfo[,thesecols]))
  # }
  return(NULL)
}
