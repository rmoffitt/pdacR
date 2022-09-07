pull_CPTAC_seq <- function(){
  library(TCGAbiolinks)
  library(readxl)
  query = GDCquery(project = 'CPTAC-3',
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "STAR - Counts")

  ## Current CPtAC structure causes error in TCGAbiolinks query workflow. Temporary fix reccomended by TCGAbiolinks is to remove duplicates for now. Once TCGAbiolinks is updated, this will be removed.
  query$results[[1]] <- query$results[[1]][!duplicated(query$results[[1]]$sample.submitter_id),]
  #GDCdownload(query)

  expdat <- GDCprepare(query = query,
                       save = F,
                       summarizedExperiment = F)
  ind = grep("tpm",colnames(expdat),ignore.case = T)
  CPTAC_exp = list(ex = expdat[,..ind],
                   sampInfo = data.frame(),
                   featInfo = data.frame(ENSEMBL = expdat$gene_id,
                                         SYMBOL = expdat$gene_name,
                                         gene_type = expdat$gene_type))

  colnames(CPTAC_exp$ex) = stringr::str_split_fixed(colnames(CPTAC_exp$ex),"_",3)[,3]

  CPTAC_exp$sampInfo = data.frame(samples = colnames(CPTAC_exp$ex),
                                  submitter_id = substr(colnames(CPTAC_exp$ex),1,9))
  CPTAC_exp$sampInfo$submitter_id[-grep("C3",CPTAC_exp$sampInfo$submitter_id)] = CPTAC_exp$sampInfo$samples[-grep("C3",CPTAC_exp$sampInfo$submitter_id)]

  ## Keep only the Pancreas associated files
  clin_query = GDCquery_clinic(project = "CPTAC-3")
  panc_only = clin_query[grep("pancreas", clin_query$tissue_or_organ_of_origin, ignore.case = T,value = F),]


  ## Grab information regarding if normal or tumor
  clin_query = GDCquery_clinic(project = "CPTAC-3", type = "biospecimen")
  clin_query = clin_query[,c("submitter_id","tissue_type","sample_type","composition")]



  ## Add duplicate rows to sampInfo for patients with multiple submissions
  CPTAC_exp$sampInfo = dplyr::right_join(x = CPTAC_exp$sampInfo,
                                         y = panc_only,
                                         by = "submitter_id")

  ## Remove sampInfo for non-Panc associated info (missing in >= 67%)
  CPTAC_exp$sampInfo = CPTAC_exp$sampInfo[,-which(as.numeric(apply(CPTAC_exp$sampInfo,2,FUN = function(x){sum(is.na(x))})) > nrow(CPTAC_exp$sampInfo)*.67)]

  ## Trim to panc only
  indx = which(colnames(CPTAC_exp$ex) %in% CPTAC_exp$sampInfo$samples)
  CPTAC_exp$ex = as.matrix(CPTAC_exp$ex[,..indx])

  ## Order ex cols same as sampInfo rows
  CPTAC_exp$ex = CPTAC_exp$ex[,match(colnames(CPTAC_exp$ex),
                                     table = CPTAC_exp$sampInfo$samples)]

  ## use a loop to identify which entries are tumor vs normal and blood vs tissue
  for(i in 1:nrow(clin_query)){
    subId = clin_query$submitter_id[i]
    cols = grep(subId,x = CPTAC_exp$sampInfo$samples)
    print(cols)
    if(length(cols)>0){
      CPTAC_exp$sampInfo[cols,c("tissue_type","sample_type","composition")] = clin_query[i,2:4]
    }
  }

  ## parse survival info
  CPTAC_exp$sampInfo$censorA.0yes.1no = factor(CPTAC_exp$sampInfo$vital_status,
                                               levels = levels(as.factor(CPTAC_exp$sampInfo$vital_status)),
                                               labels = c(0,1,NA))
  CPTAC_exp$sampInfo$survivalA = numeric(length = nrow(CPTAC_exp$sampInfo))
  CPTAC_exp$sampInfo$survivalA[CPTAC_exp$sampInfo$vital_status == "Alive"] = CPTAC_exp$sampInfo$days_to_last_follow_up[CPTAC_exp$sampInfo$vital_status == "Alive"]
  CPTAC_exp$sampInfo$survivalA[CPTAC_exp$sampInfo$vital_status == "Dead"] = CPTAC_exp$sampInfo$days_to_death[CPTAC_exp$sampInfo$vital_status == "Dead"]

  CPTAC_exp$sampInfo$censorB.0yes.1no = factor(CPTAC_exp$sampInfo$progression_or_recurrence,
                                               levels = levels(as.factor(CPTAC_exp$sampInfo$progression_or_recurrence)),
                                               labels = c(0,NA,1))
  CPTAC_exp$sampInfo$survivalB = numeric(length = nrow(CPTAC_exp$sampInfo))
  CPTAC_exp$sampInfo$survivalB[CPTAC_exp$sampInfo$progression_or_recurrence == "yes"] = CPTAC_exp$sampInfo$days_to_recurrence[CPTAC_exp$sampInfo$progression_or_recurrence == "yes"]
  CPTAC_exp$sampInfo$survivalB[CPTAC_exp$sampInfo$progression_or_recurrence == "no"] = CPTAC_exp$sampInfo$days_to_last_follow_up[CPTAC_exp$sampInfo$progression_or_recurrence == "no"]


  ## Make metadata
  CPTAC_exp$metadata = list(log.transformed = FALSE,
                            survivalA = "overall survival days",
                            survivalB = "recurrence-specific survival days",
                            exp.type = "RNAseq",
                            reference = "Cao L et al, Cell, 2021, PMD: 34534465",
                            description = "213 solid tissue samples from 161 patients. 52 are adjacent normal. Of the 161 tumor samples, 105 were cleared as having sufficient cellularity by KRAS VAF",
                            accession = "GDC: CPTAC-3 PDA")

  informative = apply(CPTAC_exp$sampInfo,2,FUN = function(x){length(unique(x))>1})
  CPTAC_exp$sampInfo = CPTAC_exp$sampInfo[,informative]
  CPTAC_exp$sampInfo = CPTAC_exp$sampInfo[,
                                          c(1:3,5:6,8:10,18:19,22:24,28:30,38:ncol(CPTAC_exp$sampInfo))]

  ## Add supplementary info from publication
  CPTAC_sampleInfo <- read_excel("inst/extdata/CPTAC/CPTAC_sampleInfo.xlsx", sheet = "Molecular_phenotype_data")
  colnames(CPTAC_sampleInfo)[1] = "submitter_id"

  CPTAC_exp$sampInfo = dplyr::right_join(CPTAC_exp$sampInfo, CPTAC_sampleInfo, by = "submitter_id")
  CPTAC_exp$sampInfo[which(CPTAC_exp$sampInfo$tissue_type == "Normal"),29:ncol(CPTAC_exp$sampInfo)] = NA

  ## Make sure continuous variables are represented as numeric
  CPTAC_exp$sampInfo[,3:ncol(CPTAC_exp$sampInfo)] = as.data.frame(lapply(CPTAC_exp$sampInfo[,3:ncol(CPTAC_exp$sampInfo)],
                                                                         FUN = function(x){
                                                                           if(length(levels(as.factor(x))) > 10) {
                                                                             as.numeric(x)
                                                                           } else {
                                                                             as.character(x)
                                                                           }
                                                                         }
  ))

  CPTAC_exp$sampInfo$cellularity_call_from_VAF = factor(ifelse(CPTAC_exp$sampInfo$KRAS_VAF >= 0.075, "Acceptable", "LowPurity"))
  CPTAC_exp$sampInfo = CPTAC_exp$sampInfo[,order(colnames(CPTAC_exp$sampInfo))]
  CPTAC_exp$sampInfo = CPTAC_exp$sampInfo[,c(which(colnames(CPTAC_exp$sampInfo) == "submitter_id"),
                                             which(!(colnames(CPTAC_exp$sampInfo) == "submitter_id")))]
  CPTAC_exp$sampInfo = CPTAC_exp$sampInfo[which(CPTAC_exp$sampInfo$samples %in% colnames(CPTAC_exp$ex)),]

  ## CPTAC has released more seq files, but original supplement covers initial 179 samples. Generate a public version of publication file, save additional for now.
  CPTAC_exp$ex = CPTAC_exp$ex[,which(colnames(CPTAC_exp$ex) %in% CPTAC_exp$sampInfo$samples)]

  ## =============================
  # Generate default filter
  ## =============================
  CPTAC_exp$metadata$default_selections = list(filter_column = c("cellularity_call_from_VAF","sample_type"),
                                               filter_levels = c("cellularity_call_from_VAF:LowPurity", "sample_type:Solid Tissue Normal"))


  saveRDS(CPTAC_exp, file = "./data/CPTAC_exp.rds", compress = T)
  return(NULL)
}
