#' generate_mouse_data_set_list
#'
#' @export

generate_mouse_data_set_list <- function() {
  mouse_data_set_list <- data.frame(labels=c(
    "Moffitt seq, 2015"),
    variablenames =c("Moffitt_S2.Mm"))

  saveRDS(mouse_data_set_list,
          file = "./data/mouse_data_set_list.rds")
  return(NULL)
}
