
#' round_df() rounds the numeric columns of the input data frame to the specified number of digits.
#'
#' @param data data frame to round
#' @param digits Number of digits to round to

#' @export
round_df <- function(data, digits) {
  num_data<-data %>% dplyr::select(where(is.numeric))
  not_num_data<-data %>% dplyr::select(-where(is.numeric))
  
  cbind(not_num_data, round(num_data, digits = digits))
}