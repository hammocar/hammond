#' Calculate target-level breakdown of all replicates of a interferent/sample type combination
#'
#' @param dtf Data frame containing `interferent` column including controls
#' @param interferent_str String indicating the interferent of interest
#' @param sample_type_str String indicating the sample type of interest
#'
#' @return Data frame containing positive calls per target and total number of replicates
#' @export
#'
#' @examples
vl_positive_calls_by_target <- function(dtf, interferent_str, sample_type_str){
  dtf |> 
    dplyr::filter(interferent == interferent_str,
                  sample_type == sample_type_str) |> 
    dplyr::select(interferent, sample_type, replicate_number, target_number, target_call) |> 
    dplyr::group_by(target_number) |> 
    dplyr::summarise(positive_calls = sum(target_call == "POSITIVE"), 
                     n_replicates = n(), 
                     .groups = "drop") |> 
    dplyr::mutate(discrepant_call = dplyr::case_when(
      positive_calls == 0 ~ FALSE,
      positive_calls == n_replicates ~ FALSE,
      TRUE ~ TRUE
    )) 
}

#' Calculate number of positive calls per target in a wide format
#'
#' @param dtf Data frame containing `interferent` column including controls
#' @param interferent_str String indicating the interferent of interest
#' @param sample_type_str String indicating the sample type of interest
#'
#' @return Data frame containing positive calls per target with targets 1-16 as columns 
#' @export
#'
#' @examples
vl_positive_calls_wide <- function(dtf, interferent_str, sample_type_str){
  vl_positive_calls_by_target(dtf, interferent_str, sample_type_str) |> 
    dplyr::select(target_number, positive_calls) |> 
    tidyr::pivot_wider(names_from = target_number, values_from = positive_calls, 
                       names_prefix="target_") |> 
    dplyr::mutate(interferent = interferent_str, 
                  sample_type = sample_type_str) |> 
    dplyr::relocate(interferent, sample_type, .before = 1)
}

# vl_discrepant_calls_total <- function(dtf, interferent_str, sample_type_str){
#   vl_positive_calls_by_target(dtf, interferent_str, sample_type_str) |>  
#     dplyr::pull(discrepant_call) |> 
#     sum() |> 
#     tibble::enframe(name = NULL, value = "total_discrepant_calls") |> 
#     dplyr::mutate(interferent = interferent_str, 
#                   sample_type = sample_type_str) |> 
#     dplyr::relocate(interferent, sample_type, .before = 1)
# }
