#' Calculate percentage of targets called negative for all interferents and the control.
#'
#' @param dtf Data frame containing `interferent` column including controls
#' @param sample_type_str String indicating the sample type of interest
#'
#' @return Data frame including total negative calls, total targets, and percent negative calls
#' @export
#'
#' @examples
vl_pct_negative_calls <- function(dtf, sample_type_str){
  dtf |> 
    dplyr::mutate(interferent = forcats::as_factor(interferent),
                  interferent = forcats::fct_relevel(interferent, "control", after = Inf)) |> 
    dplyr::filter(sample_type == sample_type_str) |> 
    dplyr::group_by(interferent) |> 
    dplyr::summarise(total_negative_calls = sum(target_call == "NEGATIVE"),
                     total_targets = n(), 
                     percent_negative_calls = 
                       total_negative_calls / total_targets * 100, 
                     .groups = "drop") 
}

#' Calculate percentage of targets called positive for all interferents and the control.
#'
#' @param dtf Data frame containing `interferent` column including controls
#' @param sample_type_str String indicating the sample type of interest
#'
#' @return Data frame including total positive calls, total targets, and percent positive calls
#' @export
#'
#' @examples
vl_pct_positive_calls <- function(dtf, sample_type_str){
  dtf |> 
    dplyr::mutate(interferent = forcats::as_factor(interferent),
                  interferent = forcats::fct_relevel(interferent, "control", after = Inf)) |> 
    dplyr::filter(sample_type == sample_type_str) |> 
    dplyr::group_by(interferent) |> 
    dplyr::summarise(total_positive_calls = sum(target_call == "POSITIVE"),
                     total_targets = n(), 
                     percent_positive_calls = 
                       total_positive_calls / total_targets * 100, 
                     .groups = "drop") 
}
