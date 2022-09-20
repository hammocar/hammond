#' Create histograms of all interferents and controls for both sample types
#'
#' @param dtf Data frame containing `interferent` column including controls
#' @param metric_str String indicating the metric of interest. Must be a column in `dtf`
#'
#' @return Grid-faceted histograms with means indicated by red dashed lines
#' @export
plot_metric_histograms <- function(dtf, metric_str){
  metric_sym <- rlang::ensym(metric_str)
  dtf <- 
    dtf |> 
    dplyr::mutate(interferent = forcats::as_factor(interferent),
                  interferent = forcats::fct_relevel(interferent, "control", after = Inf)) 
  
  interferent_means <- 
    dtf |> 
    dplyr::group_by(interferent, sample_type) |> 
    dplyr::summarise(interferent_mean = mean({{metric_sym}}), .groups = "drop")
  
  dtf |> 
    ggplot2::ggplot(aes(x = !!metric_sym)) + 
    ggplot2::geom_histogram(bins = 10) + 
    ggplot2::facet_grid(interferent ~ sample_type) + 
    ggplot2::geom_vline(data = interferent_means, 
                        mapping = aes(xintercept = interferent_mean), 
                        color = "red",
                        linetype = "dashed") +
    ggplot2::theme_bw()
}