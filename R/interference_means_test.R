#' Perform difference of means test between control replicates and interferent replicates given sample type
#'
#' @param dtf Data frame containing `interferent` column including controls
#' @param interferent_str String indicating the interferent of interest
#' @param metric_str String indicating the metric of interest. Must be a column of `dtf`
#' @param sample_type_str String indicating the sample type of interest
#' @param conf_level Two-sided confidence level for the t.test 
#'
#' @return Tibble containing estimated means for control and interferent replicate samples, the point estimate for difference in means, and the confidence interval for difference in means.
#' @export
#'
#' @examples
test_interferent_means <- 
  function(dtf, interferent_str, metric_str, sample_type_str, conf_level = 0.95){
    metric_sym <- rlang::ensym(metric_str)
    
    control_metric <- 
      dtf |> 
      dplyr::filter(sample_type == sample_type_str,
                    interferent == "control") |> 
      dplyr::pull({{metric_sym}})
    
    interferent_metric <- 
      dtf |> 
      dplyr::filter(sample_type == sample_type_str,
                    interferent == interferent_str) |> 
      dplyr::pull({{metric_sym}})
    
    mdl <- 
      t.test(interferent_metric, 
             control_metric, 
             alternative = "two.sided", 
             conf.level = conf_level)
    
    tibble::tibble(
      interferent = interferent_str, 
      sample_type = sample_type_str, 
      metric = metric_str, 
      mean_interferent = mdl$estimate[1],
      mean_control = mdl$estimate[2],
      mean_difference = mean_interferent - mean_control, 
      conf_level = conf_level, 
      ci_lower = mdl$conf.int[1], 
      ci_upper = mdl$conf.int[2], 
      p_value = mdl$p.value
    )
  }

#' Produce well formated tables of results from `test_interferent_means()`
#'
#' @param dtf Data frame result from `test_interferent_means()`
#'
#' @return Tables of confidence intervals. Must be used in .Rmd file in a chunk with option `results = "asis"`
#' @export
#'
#' @examples
table_interferent_cis <- function(dtf){
  dtf |> 
    dplyr::group_by(sample_type) |> 
    dplyr::group_walk(.f = \(.x, .y){
      metric <- .x$metric[1]
      
      t <- 
        .x |> 
        dplyr::select(-metric) |> 
        kableExtra::kable(digits = 3, booktabs = TRUE, caption = .y[[1]]) |> 
        kableExtra::kable_styling(full_width=FALSE, bootstrap_options="striped") 
      
      print(t)
    }) 
}

#' Produce plots of confidence intervals from `test_interferent_means()`
#'
#' @param dtf 
#'
#' @return
#' @export
#'
#' @examples
plot_interferent_cis <- function(dtf){
  dtf |> 
    dplyr::mutate(interferent = forcats::fct_reorder(interferent, dplyr::desc(interferent))) |> 
    ggplot(aes(xmin = ci_lower, xmax = ci_upper, y = interferent)) + 
    ggplot2::geom_errorbar() + 
    ggplot2::geom_point(aes(x = mean_difference, y = interferent)) + 
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") + 
    ggplot2::facet_grid(. ~ sample_type) + 
    ggplot2::labs(title = dtf$metric) + 
    ggplot2::theme_bw()
}