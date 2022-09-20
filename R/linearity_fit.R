#' Fits WLS regression model for a linearity study as per CLSI EP06 (Edition 2) and returns the data summarized by indended value with mean, sd, fitted values, deviations from linearity, and percent deviations from linearity
#'
#' @param data Specifies the dataset to be used for the analysis.
#' @param intended_col A character string containing the name of the variable for the intended measureand levels.
#' @param measured_col A character string containing the name of the variable for the observed measured values.
#' @param weight_col A character string containing the name of the variable for the weights to be used in WLS regression. If left empty, will use # replicates / fitted variance as weights.
#' @return Returns: a data frame summarized by indended measureand level with mean, sd, fitted values, deviations from linearity, and percent deviations from linearity
#' @export
#' 
linearity_fit <- function(data, intended_col, measured_col, weight_col, weight_var) {
  
  # Change to data frame so the column names play nice  
  data<-as.data.frame(data)
  
  # Assign intended and measured column names (for generalizability)
  data$intended_col<-data[,intended_col]
  data$measured_col<-data[,measured_col]
  data$weight_col<-data[,weight_col]
  data$weight_var<-data[,weight_var]
  
  # Summarize data by intended, preserving weight column 
  linearity_data<-
    data %>% 
    group_by(intended_col, weight_col, weight_var) %>% 
    dplyr::summarise(sd_at_intended = sd(measured_col),
                     n = n(),
                     mean = mean(measured_col))
  

  # Fit the expected vs measured WLS model
  linearity_mod<-lm(mean ~ intended_col, data = linearity_data, weights = weight_col)
  
  # Add fitted values to summarized data
  linearity_data$fitted<-linearity_mod$fitted.values
  
  # Create deviations from linearity
  linearity_data$deviation<-linearity_data$mean - linearity_data$fitted
  
  # Create percent deviations
  linearity_data$percent_deviation<- linearity_data$deviation / linearity_data$fitted
  
  # Output data
  linearity_data
  
}
