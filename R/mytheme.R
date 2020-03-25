#' My plot theme
#'
#' @param title_size Size of the title text
#' @param axis_and_legend_title_size Size of the legend title, axis titles, and subtitle text.
#' @param text_size Size of axis and legend text
#'
#' @return a ggplot theme object
#' @export
#'
#' @examples
mytheme <- function(title_size = 22, axis_and_legend_title_size = 15, text_size = 12) {

theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size= title_size, hjust=0),
      # Same size for subtitle, axis titles, and legend title, and strip (for facet wraps)
      plot.subtitle  = element_text(family = "Trebuchet MS", color="#666666", face="bold", size= axis_and_legend_title_size, hjust=0),
      axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size= axis_and_legend_title_size),
      legend.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size= axis_and_legend_title_size),
      strip.text = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=axis_and_legend_title_size),
      # Same size for legend and axis text
      legend.text = element_text(family = "Trebuchet MS", color="#666666", face="bold", size= text_size),
      axis.text = element_text(family = "Trebuchet MS", color="#666666", face="bold", size= text_size),
      # get rid of ugly grey background
      legend.key = element_rect(colour = NA, fill = NA),
      panel.background = element_rect(fill = "white", colour = "#666666"),
      panel.grid.major = element_line(size = 0.0015, linetype = 'solid',colour = "#666666"),
      panel.grid.minor = element_line(size = 0.0005, linetype = 'solid',colour = "#666666"),
      strip.background = element_rect(color = "#666666", fill = "white"),
      legend.key.size = unit(2.5, "cm"))

}
#' @import gridExtra

#' @import grid

#' @import tibble

#' @import ggplot2
NULL
