#' @import gridExtra

#' @import grid

#' @import tibble


#' Plot x vs. y with boxplots outside of each axis
#'
#' @param data  A dataframe with the columns "x" and "y" (numeric) to be plotted.
#' @param x A character string "x" with the column name to be plotted on the x-axis
#' @param y A character string "y" with the column name to be plotted on the y-axis
#' @param ... Optional arguments `color` and `color_labels`. `color` must be a factor, `color_labels` must be a character string of the same length as the number of levels in factor `color`
#'
#' @return a ggplot object
#' @export
#'
#' @examples
# a<-rnorm(100, 20, 5)
# b<-rpois(100, 50)
# c<-rbernoulli(100,.7)
# data<-tibble(a = a,
#              b = b,
#              c = c)
# xy_boxplots(data,
#             x = "a",
#             y = "b",
#             color = c)
#

xy_boxplots <- function(data, x, y, ...) {

  data$x <- unlist(data[, x])
  data$y <- unlist(data[, y])

  args<- list(...)

  z<-theme(strip.background = element_rect(fill = "white", colour = NA),
           strip.text = element_text(colour = "black", face = "bold", size = rel(0.8)))

  main<-ggplot(data,
               aes( x = x, y = y,  color = args$color))+
    geom_point(size = 1.5)+
    labs(x = x, y = y, color = "")+
    scale_color_discrete(labels = args$color_labels)+
    theme(legend.position = "bottom")+
    theme(title = element_text(face = "bold", size = 15),
          rect = element_rect(fill = "white",
                              colour = "black",
                              size = 0.5,
                              linetype = 1),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.border = element_rect(colour = "black", linetype ="solid", fill = NA),
          panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92", size = 0.25))+
    z


  x_boxplot<-ggplot(data, aes( y = x))+
    geom_boxplot(aes(color = args$color), width = 3)+
    theme(title = element_text(face = "bold", size = 12),
          axis.text  = element_text(face = "bold", size = 10),
          legend.text  = element_text(size = 12),
          legend.position = "bottom" ,
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92", size = 0.25))+
    z


  y_boxplot<-ggplot(data, aes( y = y, color = args$color))+
    geom_boxplot(aes(color = args$color), width = 3)+
    theme(title = element_text(face = "bold", size = 12),
          axis.text  = element_text(face = "bold", size = 10),
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92", size = 0.25))+
    z

  plots <-list(main +
                 theme(panel.background = element_rect(fill = "white", colour = "black", linetype ="solid")),
               y_boxplot+
                 theme(panel.background = element_rect(fill = "white", colour = "NA"),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       legend.position = "none") ,
               x_boxplot +
                 theme(panel.background = element_rect(fill = "white", colour = "NA"),
                       axis.title.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank(),
                       legend.position = "none") +
                 coord_flip())

  grobs<-list()
  widths <- list()
  heights<-list()

  for (i in 1:length(plots)){
    grobs[[i]] <- ggplotGrob(plots[[i]])

    widths[[i]] <- grobs[[i]]$widths[2:5]
    heights[[i]] <- grobs[[i]]$heights[c(8,11)]


  }

  maxwidth <- do.call(grid::unit.pmax, widths)


  for (i in 1:length(grobs)){
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }

  grobs[[2]]$heights[c(8,10,9,11)]<-grobs[[1]]$heights[c(8,9,10,11)]

  grid.arrange(grobs = grobs, layout_matrix = rbind(c(3,3,3,3,NA),
                                                    c(1,1,1,1,2),
                                                    c(1,1,1,1,2),
                                                    c(1,1,1,1,2)))

}

# xy_boxplots(tibble(a = c(1,2,3),
#                    b = c(4,5,6),
#                    c = c(0,0,1)),
#              x = "a",
#             y = "b",
#             color = "c")
