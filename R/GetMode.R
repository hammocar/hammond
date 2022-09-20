#' Finds the mode(s) of a vector.
#'
#' @param x A vector of values for which you want to find the modes.
#' @return Returns the mode of a vector. If there are multiple modes it returns a vector of the modes.
#' @examples
#' V = c("A","A","B","C")
#' GetMode(V)
#' @export
#' 

GetMode = function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}
