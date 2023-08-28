#' split something
#'
#' @param x the thing to be split
#' @param split the things where it is split
#'
#' @description
#' does somehting nice
#'
#'
#' @return returns the splitted things
#' @export
#'
#' @examples strsplit1("asda sdsa", split = " ")

strsplit1 <- function(x, split) {
  strsplit(x, split = split)[[1]]
}
