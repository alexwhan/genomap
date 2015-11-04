#' Convert trinomial distribution into x y values
#'
#' Takes a data.frame of three columns where each row represents a trinomial
#' distribution. Returns xy coordinates for plotting in a triangle.\
#'
#' @param df A data.frame of three columns. Where a row does not sum to 1, it
#'   will be rescaled.
trixy <- function(df) {
  halfheight <- cos(pi/4)/2
  xy <- data.frame(x = (df[,1] * (-0.5) + df[,3] * 0.5),
                   y = (df[,1] * (-halfheight) + df[,2] * halfheight + df[,3] * (-halfheight)))
  return(xy)
}

#' Produce a plot of a trinomial distribution in a Holmans' triangle ' ' Takes a
#' data.frame of three columns and returns a triangle plot with a point for each
#' row, and positioned between the vertices based on the trinomial distribution.
#' After broman::triplot
#'
#' @param df A data.frame of three columns. Where a row does not sum to 1, it
#'   will be rescaled.
#' @param colour A vector to be passed to geom_point to colour points. Needs to
#'   be sorted the same as df, since it is joined by /code{cbind()}.
#' @param labelPoints A vector of labels for points, if the points in the plot
#'   should be labelled.
#' @export
ggholman <- function(df, colour = NULL, labelPoints = NULL) {
  #check there are three columns
  if(ncol(df) != 3) stop("The df needs to have three columns")
  if(any(rowSums(df) - 1 > 1e-06)) {
    df <- df/rowSums(df)
    warning("Rows do not sum to 1 - rescaling")
  }
  halfheight <- cos(pi/4)/2

  vert <- trixy(data.frame(a=c(1, 0, 0, 1),
                           b=c(0, 1, 0, 0),
                           c=c(0, 0, 1, 0)))
  labelpoint <- vert[1:3,] * 1.1
  labelpoint$label <- names(df)
  points <- trixy(df)
  if(!is.null(labelPoints)) points$pointLabel <- labelPoints

  if(!is.null(colour)) points <- cbind(points, colour)
  names(points)[grepl("colour", names(points))] <- "P.value"

  tri <- ggplot(vert, aes(x, y)) + geom_path() + geom_text(data = labelpoint, aes(label = label)) +
    coord_fixed(ratio = 1/cos(pi/6)) + ggholman.theme

  if(!is.null(colour)) tri <- tri + geom_point(data = points, aes(colour = colour)) else
    tri <- tri + geom_point(data = points)

  if(!is.null(labelPoints)) tri <- tri + geom_text(data = points, aes(label = pointLabel))
  return(tri)
}

#' A theme for trinomial plots
ggholman.theme <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                        axis.title = ggplot2::element_blank(),
                        axis.text = ggplot2::element_blank(),
                        axis.ticks = ggplot2::element_blank())
