#' Convert trinomial distribution into x y values
#' Takes a data.frame of three columns where each row represents a trinomial
#' distribution. Returns xy coordinates for plotting in a triangle.
#'
#' @param df A data.frame of three columns. Where a row does not sum to 1, it
#'   will be rescaled.
trixy <- function(df) {
  halfheight <- cos(pi/4)/2
  xy <- data.frame(x = (df[[1]] * (-0.5) + df[[3]] * 0.5),
                   y = (df[[1]] * (-halfheight) + df[[2]] * halfheight + df[[3]] * (-halfheight)))
  return(xy)
}

#' Produce a plot of a trinomial distribution in a Holmans' triangle ' ' Takes a
#' data.frame of three columns and returns a triangle plot with a point for each
#' row, and positioned between the vertices based on the trinomial distribution.
#' After broman::triplot
#'
#' @param df A data.frame of three columns. Where a row does not sum to 1, it
#'   will be rescaled. Negative values will cause an error.
#' @param colour A vector to be passed to geom_point to colour points. Needs to
#'   be sorted the same as df, since it is joined by /code{cbind()}.
#' @param colourLegend A character string to be used as a legend label if points
#'   are coloured.
#' @param labelPoints A vector of labels for points, if the points in the plot
#'   should be labelled.
#' @export
#' @importFrom ggplot2 aes_string
ggholman <- function(df, colour = NULL, colourLegend = NULL, labelPoints = NULL) {
  #Check there aren't any negative values
  if(any(df < 0, na.rm = TRUE)) stop("The data cannot contain negative values")
  #check there are three columns
  if(ncol(df) != 3) stop("The df needs to have three columns")
  if(any(rowSums(df) - 1 > 1e-06, na.rm = TRUE)) {
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

  if(!is.null(colour)) points$colour <- colour
  if(!is.null(colourLegend)) {
    names(points)[grepl("colour", names(points))] <- colourLegend
  }

  tri <- ggplot2::ggplot(vert, aes_string("x", "y")) + ggplot2::geom_path() +
    ggplot2::geom_text(data = labelpoint, aes_string(label = "label")) +
    ggplot2::coord_fixed(ratio = 1/cos(pi/6)) + ggholman.theme

  if(!is.null(colour)) tri <- tri + ggplot2::geom_point(data = points, aes_string(colour = colourLegend)) else
    tri <- tri + ggplot2::geom_point(data = points)

  if(!is.null(labelPoints)) tri <- tri +
    ggplot2::geom_text(data = points, aes_string(label = "pointLabel"))
  return(tri)
}

#' A theme for trinomial plots
ggholman.theme <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                        axis.title = ggplot2::element_blank(),
                        axis.text = ggplot2::element_blank(),
                        axis.ticks = ggplot2::element_blank())

#' Plot genotype data for a linkage group
#'
#' @param df Sorted data.frame output from /code{genoComp}
#' @export
#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes_string
genoPlot <- function(df) {
  df <- df[!is.na(df$score),]
  df$score <- factor(df$score)
  df.p <- ggplot2::ggplot(df, aes_string("markerName", "Gen.sort")) +
    ggplot2::geom_tile(aes_string(fill = "score")) +
    ggplot2::labs(y = "Individual", x = "Marker") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90), panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10))

  return(df.p)
}

#' Compares segregation of markers based on agreement of two markers (population parents)
#' @param df A data.frame
#' @param genos Not currently used
#' @param parent1 The first parent (usually female) - an unquoted string
#' @param parent2 The second parent (usually male) - an unquoted string
#' @param missingString A string value defining missing genotype scores
#' @importFrom magrittr %>%
#' @importFrom ggplot2 aes_string
genoSeg <- function(df, genos, parent1, parent2, missingString = "--") {
  sort_var <- "ord"
  df.o <- df %>%
    dplyr::mutate_("ord" = parent1 == parent2) %>%
    dplyr::arrange_(lazyeval::interp(~desc(val), val = as.name(sort_var)),
                    deparse(substitute(parent1)))

  dfg <- df %>%
    tidyr::gather_("Genotype", "score",
                   dplyr::select_vars_(names(df),
                                       names(df),
                                       exclude = "markerName"))

  dfg$marker.sort <- factor(dfg$markerName, levels = rev(df.o$markerName), ordered = TRUE)
  dfg <- dfg[!is.na(dfg$score) & dfg$score != missingString]
  dfg$score <- factor(dfg$score)

  df.p <- ggplot2::ggplot(aes_string("marker.sort", "Genotype")) +
    ggplot2::geom_tile(aes_string(fill = "score")) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
    ggplot2::xlab("Marker")

  return(df.p)
}
