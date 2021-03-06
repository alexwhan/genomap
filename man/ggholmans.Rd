% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genoplotters.R
\name{ggholmans}
\alias{ggholmans}
\title{Produce a plot of a trinomial distribution in a Holmans' triangle ' ' Takes a
data.frame of three columns and returns a triangle plot with a point for each
row, and positioned between the vertices based on the trinomial distribution.
After broman::triplot}
\usage{
ggholmans(df, colour = NULL, colourLegend = NULL, labelPoints = NULL)
}
\arguments{
\item{df}{A data.frame of three columns. Where a row does not sum to 1, it
will be rescaled. Negative values will cause an error.}

\item{colour}{A vector to be passed to geom_point to colour points. Needs to
be sorted the same as df, since it is joined by /code{cbind()}.}

\item{colourLegend}{A character string to be used as a legend label if points
are coloured.}

\item{labelPoints}{A vector of labels for points, if the points in the plot
should be labelled.}
}
\description{
Produce a plot of a trinomial distribution in a Holmans' triangle ' ' Takes a
data.frame of three columns and returns a triangle plot with a point for each
row, and positioned between the vertices based on the trinomial distribution.
After broman::triplot
}
