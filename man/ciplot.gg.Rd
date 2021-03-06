% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ciplot.gg.R
\name{ciplot.gg}
\alias{ciplot.gg}
\title{Draw 95\% Confidence Interval}
\usage{
ciplot.gg(d, exp = FALSE, yrange = c(-5, 5), group = c("Right", "Left"),
  group.label = c("Right Hemisphere", "Left Hemisphere"))
}
\arguments{
\item{d}{Bootstrapped coefficients matrix.}

\item{exp}{Whether take exponential for the estimated CIs.}

\item{yrange}{set up range for the 95\% CIs.}
}
\description{
Draw 95\% Confidence Interval
}
\examples{
ciplot.gg(d)
}
\references{
TBA
}
\author{
Liwen Wu, Seonjoo Lee, \email{sl3670@cumc.columbia.edu}
}
\keyword{confidence}
\keyword{interval}
