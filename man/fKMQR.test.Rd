% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fKMQR.R
\name{fKMQR.test}
\alias{fKMQR.test}
\title{P-value computation for method 'Pearson'}
\usage{
fKMQR.test(Y, X, tau, m, test_num, score = NULL, K = NULL)
}
\arguments{
\item{Y}{n*1 vector, the observed outcome for inference}

\item{X}{n*p matrix, the observed covariates for inference}

\item{tau}{k*1 vector, a grid of quantile levels; e.g., 0.5 for the median, 0.75 for the 3rd quartile; default is c(0.1, 0.25, 0.5, 0.75, 0.9)}

\item{m}{numeric variable, the order of B-spline function; default is 3}

\item{test_num}{a vector, representing the test corresponds to which covariate(s) in X.}

\item{score}{'null' as default}

\item{K}{'null' as default}
}
\value{
p-value for method 'Pearson'
}
\description{
P-value computation for method 'Pearson'
}
\keyword{internal}
