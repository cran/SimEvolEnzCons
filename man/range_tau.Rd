% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/range_tau.R
\name{range_tau}
\alias{range_tau}
\title{Bounds of \emph{tau}}
\usage{
range_tau(E_ini_fun,B_fun)
}
\arguments{
\item{E_ini_fun}{Numeric vector of initial concentrations}

\item{B_fun}{Numeric vector of global co-regulation coefficients. Same length as \code{E_ini_fun}.}
}
\value{
Numeric vector of the inferior and the superior bounds of driving variable \eqn{\tau}
}
\description{
Computes the bounds of the driving variable (or position) \eqn{\tau} on the straight line of relative concentrations such as thereof are between 0 and 1
}
\details{
Relative concentrations is defined in [0,1].
This function \code{range_tau} computes the bounds of \eqn{\tau} such as all relative concentrations are between 0 and 1.

The inferior (resp. superior) bound of \eqn{\tau} corresponds to minimal (resp. maximal) value of \eqn{\tau}
such as all relative concentrations are superior or equal to 0 \bold{and} inferior or equal to 1,
with at least one relative concentration equal to 0 or 1.
}
\examples{
beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
B <- apply(beta,1,sumbis)
E0 <- c(30,30,30)

range_tau(E0,B)

}
\seealso{
See details of function \code{\link{droite_e}} for detailed explanations on \eqn{\tau}.
}
