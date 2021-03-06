% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/droites.R
\name{droites}
\alias{droites}
\alias{droite_e}
\alias{droite_E.Reg}
\alias{droite_E.CR}
\alias{droite_tau}
\title{Line of relative enzyme concentrations (co-regulations cases)}
\usage{
droite_e(tau_fun,E_ini_fun,B_fun)

droite_E.Reg(tau_fun,E_ini_fun,B_fun)

droite_E.CR(tau_fun,E_ini_fun,B_fun)

droite_tau(E_fun,E_ini_fun,B_fun)
}
\arguments{
\item{tau_fun}{Numeric value of the position of relative enzyme concentrations on the line}

\item{E_ini_fun}{Numeric vector of initial concentrations}

\item{B_fun}{Numeric vector of global co-regulation coefficients. Same length as \code{E_ini_fun}.}

\item{E_fun}{numeric vector of current concentrations on the line. Same length as \code{E_ini_fun}.}
}
\value{
\code{droite_e} returns a numeric vector of relative concentrations

\code{droite_E.Reg} returns a numeric vector of absolute concentrations

\code{droite_E.CR} returns a numeric vector of absolute concentrations

\code{droite_tau} returns a numeric value giving the position on the line
}
\description{
Computes the position of the point of relative concentrations on the line on which they move on in the cases of co-regulations
}
\details{
In the cases of co-regulations, relative enzymes concentrations evolve along a straight line.
This line is determined by to factors: initial enzyme concentrations and global co-regulation factors.
The driving variable \eqn{\tau} is a parameter indicating the position of the relative enzyme concentrations \emph{e} on this line.


Function \code{droite_e} gives the relative concentrations corresponding to the input position \code{tau_fun}.

Function \code{droite_E.Reg} gives the absolute concentrations corresponding to the input position \code{tau_fun} in case of \bold{regulation only} (constraints abbreviation \code{"RegPos"} or \code{"RegNeg"}).     

Function \code{droite_E.CR} gives the absolute concentrations corresponding to the input position \code{tau_fun} in case of \bold{competition and regulation} (constraints abbreviation \code{"CRPos"} or \code{"CRNeg"}).

Function \code{droite_tau} gives the position \eqn{\tau} corresponding to the input enzyme concentrations \code{E_fun}.


Note that if initial relative concentrations \code{E_ini_fun} is a multiple of \code{1/B_fun}, the line becomes a point and position \eqn{\tau} does not exist.
}
\section{Special results}{

\bold{\emph{Initial point} }

If \code{tau_fun} is equal to 0, \code{droite_e} returns the value of initial relative concentrations, \emph{i.e.} value of \code{E_ini_fun} divided by \code{sum(E_ini_fun)};
\code{droite_E.Reg} and \code{droite_E.CR} returns the value of \code{E_ini_fun}.

If \code{E_fun} is a multiple of \code{E_ini_fun}, function \code{droite_tau} returns \code{0} \emph{(initial point)}.

\bold{\emph{End point} }

If \code{tau_fun} is equal to 1, \code{droite_e} returns the reverse value of \code{B_fun};
\code{droite_E.CR} returns a multiple of the reverse value of \code{B_fun};
\code{droite_E.Reg} returns \code{Inf}.

If \code{E_fun} is a multiple of \code{1/B_fun}, function \code{droite_tau} returns \code{1} \emph{(end point)}.

\bold{\emph{Line becomes a point} }

If \code{E_ini_fun} is a multiple of \code{1/B_fun},\code{droite_e} returns the value of \code{1/B_fun};
\code{droite_E.Reg} returns an error, because concentrations \code{E} can be any multiple of \code{1/B_fun} without variation of relative concentrations;
\code{droite_E.CR} returns \code{E_ini_fun};
\code{droite_tau} returns an error because \eqn{\tau} does not exist in this case.

\bold{\emph{Line does not exist} }

If there only one enzyme (length of \code{E_ini_fun} is equal to 1), relative concentrations is always equal to 1.
\code{droite_e} should return 1;
\code{droite_E.Reg}, \code{droite_E.CR} and \code{droite_tau} should return an error.
}

\examples{
beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
B <- apply(beta,1,sumbis)
E0 <- c(30,30,30)
tau <- 0.5

droite_e(tau,E0,B)

E <- droite_E.Reg(tau,E0,B)
droite_tau(E,E0,B)

E <- droite_E.CR(tau,E0,B)
droite_tau(E,E0,B)




}
\seealso{
To compute global co-regulation coefficients \code{B_fun} from co-regulation matrix \code{beta_fun}, see the example or use function \code{\link{compute.B.from.beta}}.
}
