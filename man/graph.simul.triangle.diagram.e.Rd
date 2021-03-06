% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.simul.triangle.diagram.R
\name{graph.simul.triangle.diagram.e}
\alias{graph.simul.triangle.diagram.e}
\title{Triangular diagram of relative concentrations for simulations of enzyme evolution}
\usage{
graph.simul.triangle.diagram.e(all_res_sim,which.enz=c(1,2,3),which.sim=NULL,
new.window=FALSE,add.eq=TRUE,add.line.eq.eff=FALSE,add.curve.lines=FALSE,
nbniv=9,niv.palette=NULL,posi.legend="topleft",cex.ini=1,cex.th=1.2,cex.eff=0.6)
}
\arguments{
\item{all_res_sim}{List, the output of function \code{\link{simul.evol.enz.multiple}} (results of evolution simulation).}

\item{which.enz}{Numeric. Which enzymes would be represented? Default is the first three enzymes, i.e. \code{c(1,2,3)}.}

\item{which.sim}{Numeric vector containing integer numbers between 1 and \code{nsim}. Which simulations would you represent? If \code{NULL} (default), all simulations would be represented.}

\item{new.window}{Logical. Do graphics appear in a new window?}

\item{add.eq}{Logical. Do equilibrium appear on graph?}

\item{add.line.eq.eff}{Logical. Add line of effective equilibrium for all initial concentrations? Default is \code{FALSE}. \emph{See details.}}

\item{add.curve.lines}{Logical. Add curve lines for flux dome? Default is \code{FALSE}. \emph{See details.}}

\item{nbniv}{Numeric. Number of contour lines, between 3 and 11. Default is \code{9}}

\item{niv.palette}{Character vector. Color palette to be passed to contour lines.}

\item{posi.legend}{Contour line legend position. See \emph{details}.}

\item{cex.ini, cex.th, cex.eff}{Numeric. Size of remarkable points (respectively initial, theoretical and effective relative concentrations).}
}
\value{
Invisible matrix \code{why_no_lines} who explain why line of effective equilibrium or curve lines not appear even the corresponding options are set to \code{TRUE}.
If \code{why_no_lines} is \code{NULL}, options \code{add.line.eq.eff} and \code{add.curve.lines} have been set to \code{FALSE} (default).

Function \code{graph.simul.triangle.diagram.e} returns point coordinates on triangular diagram.
}
\description{
Graphics of enzyme evolution simulations obtained by function \code{\link{simul.evol.enz.multiple}}.
Function \code{graph.simul.triangl.diagram.e} gives a triangular diagram of relative concentrations.
}
\details{
\bold{WARNING! If there is more than three enzymes in simulations (\eqn{n>3}), be careful for interpretations.}
Indeed, triangular diagram is projection of a plane for which sum of coordinates is equal to 1. 
In other words, if \eqn{n>3}, represented points are not strictly relative concentrations such as \eqn{e_i = E_i / Etot},
but are \eqn{x_i = E_i / (E_i + E_j + E_k)} for the three represented enzymes \code{i}, \code{j} and \code{k}.

Triangular diagram cannot be developed if there is less than three enzymes (\eqn{n<3} or \code{length(which.enz)<3}).
Indeed, it is a three-dimensional graph. 

If there is only one simulation (length of \code{which.sim} equal to 1), color is black. Else, color are taken in palette \code{rainbow}.


Option \code{add.line.eq.eff=TRUE} adds a line of effective equilibrium for all possible initial concentrations.
This line is only available for certain constraints, which are \code{"CRPos"} and \code{"CRNeg"}.
Code is written only for three enzymes \code{n=3}, and requires to use same initial kinetic parameters for all simulations (set \code{same.kin0=TRUE}) and no mutation of activities (set \code{pmutA=0}).
Also, same initial total concentration for all enzyme is required. Easiest way is to set \code{is.random.E0=TRUE}.
See function \code{\link{predict_eff_allE0}}.

Option \code{add.curve.lines=TRUE} adds contour lines of flux dome.
Contour lines are only available in case of competition, which are \code{"Comp"}, \code{"CRPos"} and \code{"CRNeg"}.
It is only available for three enzymes \code{n=3}, and other conditions are the same than option \code{add.line.eq.eff=TRUE}.
\code{nbniv} and \code{niv.palette} adjust number and color of contour lines respectively.
See function \code{\link{flux.dome.projections}}.


In triangular diagram, left axis correspond to first value of \code{which.enz} for enzyme number, bottom axis to its second value and right axis to its third value.
If options \code{add.line.eq.eff=TRUE} or \code{add.curve.lines=TRUE} are chosen and if all tests have been successfully passed, \code{which.enz} is automatically set to \code{c(1,2,3)}.
}
\examples{

 \donttest{
 #Saved simulation
 data(data_sim_RegNeg)
 graph.simul.triangle.diagram.e(data_sim_RegNeg, new.window=TRUE, add.line.eq.eff=TRUE)
 
 
 #add curve lines
 data(data_sim_Comp)
 graph.simul.triangle.diagram.e(data_sim_Comp,new.window=TRUE,add.curve.lines=TRUE)
 
 #all options
 data(data_sim_CRNeg)
 graph.simul.triangle.diagram.e(data_sim_CRNeg,new.window=TRUE,add.curve.lines=TRUE,
 add.line.eq.eff=TRUE)
 }
 
 
 \donttest{
 #New simulation
# case for 3 enzymes
n <- 3
E0 <- c(30,30,30)
kin <- c(1,10,30)
Keq <- c(1,1,1)
nsim <- 2 # 2 simulations
N <- 1000
beta <- diag(1,n)
beta[upper.tri(beta)] <- c(0.32,0.32*(-0.43),-0.43)
#beta_12 = 0.32, beta_13 = beta_12 x beta_23, beta_23 = -0.43
t_beta <- t(beta) #because R fills matrix column by column
beta[lower.tri(beta)] <- 1/t_beta[lower.tri(t_beta)] #beta_ji = 1/beta_ij
if (n==3) {beta[lower.tri(beta)] <- 1/beta[upper.tri(beta)]} #only available if n=3
correl <- "RegNeg"

evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta,npt=250,
is.random.E0=TRUE,same.E0=FALSE)
graph.simul.triangle.diagram.e(evol_sim, new.window=TRUE, add.line.eq.eff=TRUE)

#add curve lines 
nsim <- 3
correl <- "Comp"
evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta,npt=250,
is.random.E0=TRUE,same.E0=FALSE)
graph.simul.triangle.diagram.e(evol_sim,new.window=TRUE,add.curve.lines=TRUE)

#all options
correl <- "CRNeg"
evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta,npt=250,
is.random.E0=TRUE,same.E0=FALSE)
graph.simul.triangle.diagram.e(evol_sim,new.window=TRUE,add.curve.lines=TRUE,add.line.eq.eff=TRUE)

#several enzyme
n <- 5
E0 <- c(30,30,30,30,30)
kin <- c(1,10,30,100,1000)
Keq <- c(1,1,1,1,1)
correl <- "SC"
evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta,npt=250)
graph.simul.triangle.diagram.e(evol_sim,new.window=TRUE)
}


}
\seealso{
See \code{\link{simul.evol.graph.methods}} for other plots of enzyme evolution simulation.

See \code{\link{flux.dome.projections}} to have details about contour lines.
}
