#' SimEvolEnzCons: Simulation of Enzyme Evolution Under Constraints
#' 
#' @description 
#' Simulate the evolution of enzyme concentrations under constraints in a metabolic pathway,
#' in accordance to the theoretical background developed in Coton \emph{et al.} (2021).
#' 
#' The functions are divided in five sections. Succinct definitions are given below.
#' There is also a section about the correct use of recurrent function parameters.
#' 
#' In version 2.0.0 and more, some functions have been added or modified to take account of regulation groups.
#' Most important modifications are listed in section \emph{Regulation groups}.
#' 
#' 
#' @section Informations about recurrent parameters:
#' 
#' \bold{Integer numbers}
#' 
#' \code{n_fun} (or \code{n}) is the number of enzymes in he pathway.
#' 
#' \code{nsim} is the number of simulations.
#' 
#' \bold{Allowed constraints}
#' 
#' \code{correl_fun} is used to indicate the constraint applied on the system.
#' Possible constraints and their corresponding abbreviation to pass in \code{correl_fun} are listed below:
#' \itemize{
#'    \item \code{"SC"}: independence between all enzymes 
#'    \item \code{"Comp"}: competition for resources 
#'    \item \code{"RegPos"}: positive regulation
#'    \item \code{"RegNeg"}: negative regulation 
#'    \item \code{"CRPos"}: competition plus positive regulation 
#'    \item \code{"CRNeg"}: competition plus negative regulation
#' }
#' \code{\link{name.correl}} helps to find the correct abbreviation, and \code{\link{is.correl.authorized}} verifies the abbreviation.
#' 
#' \bold{Co-regulation coefficients}
#' 
#' For cases with co-regulations (i.e. \code{correl_fun} value is \code{"RegPos"}, \code{"RegNeg"}, \code{"CRPos"} or \code{"CRNeg"}), \code{beta_fun} or \code{B_fun} is obligatory.
#' In other cases (i.e. \code{correl_fun} value is \code{"SC"} or \code{"Comp"}), \code{beta_fun} and \code{B_fun} are ignored, that is why default is \code{NULL}.
#' 
#' \code{beta_fun} is a square matrix of size \code{n*n}, indicating co-regulation coefficients. \code{B_fun} is vector of length \code{n}, indicating global co-regulation coefficients.
#' 
#' \code{\link{compute.beta.from.B}} and \code{\link{compute.B.from.beta}} help to compute \code{beta_fun} and \code{B_fun} from another.
#' 
#' Initial concentrations \code{E_ini_fun} is used in the same way as \code{B_fun}
#' 
# 
#' 
#' 
#' @section Basic functions:
#' 
#' Basic functions used for other functions.
#' 
#'  \itemize{
#'    \item \code{\link{flux}}: computes flux depending on A, E and X
#'    \item \code{\link{activities}}: computes pseudo-activities depending on kinetic parameters and equilibrium constants
#'    \item \code{\link{alpha_ij}}: computes redistribution coefficients
#'    \item \code{\link{coef_rep}}: computes response coefficients
#'    \item \code{\link{coef_sel.continue}}: computes selection coefficient from expression \eqn{s_i = R_i \delta_i/E_i}
#'    \item \code{\link{coef_sel.discrete}}: computes selection coefficient from expression \eqn{s_i = (J^m - J^r)/J^r}
#'    \item \code{\link{compute.delta}}: computes the actual effect \eqn{\delta} of a mutation
#'    \item \code{\link{droite_e}}, \code{\link{droite_E.CR}}, \code{\link{droite_E.Reg}} and \code{\link{droite_tau}}: gives respectively relative concentrations \bold{e}, concentrations \bold{E} (competition + regulation), concentrations \bold{E} (regulation only) and driving variable \emph{tau} when there is regulation 
#'    \item \code{\link{range_delta}}: bounds of \eqn{\delta_i} for mutant concentrations  \emph{E_j} between 0 and \emph{Etot}
#'    \item \code{\link{range_tau}}: bounds of driving variable \emph{tau} for relative concentrations between 0 and 1
#'    \item \code{\link{name.correl}}: gives abbreviation of constraint names
#'  }
#'
#' 
#' @section Equilibrium:
#' 
#' Functions used to find equilibrium of relative concentrations.
#' 
#' \itemize{
#'    \item \code{\link{predict_th}}: computes theoretical equilibrium
#'    \item \code{\link{predict_eff}}: computes effective equilibrium
#'    \item \code{\link{predict_eff_allE0}}: computes effective equilibrium for various initial concentrations
#' }
#' 
#' @section Graphics:
#' 
#' Functions for drawing figures.
#' 
#' \itemize{
#'    \item \code{\link{flux.dome.graph3D}}: gives dome of flux in a 3D-plot
#'    \item \code{\link{flux.dome.projections}}: gives projection on plane of relative concentrations of the flux dome in a triangular plot
#'    \item \code{\link{graph.simul.by.time.by.enz}} and \code{\link{graph.simul.by.time.by.sim}}: illustrations of simulation results. Give different plots with time in x-axis, and color pattern depends on enzymes and simulation numbers respectively
#'    \item \code{\link{graph.simul.by.time.RNV}}: various plots of RNVs, computing from simulation results
#'    \item \code{\link{graph.simul.triangle.diagram.e}}: gives a triangular diagram of relative concentrations from simulation results
#'    \item \code{\link{RNV.graph.double.at.eq}}: RNV size and relative concentrations depending on activities at equilibrium
#' } 
#' 
#' @section Simulation:
#' Functions for simulations of enzyme evolution.
#' 
#' \itemize{
#'    \item \code{\link{simul.evol.enz.one}}: evolution of enzyme concentrations (one simulation)
#'    \item \code{\link{simul.evol.enz.multiple}}: evolution of enzyme concentrations (multiple simulations)
#' }
#' 
#' 
#' @section RNV:
#' Functions for computing \emph{Range of Neutral Variation} of enzyme concentrations.
#' 
#' \itemize{
#'    \item \code{\link{RNV.delta.all.enz}}: computes \eqn{\delta_i} at limits of neutral zone for each enzyme
#'    \item \code{\link{RNV.for.simul}}: computes RNV in the simulations
#'    \item \code{\link{RNV.mean.simul}}: computes mean of RNV size in the simulations
#'    \item \code{\link{RNV.ranking.order.factor}}: gives the name and the value of of the factor that influences the ranking order of RNV
#'    \item \code{\link{RNV.size.at.equilibr}}: computes RNV at equilibrium, then plots RNV size against the ranking-order factor
#'    \item \code{\link{graph.simul.by.time.RNV}}: various plots of RNVs, computing from simulation results
#' }
#' 
#' \bold{Informations about \emph{Range of Neutral Variations} (RNV)}
#' 
#' The \emph{Range of Neutral Variations} (RNV) are mutant concentration values such as coefficient selection is between \eqn{1/(2N)} and \eqn{-1/(2N)}.
#' 
#' Inferior (resp. superior) bound of RNV corresponds to selection coefficient equal to \eqn{-1/(2N)} (resp. \eqn{1/(2N)}).
#' 
#' Depending on applied constraint \code{correl_fun}, it exists 1 or 2 RNV.
#' In case of independence (\code{"SC"}) or positive regulation between all enzymes (\code{"RegPos"}), flux has no limit and there is only one RNV.
#' In other cases (competition and/or negative regulation), because flux can reach a maximum, there is two RNV: 
#' a "near" one, for small mutations, and a "far" one for big mutations that put mutants on the other side of flux dome.
#'
#' The RNV size is the absolute value of \eqn{\delta_i^sup} minus \eqn{\delta_i^inf}.
#' If there is no superior bounds but there is two RNVs, RNV size is obtained by the difference of the two \eqn{\delta_i^inf}.
#' 
#' 
#' @section Regulation groups:
#' In version 2.0.0 and more, some functions have been added or modified to take account of regulation groups.
#' Most important new functions are listed below, with function section in parenthesis:
#' 
#' \itemize{
#'    \item \code{\link{class_group}} (basic): classifies enzymes in regulation groups from the co-regulation matrix
#'    \item \code{\link{group_types}} (basic): gives types of regulation groups, aka negative group, positive group or singleton, from co-regulation matrix
#'    \item \code{\link{predict_grp}} (equilibrium): computes equilibrium for the different relative concentrations when there are regulation groups
#'    \item \code{\link{apparent.activities.Aq}} (equilibrium): computes apparent activities at equilibrium for regulation groups
#'    \item \code{\link{graph.simul.group}} (graphics): gives graphics of the different kind of relative enzyme concentrations through time when there are regulation groups
#'    \item \code{\link{extract.tabEtot}} (simulation): extracts table of enzyme concentration from simulation results and computes sum of concentrations in a group. Useful to compute the different kind of relative concentrations.
#' }
#' 
#' Other graphics function have been modified to take account of equilibrium for regulation groups.
#' 
#' Two new datasets have been added, to have examples when there are regulation groups.
#' 
#' However, \code{\link{RNV.size.at.equilibr}} and \code{\link{RNV.ranking.order.factor}} have not been modified, as we cannot determine the relative RNV at equilibrium where there are regulation groups.
#' Moreover, simulations \code{\link{simul.evol.enz.one}} do not compute equilibrium when there are regulation groups.
#' 
#' 
#'  
#' @references
#' 
#' Coton, C., Talbot, G., Le Louarn, M., Dillmann, C., de Vienne, D., 2021. Evolution of enzyme levels in metabolic pathways: A theoretical approach. bioRxiv 2021.05.04.442631. https://doi.org/10.1101/2021.05.04.442631
#' 
#' Version 2: Coton, C., Dillmann, C., de Vienne, D., 2021. Evolution of enzyme levels in metabolic pathways: A theoretical approach. Part 2.
#'
#'
#' @name SimEvolEnzCons
#' @aliases SimEvolEnzCons-package
#' @docType package
#' 
NULL

#>NULL