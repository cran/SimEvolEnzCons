#' Indirect mutation of enzyme concentrations
#'
#' Computes the mutant value of enzyme concentrations by an indirect method.
#' 
#'
#' @details 
#' This mutation method is named \emph{indirect}, because redistribution coefficient matrix \eqn{M_\alpha} and actual mutation effect are used to compute mutant values
#' rather than canonical mutation effect. Expression is : \eqn{E_j^m = E_j^r + \alpha_ij * \delta_i}
#' 
#' Constraints between enzymes are implicitly described in redistribution coefficients matrix. 
#' 
#' 
#' @usage mut.E.indirect(delta_fun,E_res,alpha_fun,i_fun)
#'
#' @param delta_fun Numeric value of the actual effect of a mutation targeting enzyme \code{i_fun}, i.e. \eqn{\delta_i}
#' @param E_res Numeric vector of resident enzyme concentrations
#' @param alpha_fun Numeric matrix of redistribution coefficients
#' @param i_fun Integer number indicating the enzyme targeted by the mutation
#'
#'
#' @return Numeric vector corresponding to mutant value of enzyme concentrations 
#' 
#' 
#' @seealso 
#' Use function \code{\link{compute.delta}} to compute the \bold{apparent} mutation effect.
#' 
#' Use function \code{\link{alpha_ij}} to compute matrix of redistribution coefficients.
#' 
#' See function \code{\link{mut.E.direct}} for a direct computation method of mutation.
#'
#'
#'
#' @examples
#' E <- c(30,30,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- apply(beta,1,sumbis)
#' correl <- "RegPos"
#' mu <- 1 #canonical size of mutation
#' 
#' alph <- alpha_ij(E,correl,beta)
#' delta <- compute.delta(mu,E,correl,B)
#' 
#' i <- 3 #enzyme directly targeted by mutation
#' mut.E.indirect(delta[i],E,alph,i)
#'
#' @export




#Compute mutant E from formula Ejm = Ejr + alpha_ij*(Eim-Eir)
mut.E.indirect <- function(delta_fun,E_res,alpha_fun,i_fun) {
  #verif param
  if (i_fun<1|i_fun>length(E_res)) {
    stop("You do not choose an enzyme in the sytem. i_fun need to be between 1 and n.")
  }
  ###
  E_mut_fun <- E_res + alpha_fun[i_fun,]*delta_fun
  return(E_mut_fun)
}


