#' Selection coefficient computation
#'
#' Computes the selection coefficient using the continuous expression \eqn{s_i = R_i*\delta_i/E_i}
#' 
#'
#' @details 
#' Computes the selection coefficient using a continuous expression \eqn{s_i = R_i*\delta_i/E_i}
#'     
#' Only mutations of concentrations are been considered.
#' 
#' \code{i_fun} is the number of the enzyme targeted by the mutation. It is an integer number between 0 and \code{n}, which is the total number of enzyme in the pathway.
#' If \code{i_fun} is between 1 and \code{n}, \code{delta_fun} needs to be a single value and function \code{coef_sel.continue} computes the selection coefficient of a mutation of actual effect \code{delta_fun} targeting \code{i_fun}.
#' If \code{i_fun} is set to 0, \code{delta_fun} needs to be a vector of same length as \code{E_res}. Each value of \code{delta_fun} is the actual effect of the mutation, and the position of this value in the vector is the target enzyme number.
#' Thus, to see the effect of a mutation of given actual effect on every enzyme, set \code{i_fun} to 0 and \code{delta_fun} has to be a vector of same length as \code{E_res}.
#' 
#' 
#' @usage coef_sel.continue(i_fun,E_res,A_fun,delta_fun,correl_fun,beta_fun=NULL)
#'     
#'
#' @inheritParams coef_rep
#' @inheritParams mut.E.indirect 
#' @param i_fun Integer number indicating the enzyme targeted by the mutation. \emph{See details}
#' @param delta_fun Numeric. Actual effect of a mutation targeting enzyme \code{i_fun}, i.e. \eqn{\delta_i}. \emph{See details}
#'
#'
#' @return Numeric value of the selection coefficient for the target enzyme.
#' 
#' If \code{i_fun} is set to 0, returns the numeric vector of the selection coefficients for the different enzyme.
#' 
#' 
#' @seealso 
#' Use function \code{\link{activities}} to compute enzyme activities.
#' 
#' @references 
#' Coton et al. (2021)
#'
#'
#'
#' @examples
#' 
#' #### Set
#' A <- c(1,10,30)
#' E <- c(30,30,30)
#' correl <- "CRPos"
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- apply(beta,1,sumbis)
#' 
#' ### Mutation
#' mu <- 1
#' i <- 3
#' delta <- compute.delta(mu,E,correl,B)
#' 
#' #for enzyme i
#' coef_sel.continue(i,E,A,delta[i],correl,beta)
#' 
#' #for all enzyme
#' coef_sel.continue(0,E,A,delta,correl,beta)
#' 
#'
#' @export





# Continuous selection coefficient
coef_sel.continue <- function(i_fun,E_res,A_fun,delta_fun,correl_fun,beta_fun=NULL) {
  
  n_fun <- length(E_res)
  
  ### Verif parameters
  is.correl.authorized(correl_fun)
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  if (length(A_fun)!=n_fun) {
    stop("You forgot some enzymes. A_fun and E_res need to have the same length.")
  }
  if (length(delta_fun) != 1 & length(delta_fun) != n_fun) {
    stop(" 'delta_fun' is not a convenient data.")
  }
  if (i_fun==0 & length(delta_fun) != n_fun) {
    stop("Set the same number of enzyme in 'delta_fun' as in 'E_res' or select a target enzyme in 'i_fun'.")
  }
  if (i_fun<0 | i_fun>n_fun) {
    stop("Select a target enzyme.")
  }
  
  
  ### Compute response coefficient
  R_fun <- coef_rep(E_res,A_fun,correl_fun,beta_fun)
  #delta_fun <- compute.delta(mu_fun,E_res,correl_fun,B_fun)
  
  ### Selection coefficient
  if (i_fun ==0) {
    s_fun <- R_fun*delta_fun/E_res
  } else {
    s_fun <- R_fun[i_fun]*delta_fun/E_res[i_fun]
  }
  
  return(s_fun)
}

