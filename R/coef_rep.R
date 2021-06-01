#' Response coefficient computation
#'
#' Computes the response coefficients for each enzyme
#' 
#'
#' @details 
#' The response coefficients are influenced by the applied constraint.
#' 
#' Response coefficients are an extension of the control coefficients defined in the Metabolic Control Theory (Kacser and Burns, 1973).
#' For further details on how the response coefficient is calculated, see Lion \emph{et al.} (2004).
#' 
#' \code{A_fun} and \code{E_fun} need to have the same length.
#' The restrictions on \code{beta_fun} are explained in \code{\link{alpha_ij}}.
#' 
#' 
#' @usage coef_rep(E_fun,A_fun,correl_fun,beta_fun=NULL)
#'     
#'
#' @param A_fun Numeric vector of activities
#' @inheritParams alpha_ij
#'
#'
#' @return Numeric vector of the response coefficients of each enzyme
#' 
#' @seealso 
#' Use function \code{\link{activities}} to compute enzyme activities.
#' 
#'
#' @references 
#' Lion, S., F. Gabriel, B. Bost, J. Fiévet, C. Dillmann, and D. De Vienne, 2004.
#' An extension to the metabolic control theory taking into account correlations between enzyme concentrations.
#' European Journal of Biochemistry 271:4375–4391.
#'
#' @examples
#' A <- c(1,10,30)
#' E <- c(30,30,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' correl <- "SC"
#' 
#' response <- coef_rep(E,A,correl,beta)
#' 
#' 
#' # result : response = c(0.88235294, 0.08823529, 0.02941176)
#'
#' @export


# Computation of response coefficients depending on constraints
coef_rep <- function(E_fun,A_fun,correl_fun,beta_fun=NULL) {
  # Number of enzymes
  n_fun <- length(E_fun)
  
  ### Verif parameters
  if (length(A_fun)!=length(E_fun)) {
    stop("You forgot some enzymes. A_fun and E_fun need to have the same length.")
  }
  
  # Computation of the alpha coefficients matrix
  Alpha <- alpha_ij(E_fun,correl_fun,beta_fun)
  
  R_fun <- rep(0,n_fun)
  #for each enzyme
  for (i_fun in 1:n_fun) {
    R_fun[i_fun] <- E_fun[i_fun]*sum(Alpha[i_fun,]/(A_fun*(E_fun^2)))/sum(1/(A_fun*E_fun))
  }
  return(R_fun)
}

