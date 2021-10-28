#' Verification of beta matrix accuracy
#'
#' Verifies if the the matrix \code{beta_fun} of co-regulation coefficients is accurate for other functions
#'  
#'
#' @details 
#' Different tests are performed on matrix \code{beta_fun} to verify its accuracy.
#' 
#' \itemize{
#'    \item Is there regulation in \code{correl_fun}? If yes, \code{beta_fun} is necessary.
#'    \item Does \code{beta_fun} have a correct size? Compare \code{nrow(beta_fun)} and \code{ncol(beta_fun)} to number of enzymes \code{n_fun}. If difference, stops.
#'    \item Is there negative regulation in \code{correl_fun}? If yes, does \code{beta_fun} include a negative regulation?  If difference, stops.
#'    \item Each element of \code{beta_fun} diagonal needs to be equal to 1.
#'    }
#' 
#' 
#' @usage is.beta.accurate(beta_fun, n_fun, correl_fun)
#'     
#' @param beta_fun Numeric matrix of co-regulation coefficients
#' @param n_fun Number of enzymes in the system
#' @param correl_fun Character string indicating the constraint applied on the system
#' 
#'
#'
#' @return Return \code{TRUE} if all conditions are respected, else stop
#' 
#'
#'
#' @seealso 
#' To verify vector of global co-regulation coefficients, see function \code{\link{is.B.accurate}}. 
#' 
#'
#'
#' @examples
#' 
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' is.beta.accurate(beta,3,"CRPos")
#' 
#'
#' @export


is.beta.accurate <- function(beta_fun, n_fun, correl_fun) {
  
  #verif if available correlation
  is.correl.authorized(correl_fun)
  
  #value of beta in these cases is not important but normally equal to 0 for all enzyme, except in the diagonal which is 1
  # this writing avoid multiple test
  if (correl_fun=="SC"|correl_fun=="Comp") {
    beta_fun <- matrix(0,nrow=n_fun,ncol=n_fun)
    diag(beta_fun) <- 1
  }
  #in other cases (presence of regulation), beta is required
  if (length(beta_fun)==0) {
    stop("beta_fun is required.")
  }
  
  # verif if a parameter is missing
  if (nrow(beta_fun)!=ncol(beta_fun)) {
    stop("beta_fun need to have same number of rows and columns.")
  }
  if (nrow(beta_fun)!=n_fun) {
    stop("An enzyme is missing in beta_fun. ")
  } 
  
  # verif adequacy between correl and B values
  if (all(beta_fun>=0)&(correl_fun=="RegNeg"|correl_fun=="CRNeg")) {
    stop("correl_fun and beta_fun are not consistent. correl_fun includes a negative regulation and not beta_fun. ")
  }
  if (any(beta_fun<0)&(correl_fun=="RegPos"|correl_fun=="CRPos")) {
    stop("correl_fun and beta_fun are not consistent. beta_fun includes a negative regulation and not correl_fun. ")
  }
  
  # verif if elements of beta diagonal are equal to 1
  if (class(all.equal(diag(beta_fun),rep(1,n_fun)))=="character") {
    stop("Diagonal elements are not equal to 1. Not convenient data. ")
  }
  

  return(TRUE)
}

