#' Verification of B accuracy
#'
#' Verifies if the vector \code{B_fun} of global co-regulation coefficients is accurate for other functions
#'  
#'
#' @details 
#' Different tests are performed on parameter \code{B_fun} to verify its accuracy.
#' 
#' \itemize{
#'    \item Is there regulation? If yes, \code{B_fun} is necessary.
#'    \item Do \code{B_fun} have a correct length? Compare \code{length(B_fun)} and number of enzymes \code{n_fun}. If difference, stops.
#'    \item Is there negative regulation in \code{correl_fun}? If yes, does \code{B_fun} include a negative regulation?  If difference, stops.
#'    \item Sum of \code{1/B_fun} need to be equal to 1.
#'    }
#' 
#' 
#' @usage is.B.accurate(B_fun, n_fun, correl_fun)
#'     
#' @param B_fun Numeric vector of global co-regulation coefficients
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
#' To verify matrix of co-regulation coefficients, see function \code{\link{is.beta.accurate}}.
#' 
#' 
#'
#'
#' @examples
#' 
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- apply(beta,1,sumbis) 
#' is.B.accurate(B,3,"CRPos")
#' 
#'
#' @export


is.B.accurate <- function(B_fun, n_fun, correl_fun) {
  
  #value of B in these cases is not important but normally equal to 1 for all enzyme
  # this writing avoid multiple test
  # sum(1/B) is equal to the number of groups of regulation
  if (correl_fun=="SC"|correl_fun=="Comp") {
    B_fun <- rep(1,n_fun)*n_fun
  }
  #in other cases (presence of regulation), B is required
  if (length(B_fun)==0) {
    stop("B_fun is required.")
  }
  
  # verif if a parameter is missing
  if (length(B_fun)!=n_fun) {
    stop("An enzyme is missing in B_fun. ")
  } 
  
  # verif adequacy between correl and B values
  if (all(B_fun>=0)&(correl_fun=="RegNeg"|correl_fun=="CRNeg")) {
    stop("correl_fun and B_fun are not consistent. correl_fun includes a negative regulation and not B_fun. ")
  }
  if (any(B_fun<0)&(correl_fun=="RegPos"|correl_fun=="CRPos")) {
    stop("correl_fun and B_fun are not consistent. B_fun includes a negative regulation and not correl_fun. ")
  }
  
  # verif if sum(1/B) is equal to 1
  # if (class(all.equal(sumbis(1/B_fun),1))=="character") {
  #   stop("sum(1/B_fun) is not equal to 1. Not convenient data. ")
  # }
  
  # verif if sum(1/B) is an integer (in regulation group, =p)
  if (class(all.equal(sumbis(1/B_fun),round(sumbis(1/B_fun))))=="character") {
    stop("sum(1/B_fun) is not an integer. Not convenient data. ")
  }
  

  return(TRUE)
}

