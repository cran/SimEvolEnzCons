#' Co-regulation coefficient computation
#'
#' Computes the matrix of co-regulation coefficients \eqn{M_\beta} from a vector of global co-regulation coefficients \code{B}
#' 
#'
#' @details 
#' Enzymes are supposed to be all co-regulated, i.e. no \eqn{\beta} value is null.
#' 
#' 
#' @usage compute.beta.from.B(B_fun)
#'     
#' 
#' @param B_fun Numeric vector of global co-regulation coefficients
#' 
#' 
#'
#'
#' @return Numeric matrix \code{n*n} of the co-regulation coefficients, where \code{n} is the length of \code{B_fun}.
#' 
#' If \code{beta_fun} is \code{NULL}, \code{compute.beta.from.B} returns \code{NULL}.
#' 
#' 
#' @seealso 
#' Use function \code{\link{is.B.accurate}} to verify \code{B_fun} conformity.
#'
#'
#'
#' @examples
#' B <- 1/c(0.5,0.2,0.3)
#' correl <- "RegPos"
#' 
#' is.B.accurate(B,3,correl)
#' 
#' beta <- compute.beta.from.B(B)
#' 
#'
#' @export




compute.beta.from.B <- function(B_fun) {
  
  #add is.B.accurate() ?
  
  # if beta_fun is NULL
  if (length(B_fun)==0) {
    beta_fun <- NULL
  }
  
  # in other cases
  if (length(B_fun)>1) {
    #B*t(1/B)
    beta_fun <- B_fun %*% t(1/B_fun)
  }
  
  return(beta_fun)
}

