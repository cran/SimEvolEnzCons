#' Verification of constraint name
#'
#' Verifies if the constraint is authorized for other functions or not
#'  
#'
#' @details 
#' Verifies if constraint is included in the six allowed cases.
#' 
#' Possible constraints are listed below:
#'  \itemize{
#'    \item \code{"SC"}: independence between all enzymes 
#'    \item \code{"Comp"}: competition for resources 
#'    \item \code{"RegPos"}: positive regulation
#'    \item \code{"RegNeg"}: negative regulation 
#'    \item \code{"CRPos"}: competition plus positive regulation 
#'    \item \code{"CRNeg"}: competition plus negative regulation
#' }
#' 
#' 
#' @usage is.correl.authorized(correl_fun)
#'     
#' @param correl_fun Character string indicating the abbreviation of the constraint applied on the system
#' 
#' 
#'
#'
#' @return Returns TRUE if \code{correl_fun} is an authorized case, else stop
#' 
#' 
#' @seealso Need help to correctly write \code{correl_fun}? Use function \code{\link{name.correl}}.
#'
#'
#' @examples
#' is.correl.authorized("SC")
#' #returns TRUE
#' 
#' 
#'
#' @export


is.correl.authorized <- function(correl_fun) {
  #liste of authorized cases
  autor_correl <- c("SC","Comp","RegPos","RegNeg","CRPos","CRNeg")
  if (all(correl_fun!=autor_correl)) {
    is.auto_fun <- FALSE
    stop(paste("correl_fun '",correl_fun,"' is not an authorized constraint."))
  } else {
    is.auto_fun <- TRUE
  }
  return(is.auto_fun)
}

  
  