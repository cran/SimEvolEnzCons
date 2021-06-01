#' Equation of null coefficient selection
#'
#' \emph{(Utilitary function)}. Gives the result of the equation corresponding to discrete coefficient selection and indirect mutation of E
#' 
#' @details 
#' Gives the result of the equation corresponding to discrete coefficient selection and indirect mutation of E. 
#' Corresponding equation is \eqn{\sum ( \frac{1}{A_j E_j} (\frac{1}{1 + s} - \frac{E_j}{E_j + \alpha_{ij} \delta_i} ) )}
#' 
#' Null this equation corresponds to search the \eqn{\delta_i} corresponding to the given resident concentrations and selection coefficient.
#' 
#' This function is used to find \eqn{\delta} at bounds of the Range of Neutral Variation, where selection coefficient \emph{s} is equal to \eqn{+/-1/N}.
#'
#' 
#' 
#' 
#' @usage odd.discrete.sel.coef(delta_fun,i_fun,E_res,A_fun,alpha_fun,sel_coef_fun)
#'     
#' @inheritParams mut.E.indirect
#' @param A_fun Numeric vector of activities
#' @param sel_coef_fun Numeric value of selection coefficient
#'
#'
#' @return A numeric value
#' 
#' 
#' @seealso 
#' See function \code{\link{alpha_ij}} to compute the redistribution coefficients.
#' 
#' 
#'
#' 
#'
#' 





###another formulation to discrete coefficient selection from indirect mutation of E
odd.discrete.sel.coef <- function(delta_fun,i_fun,E_res,A_fun,alpha_fun,sel_coef_fun) {

  eqn_sel <- sum(1/(A_fun*E_res)*(1/(1+sel_coef_fun)-E_res/(E_res+alpha_fun[i_fun,]*delta_fun)))
  
  return(eqn_sel)
}
