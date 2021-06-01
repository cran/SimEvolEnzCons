#' Real solutions of quadratic equation
#'
#' Computes real solutions of quadratic equation
#' 
#'
#' @details 
#' Quadratic equation is a second-degree polynomial equation of type \eqn{a x^2 + b x + c = 0}, 
#' where \code{a} is the quadratic coefficient, \code{b} the linear coefficient and \code{c} the free term.
#' 
#' 
#' 
#' @usage solv.2dg.polynom(a_fun,b_fun,c_fun)
#'     
#'
#' @param a_fun Numeric. The quadratic coefficient applied to \code{x^2}
#' @param b_fun Numeric. The linear coefficient applied to \code{x}
#' @param c_fun Numeric. The free term
#'
#'
#' @return Three possible vectors:
#' \itemize{
#'    \item Numeric vector of length 2 if there is two real roots
#'    \item Numeric value if there is a double-root
#'    \item \code{NULL} if there is no real solution
#' }
#'
#' @examples
#' solv.2dg.polynom(3,2,1)
#' #result : NULL
#' 
#' solv.2dg.polynom(1,2,1)
#' #result : -1
#' 
#' solv.2dg.polynom(1,0,-1)
#' #result : c(1,-1)
#'
#' @rdname solv.2dg.polynom
#' @export



#Search real solutions for a second degre polynomial
solv.2dg.polynom <- function(a_fun,b_fun,c_fun) {
  
  #computes discriminant
  Deter_fun <- b_fun^2-4*a_fun*c_fun
  
  #positive discriminant => 2 solutions
  if (Deter_fun>0) {
    x1 <- (-b_fun +sqrt(Deter_fun))/(2*a_fun)
    x2 <- (-b_fun -sqrt(Deter_fun))/(2*a_fun)
    sol_fun <- c(x1,x2)
  }
  
  #null discriminant => 1 double solution
  if (Deter_fun==0) {
    x0 <- (-b_fun)/(2*a_fun)
    sol_fun <- x0
  }
  
  #negative discriminant => no real solutions
  if (Deter_fun<0) {sol_fun <- NULL}
  
  
  return(sol_fun)
}

