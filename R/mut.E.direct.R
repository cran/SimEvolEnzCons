#' Direct mutation of enzyme concentrations
#'
#' Computes the mutant value of enzyme concentrations by a direct method.
#' 
#'
#' @details 
#' This mutation method is named \emph{direct}, because we used canonical mutation effect to compute mutant values
#' 
#' This function is the one used in evolution simulation.
#' 
#' 
#' @usage mut.E.direct(E_res_fun,i_fun,nu_fun,correl_fun,beta_fun=NULL)
#'
#' @param E_res_fun Numeric vector of enzyme concentrations (resident)
#' @param i_fun Integer number indicating the enzyme targeted by the mutation
#' @param nu_fun Numeric value of \bold{canonical} mutation effect
#' @inheritParams alpha_ij
#' @inheritParams is.correl.authorized
#'
#'
#' @return Numeric vector corresponding to mutant value of enzyme concentrations 
#' 
#' 
#' @seealso 
#' See function \code{\link{mut.E.indirect}} for an indirect computation method of mutation. 
#'
#'
#'
#'
#' @examples
#' E <- c(30,30,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' correl <- "RegPos"
#' mu <- 1 #canonical size of mutation
#' i <- 3 #enzyme directly targeted by mutation
#' 
#' mut.E.direct(E,i,mu,correl,beta)
#'
#' @export




# Mutation of concentrations for simulations
mut.E.direct <- function(E_res_fun,i_fun,nu_fun,correl_fun,beta_fun=NULL) {
  #NB : mutations giving negative concentrations do not be eliminated in this function
  
  #######
  #Number of enzymes
  n_fun <- length(E_res_fun)
  #Total concentrations
  Etot_fun <- sum(E_res_fun)
  #If no change
  E_mut_fun <- E_res_fun
  
  ######## Verif parameters
  is.correl.authorized(correl_fun)
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  if (i_fun<1|i_fun>n_fun) {
    stop("You do not choose an enzyme in the sytem. i_fun need to be between 1 and n.")
  }
  
  
  #Independency
  if (correl_fun=="SC") {
    # only mutation of target enzyme
    E_mut_fun[i_fun] <- E_res_fun[i_fun] + nu_fun
  }
  
  
  #Competition
  #relative proportion relative of ALL enzymes (after mutation of i) are constant
  if (correl_fun == "Comp") {
    E_pre_fun <- E_res_fun
    #Step 1: canonical mutation of enzyme i
    E_pre_fun[i_fun] <- E_res_fun[i_fun] + nu_fun
    #Step 2: redistribution of mutation effect between all enzymes
    E_mut_fun <- Etot_fun*E_pre_fun/sum(E_pre_fun)
  }
  
  
  #Regulation
  if (correl_fun == "RegPos"|correl_fun =="RegNeg") {
    E_mut_fun <- E_res_fun + beta_fun[i_fun,]*nu_fun
  }
  
  
  #Competition and regulation
  if (correl_fun == "CRPos"|correl_fun =="CRNeg") {
    #Step 1 : mutation on target enzyme and regulation
    E_pre_fun <- E_res_fun + beta_fun[i_fun,]*nu_fun
    #Step 2 : competition = redistribution of modifications among all enzymes
    E_mut_fun <- Etot_fun*E_pre_fun/sum(E_pre_fun)
  }
  
  return(E_mut_fun)
}

