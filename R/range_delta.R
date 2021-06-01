#' Bounds of \emph{delta_i}
#'
#' Computes the bounds of the actual mutation effect \eqn{\delta_i} such as all mutant concentrations are between 0 and total concentration, for a mutation targeting enzyme \eqn{i} 
#' 
#'
#' @details 
#' This function \code{range.delta} computes the bounds of \eqn{\delta_i} such as all mutant concentrations are between 0 and total concentration \eqn{Etot}, for a mutation targeting enzyme \code{i_fun}.
#' Mutant concentrations are equal to resident concentrations plus \eqn{\alpha_ij * \delta_i} \eqn{(see function \code{\link{mut.E.indirect}})}.
#' For any enzyme \eqn{j}, mutant value is \eqn{E_j^r + \alpha_ij * \delta_i}.
#' 
#' The inferior (resp. superior) bound of \eqn{\delta_i} corresponds to minimal (resp. maximal) value of \eqn{\delta_i}
#' such as all mutant concentrations are superior or equal to 0 \bold{and} inferior or equal to \eqn{Etot},
#' with at least one mutant concentration equal to 0 or \eqn{Etot}.
#' 
#' \code{tol_fun} is the accuracy (or allowed tolerance) for \eqn{\delta} bounds. It allows to avoid asymptote problem when computing the RNV.
#' 
#' 
#' 
#' 
#' @usage range_delta(E_res,alpha_fun,i_fun,tol_fun=0.0001)
#'     
#'
#' @inheritParams mut.E.indirect
#' @param tol_fun Numeric and positive value. Accuracy for delta bounds. Default is \code{0.0001}
#'
#'
#' @return Numeric vector of the inferior and the superior bounds of actual mutation effect \eqn{\delta_i}
#' 
#' 
#' @seealso 
#' See function \code{\link{alpha_ij}} to compute matrix of redistribution coefficients \code{alpha_fun}.
#' 
#'
#' @examples
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' Er <- c(30,30,30)
#' correl <- "CRPos"
#' alpha <- alpha_ij(Er,correl,beta)
#' 
#' #mutant enzyme
#' i <- 1
#'
#' range_delta(Er,alpha,i)
#'
#' @export





#Bounds of delta, such as Em is between 0 and Etot
range_delta <- function(E_res,alpha_fun,i_fun,tol_fun=0.0001) {

  ### Settings
  Etot_fun <- sum(E_res)
  n_fun <- length(E_res)
  
  ### Tests
  if (i_fun<=0 | i_fun>n_fun) {
    stop("Mutant enzyme number 'i_fun' is between 1 and n.")
  }
  if (tol_fun<=0) {
    stop("tol_fun is a positive numeric.")
  }
  
  #repetitive method
  # delta_lim_inf <- NULL
  # delta_lim_sup <- NULL
  # for (j_fun in 1:n_fun) {
  #   alpha_ij_fun <- alpha_fun[i_fun,j_fun]
  #   if (alpha_ij_fun>0) {
  #     delta_lim_inf <- c(delta_lim_inf,(-E_res[j_fun]/alpha_ij_fun))
  #     delta_lim_sup <- c(delta_lim_sup,(Etot_fun-E_res[j_fun])/alpha_ij_fun)
  #   }
  #   if (alpha_ij_fun<0) {
  #     delta_lim_inf <- c(delta_lim_inf,(Etot_fun-E_res[j_fun])/alpha_ij_fun)
  #     delta_lim_sup <- c(delta_lim_sup,(-E_res[j_fun]/alpha_ij_fun))
  #   }
  # }
  # delta_lim_fun <- c(max(delta_lim_inf)+tol_fun,min(delta_lim_sup)-tol_fun)
  
  
  
  
  #more elegant method
  
  #matrix of 2 columns (col 1 = inferior bounds ; col 2 = superior bounds) and 1 to n_fun rows depending on alpha_fun values
  all_delta_lim <- NULL
  
  #for each enzymes
  for (j_fun in 1:n_fun) {
    
    #value of redistribution coefficients alpha_ij between enzymes i_fun (mutant) and j_fun (other)
    alpha_ij_fun <- alpha_fun[i_fun,j_fun]
    
    #to keep bounds of delta for enzyme j_fun
    delta_lim_j <- NULL
    
    #if alpha_ij = 0, delta_i has no effect on E_res_j and E_res_j is unchanged, because E_mut_j = E_res_j + alpha_ij * delta_i
    if (alpha_ij_fun!=0) {
      # 0 < E_mut_j < - Etot and E_mut_j = E_res_j + alpha_ij * delta_i
      delta_lim_j <- c((-E_res[j_fun]/alpha_ij_fun),(Etot_fun-E_res[j_fun])/alpha_ij_fun)
    }
    
    #if alpha_ij = 0, bounds are reversed because alpha_ij is denominator
    if (alpha_ij_fun<0) {
      delta_lim_j <- rev(delta_lim_j)
    }
    
    #put bounds of delta_i for E_res_j with precedent results
    all_delta_lim <- rbind(all_delta_lim,delta_lim_j)
  }
  
  
  #more restreint zone for delta_i 
  #max(inferior bounds) and min(superior bounds) with a tolerance value
  delta_lim_fun <- c(max(all_delta_lim[,1])+tol_fun,min(all_delta_lim[,2])-tol_fun)
  
  
  #return a 2-length vector whith bounds of more restrained variation of delta_i
  return(delta_lim_fun)
}


