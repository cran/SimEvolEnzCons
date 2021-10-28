#' Which types are the different regulation groups?
#'
#' Determines the type of all regulation groups from the matrix of co-regulation coefficients
#'  
#'
#' @details 
#' Regulation groups exist in three types :
#' \itemize{
#'  \item positive group, when all regulation coefficients are positive
#'  \item negative group, when it exists at least one negative regulation coefficients in the group
#'  \item singletons, when enzyme is lone in the group, and is therefore independent from all others enzymes
#' }
#' 
#' The position of the groups is determined from the list of regulation, computed with function \code{\link{class_group}}.
#' 
#' 
#' @usage group_types(beta_fun)
#'     
#' @inheritParams class_group
#' 
#'
#' @return Return a list of three elements that contains the numbers of the regulation groups:
#' \itemize{
#'  \item \code{$grp_pos}: numeric vector containing the position of positive groups
#'  \item \code{$grp_single}: numeric vector containing the position of singletons
#'  \item \code{$grp_neg}: numeric vector containing the position of negative groups
#' }
#' If there is no group of a type, the corresponding element returns \code{NULL} instead of a numeric vector.
#' 
#' @seealso 
#' Function \code{\link{class_group}} to compute the list of regulation groups
#' 
#'
#' @examples
#' 
#' #Only one group
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' L_Phi <- class_group(beta)
#' group_types(beta) #gives c(1), NULL and NULL
#' 
#' 
#' #Two groups
#' n <- 3
#' beta <- diag(1,n) 
#' beta[1,2] <- -0.32 
#' beta[2,1] <- 1/beta[1,2]
#' 
#' L_Phi <- class_group(beta)
#' group_types(beta) #gives NULL, c(2) and c(1)
#' 
#' #with data
#' data(data_sim_RegNeg_1grpNeg1grpPos)
#' group_types(data_sim_RegNeg_1grpNeg1grpPos$param$beta)
#' 
#' 
#' @export

group_types <- function(beta_fun) {
  
  ## Parameters
  #nb enzymes
  n_fun <- nrow(beta_fun)
  #list of groups
  L_Phi_fun <- class_group(beta_fun)
  #total number of groups
  p <- length(L_Phi_fun)
  #global regulation coefficients
  B_fun <- compute.B.from.beta(beta_fun)
  
  ## Initialization
  #which positive groups (only positive regulation coefficients)
  grp_pos <- NULL
  #which negative groups (at least one negative regulation coefficients)
  grp_neg <- NULL
  #which singletons (only one enzyme in group)
  grp_single <- NULL
  
  #for each group
  for (q in 1:p) {
    #take the list of enzymes in this group q
    phi_q <- unlist(L_Phi_fun[q])
    #total number of enzyme in the group
    m_q <- length(phi_q)
    
    #if enzyme is alone in this group
    if (m_q==1) {
      #save the number of this group
      grp_single <- c(grp_single,q)
    } else { #if there is more than one enzyme
      
      #if the group has negative regulation (at least one B is negative)
      if (any(B_fun[phi_q]<0)) { #(sum(B[phi_q]<0)!=0)
        #save this number as negative group
        grp_neg <- c(grp_neg,q)
      }
      
      #if the group has all positive regulation
      if (all(B_fun[phi_q]>0)) {
        #save this number as positive group
        grp_pos <- c(grp_pos,q)
      }
    }
  #end of loop on q
  }
    
  return(list(grp_pos=grp_pos, grp_single=grp_single, grp_neg=grp_neg))
}

