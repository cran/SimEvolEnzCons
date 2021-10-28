#' Classify enzymes in regulation groups
#'
#' Classify the enzymes in their corresponding regulation group from the co-regulation matrix
#'  
#'
#' @details 
#' Enzymes are classified in regulation groups depending on the co-regulation matrix.
#' If the co-regulation coefficient is equal to zero between enzymes, the enzymes are in different regulation group, else there are in the same group.
#' 
#' Enzyme 1 is in the regulation group 1.
#' Other enzymes are classified in ascending order.
#' 
#' 
#' @usage class_group(beta_fun)
#'     
#' @param beta_fun Matrix of co-regulation coefficients
#' 
#'
#' @return Return a list classically named \code{L_Phi} of length \code{p} (the number of regulation groups). Each element \code{q} of \code{L_Phi} (the regulation group \eqn{\Phi_q}) contains the numbers of co-regulated enzymes.  
#' 
#' @seealso 
#' Function \code{\link{search_group}} to find the group of an enzyme.
#' 
#'
#' @examples
#' 
#' #Only one group
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' L_Phi <- class_group(beta)
#' 
#' #total number of regulation groups
#' p <- length(L_Phi) #gives 1
#' 
#' 
#' #Two groups
#' n <- 3
#' beta <- diag(1,n) 
#' beta[1,2] <- -0.32 
#' beta[2,1] <- 1/beta[1,2]
#' 
#' L_Phi <- class_group(beta)
#' p <- length(L_Phi) #gives 2
#' 
#' 
#' #with data
#' data(data_sim_RegNeg_1grpNeg1grpPos)
#' class_group(data_sim_RegNeg_1grpNeg1grpPos$param$beta)
#' 
#' @export


class_group <- function(beta_fun) { 
  #fonction permettant de classer les enzymes dans les groupes de régulation correspondant
  
  # total number of enzymes
  n_fun <- ncol(beta_fun)
  
  # Verif parameters
  #is.beta.accurate(beta_fun) but requires correl_fun
  if (nrow(beta_fun)!=ncol(beta_fun)) {
    stop("beta_fun need to have same number of rows and columns.")
  }
  
  
  # Initialization
  
  #L_Phi_fun = list of regulation groups, where each element q contains the numbers of enzymes in group Phi_q
  #first enzyme in group 1
  L_Phi_fun <- list(c(1))
  #number of group
  q=1
  
  #for each enzyme i to be classified, taken in rows of beta_fun (except enzyme 1)
  for (i in 2:n_fun) {
    
    #enzyme i is not yet classified, i.e. in a group
    flag=0
    
    #take each enzyme j (in columns of beta_fun) already classified (1 to i-1) to compare with enzyme i
    for (j in 1:(i-1)) {
      
      #if enzymes i and j are co-regulated and i is not already in a group
      if (beta_fun[i,j]!=0&flag==0) {
        
        #search the group k (among the previous group) that contains j
        k <- search_group(j,L_Phi_fun)
        #add i in that group k
        L_Phi_fun[k]<-list(c(unlist(L_Phi_fun[k]),i))
        #thus enzyme i is classified
        flag=1
      }
    }
    
    #after the loop, if enzyme i is not classified, i.e. is not regulated with previous enzymes
    if (flag==0) {
      #create a new group
      q <- q+1
      #put i in this new group
      L_Phi_fun[q] <- c(i)
      #now enzyme i is classified
      flag=1
    }
  
    #end of loop i  
  }
  
  
  return(L_Phi_fun)
}
