#' Old method for mutation of enzyme concentrations
#'
#' Computes the mutant value of enzyme concentrations with the old method (see Lion \emph{et al.} 2004).
#' 
#'
#' @details 
#' This mutation method is named \emph{old}, because it was used the first one used for evolution simulation.
#' Some improvements has been made since.
#' 
#' The main difference with mutation function \code{\link{mut.E.direct}} is the method to compute mutant values under competition constraint.
#'  In this function \code{mut.E.old}, canonical effect \eqn{\nu} is redistributed only on enzymes different from the target one,
#'  whereas in \code{mut.E.direct}, mutation canonical effect is redistributed between \bold{all} enzymes, including the target one.
#' 
#' Moreover, computation method in \code{mut.E.direct} is simplified.
#' 
#' Last, this \emph{old} method includes two model of mutations: additive and multiplicative.
#' 
#' \itemize{
#'    \item Additive method (\code{typ_E=1}): mutant concentration is the sum of the resident one plus size of mutation \eqn{\nu}
#'    \item Multiplicative method (\code{typ_E=2}): mutant concentration is the product of the resident one and \eqn{1+\nu}
#'    }
#'    Default is method 1.
#'    
#' In \code{mut.E.direct}, only method 1 is kept.
#' 
#' 
#' 
#' @usage mut.E.old(E_fun,i_fun,nu_fun,correl_fun,beta_fun=NULL,typ_E=1)
#'
#' @param E_fun Numeric vector of enzyme concentrations (resident)
#' @param i_fun Numeric value corresponding to number of the enzyme targeted by the mutation
#' @param nu_fun Numeric value of \bold{canonical} mutation effect
#' @inheritParams alpha_ij
#' @inheritParams is.correl.authorized
#' @param typ_E Numeric for mutation method. Authorized values: 1 or 2. Default is 1.
#'
#'
#' @return Numeric vector corresponding to mutant value of enzyme concentrations 
#' 
#' 
#' @seealso 
#' See function \code{\link{mut.E.direct}} for a direct computation method of mutation. 
#' 
#' @references 
#' Lion, S., F. Gabriel, B. Bost, J. Fiévet, C. Dillmann, and D. De Vienne, 2004. An extension to the metabolic control theory taking into account correlations between enzyme concentrations. European Journal of Biochemistry 271:4375–4391.
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
#' mut.E.old(E,i,mu,correl,beta)
#'
#' @export





# Mutation of concentrations for original simulations
mut.E.old <- function(E_fun,i_fun,nu_fun,correl_fun,beta_fun=NULL,typ_E=1){
  
  ########
  n_fun <- length(E_fun)
  Em <- E_fun
  
  #### Verif param
  is.correl.authorized(correl_fun)
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  if (i_fun<1|i_fun>n_fun) {
    stop("You do not choose an enzyme in the sytem. i_fun need to be between 1 and n.")
  }
  
  
  ####Mutation of target enzyme
  Em[i_fun] <- E_fun[i_fun] + nu_fun
  # additive
  if (typ_E == 1) {
    Em[i_fun]<-(E_fun[i_fun]+nu_fun)
    }
  # multiplicative
  if (typ_E == 2) {
    Em[i_fun]<-(E_fun[i_fun]*(1+nu_fun))
    }
  
  
  
  #Competition
  #relative proportion relative of enzymes other than target one are constant
  if (correl_fun == "Comp") {
    #for each enzymes
    for (r in 1:n_fun){
      #different from i, the target one
      if(r!=i_fun){
        #variation of enzyme concentrations in same proportions due to competition
        # Em_r = (E0_r/(Etot-E0_i))*(Etot-Em_i)
        # sum(Em) = sum(E)
        Em[r]<-(E_fun[r]/(sum(E_fun)-E_fun[i_fun]))*(sum(E_fun)-Em[i_fun])
        
        #if negative concentration
        if (Em[r]<0){
          Em[r]=0
        }
      }
    }
  }

  
  #Regulation
  if (correl_fun == "RegPos"|correl_fun =="RegNeg") {
    #for each enzymes
    for (r in 1:n_fun){
      #different from i, the target one
      if(r!=i_fun){
        #beta(ir) = (Em_r - E0_r)/(Em_i - E0_i)
        Em[r]<-E_fun[r] + beta_fun[i_fun,r]*(Em[i_fun]-E_fun[i_fun])
        
        #if negative concentration
        if (Em[r]<0){
          Em[r]=0
        }
      }
    }
  }
  
  
  #Competition and regulation
  if (correl_fun == "CRPos"|correl_fun =="CRNeg") {
    Em1 <- Em
    
    #Step 1: Regulation
    #for each enzyme
    for (r in 1:n_fun){
      #different from i, the target one
      if(r!=i_fun){
        #idem than regulation only
        Em1[r]<-E_fun[r] + beta_fun[i_fun,r]*(Em1[i_fun]-E_fun[i_fun])
      }
    }
    
    Em2<-Em1
    
    #Step 2: competition
    #for each enzyme
    for (r in 1:n_fun){
      #constant proportions between all enzymes
      Em2[r]<-Em1[r]/sum(Em1)*sum(E_fun)
      
      #if negatve concentration
      if (Em2[r]<0){
        Em2[r]=0
      }
    }
    #Note that: sum(Em2)=sum(Em)=E_tot !=sum(Em1)
    Em <- Em2
  }
  
  #Em[which(Em<0)] <- 0
  
  return(Em)
}



