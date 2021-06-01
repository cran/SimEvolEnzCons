#' Time step simulation (next resident values)
#'
#' \code{simul.next.resident} simulates a time step, i.e. computes effects of a new mutation and gives next resident values whatever this mutation is fixed or not
#' 
#'
#' @details 
#' This function gives the genotype (enzyme concentrations and activities) and phenotype (flux) values of the next resident in an haploid population after a time step.
#' Here, a time step corresponds to the apparition and fixation (or disappearance) of a mutation targeting one enzyme. 
#' Therefore a time step is the interval between appearances of two successive mutations.
#' 
#' This function is used for simulation of enzyme concentration evolution.
#' 
#' \bold{Algorithm}
#' 
#' \enumerate{
#'    \item target enzyme and mutation sign are chosen randomly with a uniform law
#'    \item mutation targets randomly concentration or kinetic parameter depending on \code{pmutA} value
#'    \item mutation size is chosen randomly between 0 and \code{max_mut_size_E} for concentrations (resp. \code{max_mut_size_A} for kinetic parameters), then multiplied by its sign
#'    \item mutation effects on all enzymes are computed with function \code{\link{mut.E.direct}} (resp. \code{\link{mut.kin}}). The input mutation method is only available for \code{\link{mut.E.old}}, but multiplicative mutation method \code{typ_E=2} is not accurate in case of regulation
#'    \item if concentration or kinetic parameter become negative, which is biologically impossible, they are set to 0
#'    \item activities, then flux are computed
#'    \item selection coefficient is also computed, considering flux as fitness
#'    \item fixation probability of this mutation is computed
#'    \item fixation of this mutation is random, depending on this fixation probability: if mutation is fixed, mutant become resident for next step, else resident is unchanged for next step
#'    \item returns value of the resident for next step
#' }
#' Algorithm is also detailed in Coton et al. (2021)
#' 
#' 
#' 
#' @importFrom stats runif
#' 
#' 
#' @usage simul.next.resident(E_res_fun, kin_fun, Keq_fun, N_fun,
#' correl_fun, beta_fun=NULL, X_fun=1, max_mut_size_E=1, max_mut_size_A=1,
#' pmutA=0, typ_E=1, typ_A=1, use.old.mut=FALSE)
#'
#' @inheritParams mut.E.direct
#' @inheritParams activities
#' @inheritParams mut.kin
#' @param N_fun Numeric. Population size
#' @param X_fun Numeric. Numerator of function \code{\link{flux}}. Default is 1
#' @param max_mut_size_E Numeric. Maximum absolute size of mutation for enzyme concentrations. Default is 1
#' @param max_mut_size_A Numeric. Maximum absolute size of mutation for kinetic parameters. Default is 1
#' @param pmutA Numeric. Mutation probability of kinetic parameters.
#' Higher \code{pmutA}, higher the mutation probability of kinetic parameters compared to enzyme concentrations.
#' Default is 0, i.e. no mutation of kinetic parameters
#' @inheritParams mut.E.old
#' @param use.old.mut Logical. If \code{FALSE} (default), use \code{\link{mut.E.direct}} mutation method, else use \code{\link{mut.E.old}} mutation method if \code{TRUE}
#'
#'
#'
#' @return Returns a list of 7 elements:
#' \itemize{
#'    \item \code{$E_next}: numeric vector (length \code{n}) of the enzyme concentrations of next resident
#'    \item \code{$kin_next}: numeric vector (length \code{n}) of the kinetic parameters of next resident
#'    \item \code{$Etot_next}: numeric, the total concentration of next resident
#'    \item \code{$kintot_next}: numeric, the total kinetic of next resident
#'    \item \code{$J_next}: numeric, flux of next resident
#'    \item \code{$size_mut}: numeric, mutation size, even if it is not fixed
#'    \item \code{$target_mut}: numeric, number of enzyme targeted by the mutation
#' }
#' 
#' Note that \code{n} is the number of enzymes, which is the length of \code{E_res_fun}.
#' 
#' 
#' 
#' @seealso 
#' See function \code{\link{mut.E.direct}} and \code{\link{mut.kin}} to see how enzymes are mutated. 
#'
#' Fitness is computed with function \code{\link{flux}}.
#'
#' @references Coton et al. (2021)
#'
#' @examples
#' E <- c(30,30,30)
#' kin <- c(53/0.29,50/0.78,29)
#' Keq <- c(1.1e+8,4.9e+3,1.1e+3)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' correl <- "RegPos"
#' N <- 1000
#' 
#' simul.next.resident(E, kin, Keq, N, correl, beta, pmutA=0.1)
#'
#' @export



#Values of new resident after mutation (apparition and fixation or disappearance of the mutation)
simul.next.resident <- function(E_res_fun, kin_fun, Keq_fun, N_fun, correl_fun, beta_fun=NULL, 
                                X_fun=1, max_mut_size_E=1, max_mut_size_A=1, pmutA=0,
                                typ_E=1, typ_A=1, use.old.mut=FALSE) {
  
  ##########
  n_fun <- length(E_res_fun)
  
  ######## Mutation apparition
  
  # Generates mutation sign by a random uniform
  has <- runif(1,min=0,max=1)
  signe <- 1*(has<=0.5)-1*(has>0.5)
  
  # Choice of target enzyme
  #uniformly random between all enzymes
  i_fun <- ceiling(runif(1,min=0,max=n_fun))
  
  # Mutation size depend of which enzyme parameters is affected
  
  
  # Mutation of E or kin by a probability pmutA
  if (runif(1,min=0,max=1) > pmutA){
    # Mutation size of E
    #chosen in a random uniform
    nu_abs <- runif(1,min=0,max=max_mut_size_E)
    nu_fun <- signe*nu_abs
    
    # Mutation of E
    if (use.old.mut==FALSE) {
      E_m <- mut.E.direct(E_res_fun,i_fun,nu_fun,correl_fun,beta_fun)
    } else {
      E_m <- mut.E.old(E_res_fun,i_fun,nu_fun,correl_fun,beta_fun,typ_E)
    }
    kin_m <- kin_fun
    
  } else {
    # Mutation size of A
    #chosen in a random uniform
    nu_abs <- runif(1,min=0,max=max_mut_size_A)
    if (typ_A == 3) {
      #because mutation type 3 include only positive value
      nu_fun <- nu_abs
    } else {
      nu_fun <- signe*nu_abs
    }
    
    #Mutation of kin
    E_m <- E_res_fun
    kin_m <- mut.kin(kin_fun,i_fun,nu_fun,typ_A)
  } 
  
  # Set negative enzyme parameters to 0
  E_m[which(E_m<0)] <- 0
  kin_m[which(kin_m<0)] <- 0
  
  
  
  ############  Set new resident
  # i.e. fixation or disappereance of the mutation
  
  # Activities computation
  A_fun <- activities(kin_fun,Keq_fun)
  A_m <- activities(kin_m,Keq_fun)  
  
  # Fitness
  w_res <- flux(E_res_fun,A_fun,X_fun)
  w_m <- flux(E_m,A_m,X_fun)
  
  # Selection coefficient
  s_fun <- (w_m-w_res)/w_res
  #s_fun <- coef_sel.discrete(E_res,A_res,E_m,A_m)
  
  
  # Fixation probability of the mutation
  #haploid
  if (s_fun==0){
    #if mutation is strictly neutral
    pfix <- 1/(2*N_fun)
  } else {
    pfix <- 2*s_fun/(1-exp(-2*N_fun*s_fun))
  }
  
  # Fixation or not of the mutation
  if(pfix < runif(1,min=0,max=1)){
    #Disapperance of mutation
    E_mut <- E_res_fun
    kin_mut <- kin_fun
    w_mut <- w_res
    
  } else {
    #Fixation of mutation
    E_mut <- E_m
    kin_mut <- kin_m
    w_mut <- w_m
  }
  

  return(list(E_next=E_mut,kin_next=kin_mut,Etot_next=sum(E_mut),kintot_next=sum(kin_mut),w_next=w_mut,size_mut=nu_fun,target_mut=i_fun)) 
}

