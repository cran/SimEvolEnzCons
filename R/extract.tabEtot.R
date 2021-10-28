#' Table of concentrations
#'
#' Creates a table of all absolute concentrations from simulation results, or any table of enzyme concentrations (enzyme by column)
#' 
#'
#' @details 
#' From the result table \code{tabR} of simulations, extract the values of absolute concentrations of each enzymes,
#' and computes the absolute concentrations \eqn{E^q} in each regulation and the total concentration at each time step.
#' 
#' The input \code{tabR} is the output \code{tabR} of function \code{\link{simul.evol.enz.multiple}}, or possibly the output \code{res_sim} of function \code{\link{simul.evol.enz.multiple}}.
#' 
#' The function also works for any table \code{tabR} that has \eqn{n} columns (enzyme concentrations in columns), and any number of rows. The last column \code{"sim"} is not an obligation, but if present, need to be the last column.
#' 
#' @usage extract.tabEtot(tabR,beta_fun)
#'     
#'
#' @inheritParams class_group
#' @param tabR Dataframe of simulation results, which is an output of \code{\link{simul.evol.enz.multiple}}
#' 
#'
#'
#' @return Dataframe of \eqn{n+p+2} columns, where \eqn{n} is the number of enzymes and \eqn{p} is the number of regulation groups.
#' The \eqn{n} first columns are the concentrations of enzymes 1 to \eqn{n} (named \code{E1, ..., En}), the following \eqn{p} are the sum of concentrations in regulations groups from group 1 to \eqn{p} (named \code{Ep1, ..., Epp}), the second last column is the total concentration (named \code{Etot}) and the last column is the simulation number (named \code{sim}).
#' The column are named: \code{E1, ..., En, Eq1, ... Eqp, Etot, sim}.
#' The row number of the output dataframe is the same as \code{tabR}.
#' The last column \code{"sim"} is present only if input \code{tabR} has a column named \code{"sim"}.
#' 
#'
#'  
#' 
#' @seealso 
#' Use function \code{\link{simul.evol.enz.multiple}} to run simulations.
#'
#' @examples
#' #one group
#' n <- 3
#' beta <- diag(1,n)
#' beta[1,2] <- 0.1 ; beta[2,1] <- 1/beta[1,2]
#' some.E <- matrix(30,ncol=n,nrow=3)
#' some.E[2,] <- c(20,50,10)
#' some.E[3,] <- runif(n,1,100)
#' tabEtot <- extract.tabEtot(some.E,beta)
#' 
#' \donttest{
#' #With saved simulation
#' data(data_sim_RegPos)
#' tabEtot <- extract.tabEtot(data_sim_RegPos$tabR,data_sim_RegPos$param$beta)
#' }
#' 
#'
#' @export




#conversion of tabR in tabEtot
extract.tabEtot <- function(tabR,beta_fun) {
  
  #Number of enzymes
  n_fun <- nrow(beta_fun)
  #beta without reg
  # if (correl_fun=="Comp"|correl_fun=="SC") {
  #   beta_fun <- diag(1,n_fun)
  # }
  #list of groups
  L_Phi_fun <- class_group(beta_fun)
  p_fun <- length(L_Phi_fun)
  
  #create names of columns for tabEtot
  noms_E <- rep(0,(n_fun+p_fun+1))
  #Ei
  for (j in 1:n_fun) {
    noms_E[j] <- c(paste("E",j,sep=""))
  }
  #Eq
  for (q in 1:p_fun){
    noms_E[n_fun+q] <- c(paste("Ep",q,sep=""))
  }
  #Etot
  noms_E[n_fun+p_fun+1] <- "Etot"
  #n° simul
  if (any(names(tabR)=="sim")) {
    #if input tabR has a column "sim"
    noms_E[n_fun+p_fun+2] <- "sim"
  }
  
  
  
  #Recup data
  tabE <- tabR[,1:n_fun] #enz
  tabEtot <- tabE
  #for each group
  for (q in 1:p_fun){
    #if there is more than one enzyme in the group
    if (length(unlist(L_Phi_fun[q]))>1){
      #new column with the sum of Ei of the group
      tabEtot<-cbind(tabEtot,rowSums(tabE[,unlist(L_Phi_fun[q])])) #col n+q = E^q
    }
    #if there is one enzyme in the group
    if (length(unlist(L_Phi_fun[q]))==1){
      #new column with Ei of the enzyme
      tabEtot<-cbind(tabEtot,tabE[,unlist(L_Phi_fun[q])]) #col n+q = E^q
    }
  }
  #add Etot as new column
  tabEtot<-cbind(tabEtot,rowSums(tabE[,1:n_fun])) #Etot at end : col n+p+1
  #add n° simul as new column
  if (any(names(tabR)=="sim")) {
    #if input tabR has a column "sim"
    tabEtot<-cbind(tabEtot,tabR$sim) # and col n+p+2 =n° simul
  }
  
  #names all columns
  colnames(tabEtot)<- noms_E
  
  return(tabEtot)
}