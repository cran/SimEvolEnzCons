#' Compare RNV mean in simulation with predicted RNV at equilibrium
#'
#' Compare RNV mean in simulation with predicted RNV at equilibrium to measure suitability of RNV mean as a proxy of RNV at equilibrium
#' With a graph and a linear model.
#'
#'
#' @details 
#' Function \code{RNV.mean.suitability} computes mean of RNV size by enzyme, using function \code{\link{RNV.mean.simul}}.
#' It computes also the RNV size around predicted equilibrium, using function \code{\link{RNV.size.at.equilibr}}.
#' Then \code{RNV.mean.suitability} plots the RNV means in relation to predicted RNVs for every enzyme and every simulation.
#'
#' \code{RNV.mean.suitability} computes also a linear model of RNV means in relation to RNV size at equilibrium, for selected simulations only (by setting \code{which.sim}).
#'   
#' Each simulation corresponds to one color. Colors for simulations are taken in palette \code{rainbow}.
#' Displayed number is the enzyme number.
#' Black line is the linear model. Dashed line is the symmetry line.
#' 
#'  Function \code{RNV.mean.suitability} is designed to compute mean RNV from simulations launched by \code{\link{simul.evol.enz.multiple}}.
#' Input \code{all_res_sim} is the output of \code{\link{simul.evol.enz.multiple}}.
#' 
#'   
#' 
#' @usage RNV.mean.suitability(all_res_sim,end.mean=TRUE,which.sim=NULL,
#' new.window=FALSE,posi.legend="topleft",...) 
#'
#'
#' @inheritParams RNV.mean.simul
#' @param posi.legend Character string. Where would you put the legend? See parameter \code{x} in function \code{legend}. If \code{NULL}, legend will not appear.
#'
#'
#' @return Invisible list of 6 elements:
#' \itemize{
#'    \item \code{$RNV_mean_simul}: numeric matrix of \code{n+2} columns and number of rows is between \code{nsim} and \code{2*nsim} (depending on RNV number).
#'    \code{n} first columns contain RNV mean for corresponding enzyme. Column \code{n+1} indicates simulation number and column \code{n+2} the RNV number (1 for near RNV and 2 for far RNV).
#'    \item \code{$list_eq_all}: list of \code{nsim} elements. Each element \code{s} is the output of function \code{\link{RNV.size.at.equilibr}} for simulation \code{s}
#'    \item \code{$RNV_at_eq}: numeric matrix of \code{n+1} columns and \code{nsim} rows.
#'    \code{n} first columns correspond to RNV size at predicted equilibrium for corresponding enzyme. Column \code{n+1} indicates simulation number.
#'    \item \code{$lm_compar_RNV}: object of \code{class "lm"}. Linear model of RNV mean for near RNV (RNV number 1) in relation to RNV at predicted equilibrium.
#'    \item \code{$A_mean}: numeric matrix of \code{n} columns (enzyme) and \code{nsim} rows (simulation). Each cell is the activity mean for each enzyme (in column) and each simulation (in row).
#'    \item \code{$Etot_mean}: numeric vector of length \code{nsim}. Each value is the total concentration mean for each simulation.
#' } 
#' 
#' @import graphics
#' @importFrom grDevices rainbow
#' 
#' @seealso 
#' Function \code{\link{RNV.for.simul}} is used to compute RNV. 
#' 
#' Function \code{\link{RNV.mean.simul}} is used to compute RNV mean for each simulation.
#' 
#' Function \code{\link{RNV.size.at.equilibr}} is used to compute RNV at equilibrium.
#' 
#' 
#' 
#'
#' @examples
#' 
#' # With saved simulation
#' data(data_sim_SC)
#' RNV.mean.suitability(data_sim_SC,new.window=TRUE,which.sim=c(1,4,6))
#' 
#' \donttest{
#' data(data_sim_RegNeg)
#' RNV.mean.suitability(data_sim_RegNeg,new.window=TRUE)
#' }
#'
#'
#' @export



RNV.mean.suitability <- function(all_res_sim,end.mean=TRUE,which.sim=NULL,new.window=FALSE,posi.legend="topleft",...) {
  
  ####### Take interessing parameters
  
  #input parameters
  n <- all_res_sim$param$n #number of enzymes
  nsim <- all_res_sim$param$nsim #number of simulations
  Keq_fun <- all_res_sim$param$Keq #equilibrium constants
  X_fun <- all_res_sim$param$X #numerator of flux
  beta_fun <- all_res_sim$param$beta #co-regulation coefficients
  B_fun <- all_res_sim$param$B #global co-regulation coefficients
  correl_fun <- all_res_sim$param$correl #constraints on system
  N_fun <- all_res_sim$param$N #population size
  pasobs <- all_res_sim$param$pasobs
  npt <- all_res_sim$param$npt #number of observations
  ngene <- npt*pasobs #number of generations
  pmutA <- all_res_sim$param$pmutA #proba mutation of kin
  same.E0 <- all_res_sim$param$same.E0
  is.random.E0 <- all_res_sim$param$is.random.E0
  same.kin0 <- all_res_sim$param$same.kin0
  is.random.kin0 <- all_res_sim$param$is.random.kin0
  
  
  #simulation results
  tabR <- all_res_sim$tabR
  tabP_e <- all_res_sim$tabP_e
  tabP_c <- all_res_sim$tabP_r
  
  #Which simulations if no chosen one
  if (length(which.sim)==0) {
    which.sim <- 1:nsim
  }
  
  ######## RNV
  #nearest RNV only
  which_RNV <- 1
  
  ## RNV in simul
  RNV_simul_mean <- RNV.mean.simul(all_res_sim,end.mean,which.sim,add.lm=FALSE,show.plot=FALSE,new.window=FALSE)
  #RNV_mean_size <- RNV_simul_mean$RNV_mean_size
  
  #nearest RNV only
  near.RNV <- which(RNV_simul_mean$RNV_mean_size[,(n+2)]==which_RNV)
  RNV_mean_size_near <- RNV_simul_mean$RNV_mean_size[near.RNV,]
  
  
  ## RNV at equilibrium
  list_eq_all <- list(NULL)
  RNV_at_eq <- NULL
  #for nsim simul and n enzymes
  A_mean <- matrix(NA,ncol=n,nrow=nsim)
  Etot_mean <- rep(NA,nsim) #numeric(length = nsim)
  
  for (s in which.sim) {
     #take current simulation result
    last_res <- tabR[tabR$sim==s,]
    
    # which generations are taken for mean computation ?
    if (end.mean==TRUE) {
      #mean of RNV size for half end of simul
      which_gen <- round(npt/2):npt
    } else { #mean for all resident
      which_gen <- 1:npt
    }
    ## compute mean of A and Etot to compare with mean RNV
    A_mean[s,] <- apply(last_res[which_gen,(2*n+4):(3*n+3)],2,mean,na.rm=TRUE)
    Etot_mean[s] <- mean(last_res[which_gen,(2*n+1)])
      
    #compute RNV at eq
    list_eq_all[[s]] <- RNV.size.at.equilibr(n,correl_fun,N_fun,A_fun=A_mean[s,],B_fun=B_fun,E0_fun=all_res_sim$list_init$E0[s,],Etot_eq=Etot_mean[s],show.plot=FALSE)
    
    #save RNV at eq
    #add simul number 's' (and RNV number)
    RNV_proeq <- cbind(list_eq_all[[s]]$RNV_size,s)#,1:nrow(list_eq_all[[s]]$RNV_size))
    RNV_at_eq <- rbind(RNV_at_eq,RNV_proeq[which_RNV,])
  }
  
  
  ########### Linear model
  lm_compar_RNV <- lm(as.vector(RNV_mean_size_near[,1:n])~as.vector(RNV_at_eq[,1:n]))
  
  
  ############# Graphics
  #color vector for simulations
  mycol_sim <- rainbow(nsim)
  
  #empty plot
  plot(0,0,type="n",xlab="RNV at equilibrium",ylab="RNV mean",xlim=c(0,max(RNV_at_eq[,1:n],RNV_mean_size_near[,1:n])),ylim=c(0,max(RNV_at_eq[,1:n],RNV_mean_size_near[,1:n])),...)
  #linear model
  abline(lm_compar_RNV)
  
  #if perfect correlation
  abline(0,1,lty=2)
  
  for (s in which.sim) {
    #take rows such as col n+1=s (simul s) and col n+2=1 (RNV 1), and take columns 1 to n
    #which.rows <- which(RNV_mean_size[,(n+1)]==s&RNV_mean_size[,(n+2)]==which_RNV)
    #take current simulation lines
    which.rows <- RNV_mean_size_near[,n+1]==s
    points(RNV_at_eq[which.rows,1:n],RNV_mean_size_near[which.rows,1:n],col=mycol_sim[s],type="o",pch=16)
    #add text corresponding to enzymes
    text(RNV_at_eq[which.rows,1:n],RNV_mean_size_near[which.rows,1:n],col=mycol_sim[s],labels=1:n,pos=3)
  }
  if (length(posi.legend)!=0) {
    legend(x=posi.legend,legend=c("symmetry","lm"),lty=c(2,1),lwd=1)
  }
  
  
  return(invisible(list(RNV_mean_simul=RNV_simul_mean$RNV_mean_size,RNV_at_eq=RNV_at_eq,lm_compar_RNV=lm_compar_RNV,list_all_eq=list_eq_all,A_mean=A_mean,Etot_mean=Etot_mean)))
}



