#' Mean of RNV size in simulation
#'
#' Computes mean of RNV size from simulation results and gives a graph of this RNV mean in relation to the RNV-ranking-order factor.
#' Computes also a linear model of RNV mean in relation to RNV-ranking-order factor.
#'
#' @details 
#'
#' 
#' \code{RNV.mean.simul} works in three parts: \enumerate{
#'    \item computing mean of RNV size
#'    \item plotting RNV mean in relation to various variables
#'    \item computing RNV mean against RNV-ranking-order factor 
#' }
#'   
#' \bold{About RNV mean computing}
#' 
#' Function \code{RNV.mean.simul} is designed to compute mean of RNV size in simulations launched by \code{\link{simul.evol.enz.multiple}}.
#' Input \code{all_res_sim} is the output of \code{\link{simul.evol.enz.multiple}}.
#' 
#' RNV mean is computed by enzyme and by simulation, and a general mean for each enzyme (between selected simulations) is also computed.
#' 
#' RNV mean is made on all rows in simulation results (\code{end.mean=FALSE}) or only on last half rows of each simulation (\code{end.mean=TRUE}), i.e. when equilibrium is reached.
#' 
#' 
#' \bold{About graphics}
#' 
#' Function \code{RNV.mean.simul} gives three graphics: \enumerate{
#'    \item RNV mean in relation to activities
#'    \item RNV mean in relation to the RNV-ranking-order factor (see \code{\link{RNV.ranking.order.factor}})
#'    \item RNV mean in relation to a soft value of the RNV-ranking-order factor (activities A for \code{"SC"} and \code{"Comp"}; global co-regulations coefficients B for \code{"RegPos"} and \code{"RegNeg"}; the hard value of the RNV-ranking-order factor for \code{"CRPos"} and \code{"CRNeg"}). Squares correspond to simulation.
#' }
#' 
#' Each simulation corresponds to one color. Colors for simulations are taken in palette \code{rainbow}.
#' 
#' 
#' \bold{About linear model}
#' 
#' Function \code{RNV.mean.simul} computes also a linear model of RNV mean in relation to the RNV-ranking-order factor between \bold{all} simulations (and not only selected ones by \code{which.sim}).
#' If wanted, linear model can be put on graphics.
#' 
#' 
#' 
#' \bold{About logical parameters}
#' 
#' Last graphic is RNV mean in relation to an interest variable, which is activities \bold{A} (\code{"SC"} or \code{"Comp"} cases), global co-regulation coefficients \bold{B} (\code{"RegPos"} or \code{"RegNeg"} cases) or the RNV-ranking-order factor (\code{"CRPos"} or \code{"CRNeg"} cases, see above).
#' 
#' In this last graphic, \code{add.mean=TRUE} adds the mean of RNV size between selected simulations, with black squares and line.
#' 
#' Also in this last graphic, \code{add.pred.e=TRUE} adds mean (between selected simulations) of the predicted relative concentrations at equilibrium, with grey axis, grey dots and grey dashed line.
#' 
#' \code{add.lm=TRUE} adds the linear model (black line) in the second graph (RNV mean against the RNV-ranking-order-factor).
#' 
#' 
#' @usage RNV.mean.simul(all_res_sim,end.mean=TRUE,which.sim=NULL,add.lm=TRUE,
#' add.mean=TRUE,add.pred.e=FALSE,
#' show.plot=TRUE,new.window=FALSE,cex.lab=1,mar.lab=2.5,...) 
#'
#'
#' @inheritParams RNV.for.simul
#' @inheritParams graph.simul.by.time.RNV
#' @param add.lm Logical. Add line of linear model in graphics?
#' @param show.plot Logical. Are plots visible?
#' @param add.mean Logical. Add mean of RNV size between all selected simulation for each enzyme?
#' @param add.pred.e Logical. Add predicted relative concentrations mean between selected simulations? \emph{See concerned graphs in details }
#' @inheritParams RNV.graph.double.at.eq
#' 
#'
#' @return Invisible list of 5 elements:
#' \itemize{
#'    \item \code{$RNV_all_sim}: list of \code{nsim} elements (which is the number of simulation). Each element \code{i} contains the output of function \code{\link{RNV.for.simul}} for corresponding simulation \code{i}.
#'    If simulation is not selected by \code{which.sim}, corresponding element is \code{NULL}. 
#'    \item \code{$RNV_mean_size}: numeric matrix of \code{n+2} columns and number of rows is between \code{length(which.sim)} and \code{2*length(which.sim)} (depending if there is one or two RNVs).
#'    Each of the \code{n} first columns is the RNV mean size for corresponding enzyme. Column \code{n+1} indicates the simulation number and column \code{n+2} the RNV number (1 for near RNV and 2 for far RNV).
#'    \item \code{$rank_var_value}: numeric matrix of \code{n} columns and \code{length(which.sim)} rows. Each cell is the value of the RNV-ranking-order factor to which RNV mean size is compared, for each simulation (in row) (in column) and each enzyme.
#'    \item \code{$rank_var_name}: character string, indicating the name of the RNV-ranking-order factor.
#'    \item \code{$lm_RNV}: object of \code{class "lm"}. Linear model of RNV the mean (only for near RNV or RNV number 1) in relation to the RNV-ranking-order factor.
#' } 
#' 
#' @import graphics
#' @importFrom grDevices rainbow
#' 
#' @seealso 
#' RNV is computed with function \code{\link{RNV.for.simul}}. 
#' 
#' Use function \code{\link{graph.simul.by.time.RNV}} to have other representations of RNV.
#' 
#' See function \code{\link{RNV.ranking.order.factor}} for details about RNV-ranking-order factor.
#'
#' @examples
#' 
#' # With saved simulation
#' data(data_sim_SC)
#' RNV.mean.simul(data_sim_SC,new.window=TRUE,which.sim=c(1,5,10))
#' 
#' 
#' \donttest{
#' # case for 3 enzymes
#' n <- 3
#' E0 <- c(30,30,30)
#' kin <- c(1,10,30)
#' Keq <- c(1,1,1)
#' nsim <- 2 # 2 simulations
#' N <- 1000
#' correl <- "SC"
#' 
#' evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,npt=250)
#' 
#' RNV.mean.simul(evol_sim,new.window=TRUE)
#' }
#'
#'
#' @export





RNV.mean.simul <- function(all_res_sim,end.mean=TRUE,which.sim=NULL,add.lm=TRUE,
                           add.mean=TRUE,add.pred.e=FALSE,show.plot=TRUE,new.window=FALSE,cex.lab=1,mar.lab=2.5,...) {
  
  
  #here, variable 's' is always used to indicate simulation number and 'i' for mutant enzyme number
  
  ####### Take interessing parameters
  #save par() of user
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
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
  
  #initial and final values of E, kin, A and relative concentrations e
  #matrix such as: each simulation in row, initial values from column 1 to n and final values from column n+1 to 2n
  leg_E <- cbind(all_res_sim$list_init$E0,all_res_sim$list_final$E_f)
  leg_kin <- cbind(all_res_sim$list_init$kin0,all_res_sim$list_final$kin_f)
  leg_A <- cbind(all_res_sim$list_init$A0,all_res_sim$list_final$A_f)
  leg_e <- leg_E
  for (i in 1:nsim){
    #initial relative concentrations e0
    leg_e[i,1:n] <- leg_E[i,1:n]/sum(leg_E[i,1:n])
    #final relative concentrations e_f
    leg_e[i,((n+1):(2*n))] <- leg_E[i,((n+1):(2*n))]/sum(leg_E[i,((n+1):(2*n))])
  }
  
  
  ##########
  # Compute RNV
  #######"
  
  #which simulations are taken
  if (length(which.sim)==0) {which.sim <- 1:nsim}
  
  #save RNV for all simulations
  RNV_all_sim <- vector("list",length=nsim)
  #take mean RNV size for all simul
  RNV_mean_size <- NULL
  
  for (s in which.sim) {
    
    #Take result for simulation s
    res_sim <- tabR[tabR$sim==s,]
    #compute RNV elements
    RNV_all <- RNV.for.simul(res_sim,n,N_fun,correl_fun,beta_fun,end.mean)
    #save
    RNV_all_sim[[s]] <- RNV_all
    
    nb_RNV <- RNV_all$nb_RNV
    
    RNV_dis_delta_graph <- RNV_all$RNV_delta
    RNV_dis_enz_graph <- RNV_all$RNV_enz
    RNV_dis_size <- RNV_all$RNV_size
    RNV_dis_size_divEtot <- RNV_all$RNV_size_divEtot
    RNV_proxy <- RNV_all$RNV_proxy
    RNV_proxy_divEtot <- RNV_all$RNV_proxy_divEtot
    RNV_dis_flux <- RNV_all$RNV_flux
    limits_NZ <- RNV_all$limits_NZ
    
    #save mean RNV
    #add simul number 's' and RNV number
    RNV_promean <- cbind(RNV_proxy,s,1:nb_RNV)
    RNV_mean_size <- rbind(RNV_mean_size,RNV_promean)
    
  }
  
  #nearest RNV only
  near.RNV <- which(RNV_mean_size[,(n+2)]==1)
  
  #mean for all selected simulations by enzyme, only for nearest RNV
  RNV_mean_by_enz <- apply(RNV_mean_size[near.RNV,1:n],2,mean)
  
  
  
  
    ########### RNV-ranking-order variable
  
  #for each chosen simulation, ranking variable is different
  rank_var_sim <- matrix(0,nrow=length(which.sim),ncol=n)
  #name depends only on constraints (correl_fun) but defined also by B_fun
  rank_var_name <- RNV.ranking.order.factor(rep(1,n),correl_fun,rep(1,n),B_fun)$name
  #incr -> easier than which(which.sim==s) to follow how many simulations have been seen
  incr <- 0
  #each simul
  for (s in which.sim) {
    incr <- incr +1
    rank_var_sim[incr,] <- RNV.ranking.order.factor(leg_A[s,1:n],correl_fun,leg_E[s,1:n],B_fun)$value
  }

  # int_var <- leg_A[,1:n] #A initiaux
  
  
  #interest variable : matrix of nsim rows and n columns correpsonding to which RNv mean need to be compared
  int_var_sim <- matrix(0,nrow=nsim,ncol=n)
  if (correl_fun=="SC"|correl_fun=="Comp") {
    int_var_sim <- leg_A[,1:n]
    int_var_name <- "Pseudo-activity A"
  }
  if (correl_fun=="RegPos"|correl_fun=="RegNeg") {
    int_var_sim <- matrix(B_fun,byrow=TRUE,nrow=nsim,ncol=n)
    int_var_name <- "Global co-regulation coefficient B"
  }
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    B_sim <- matrix(1/B_fun,byrow=TRUE,nrow=nsim,ncol=n)
    int_var_sim <- abs(B_sim - leg_e[,1:n])
    int_var_name <- "Absolute value of (1/B - e0)"
  }
  if (add.mean==TRUE) {
    #mean of interest variable between all selected simulation
    int_var_mean <- apply(int_var_sim[which.sim,], 2, mean)
    #take mean between selected simulations of LAST value for predict e 
    #possibly, use of end.mean for mean of predict e rather than last value
    pred_e_mean <- apply(tabP_e[which.sim*npt,1:n], 2, mean)
    #combine all mean for ordered them
    RNV_and_var_mean <- cbind.data.frame(int=int_var_mean,RNV=RNV_mean_by_enz,e=pred_e_mean)
    #order depending on interest variable
    RNV_and_var_mean_ord <- RNV_and_var_mean[order(RNV_and_var_mean$int),]
  }
  
  
  # #interest variable : matrix of nsim rows and n columns correpsonding to which RNv mean need to be compared
  # if (correl_fun=="SC") {
  #   int_var <- leg_A[,1:n]^(1/3)
  #   int_var_name <- "Activities^(1/3)"
  # }
  # if (correl_fun=="Comp") {
  #   int_var <- leg_A[,1:n]^(-1/4)
  #   int_var_name <- "Activities^(-1/4)"
  # }
  # if (correl_fun=="RegPos"|correl_fun=="RegNeg") {
  #   int_var <- matrix(abs(1/B_fun),byrow=TRUE,nrow=nsim,ncol=n)
  #   int_var_name <- "Absolute of 1/B"
  # }
  # if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
  #   B_sim <- matrix(1/B_fun,byrow=TRUE,nrow=nsim,ncol=n)
  #   int_var <- abs(B_sim - leg_e[,1:n])
  #   int_var_name <- "Absolute of (1/B - e0)"
  # }
  

  
  
  #########################
  #######""" Graphics
  #############################
  
  #color vector for simulations
  mycol_sim <- rainbow(nsim)
  #color vector for enzymes
  mycol_E <- matrix(0, ncol=2, nrow=n)
  mycol_E <- data.frame(mycol_E)
  for (j in 1:n) {
    #column 1 = enzyme nam "Enzyme j"
    mycol_E[j,1] <- paste("Enzyme ",j)
    #column 2 = enzyme color
    mycol_E[j,2] <- j+1
  }
  
  
  if (show.plot==TRUE) {
    
    #RNV in relation to A
    if (new.window==TRUE) {dev.new()}
    plot(leg_A[which.sim,1:n],RNV_mean_size[near.RNV,1:n],xlab="Activities",ylab="RNV mean",type='n',cex.lab=cex.lab,...)
    for (s in which.sim) {
      #take rows such as col n+1=s (simul s) and col n+2=1 (RNV 1)
      #RNV_mean_size[,(n+1)]==s gives a vector of length(which.sim) with TRUE or FALSE
      #which give number in vector such as the two conditions are TRUE
      which.rows <- which(RNV_mean_size[,(n+1)]==s&RNV_mean_size[,(n+2)]==1)
      points(leg_A[s,1:n],RNV_mean_size[which.rows,1:n],col=mycol_sim[s],pch=16,...)
    }
    
    #RNV size in relation to RNV-ranking-order variable
    if (new.window==TRUE) {dev.new()}
    plot(rank_var_sim[,1:n],RNV_mean_size[near.RNV,1:n],xlab=rank_var_name,ylab="RNV mean",type='n',cex.lab=cex.lab,...)
    #plot(int_var[,1:n],RNV_mean_size[near.RNV,1:n],xlab=int_var_name,ylab="RNV mean",type='n')
    for (s in which.sim) {
      #take rows such as col n+1=s (simul s) and col n+2=1 (RNV 1), and take columns 1 to n
      which.rows <- which(RNV_mean_size[,(n+1)]==s&RNV_mean_size[,(n+2)]==1)
      points(rank_var_sim[which(which.sim==s),],RNV_mean_size[which.rows,1:n],col=mycol_sim[s],type="o",pch=16,...)
      #points(int_var[s,],RNV_mean_size[which.rows,1:n],col=mycol_sim[s])
    }
  }
  
  #### Correlation between RNV mean and RNV.ranking.order.factor
  lm_RNV <- lm(as.vector(RNV_mean_size[near.RNV,1:n])~as.vector(rank_var_sim))
  #lm_RNV <- lm(as.vector(RNV_mean_size[near.RNV,1:n])~as.vector(rank_var_sim))
  
  #summary(lm_RNV)
  if (add.lm==TRUE&show.plot==TRUE) {
    #add line of lm in plot
    abline(lm_RNV,col=1,lwd=2)
  }
  
  if (show.plot==TRUE) {
    
    #enough plac at right to add right axis for relative concentrations
    if(add.pred.e==TRUE) {par(mar=c(5,5,3,5))}
    
    #RNV size in relation to A, B or |1/B-e0| depending on applied constraint
    if (new.window==TRUE) {dev.new()}
    plot(int_var_sim[which.sim,1:n],RNV_mean_size[near.RNV,1:n],xlab=int_var_name,ylab="RNV mean",type='n',cex.lab=cex.lab,...)
    for (s in which.sim) {
      #take rows such as col n+1=s (simul s) and col n+2=1 (RNV 1), and take columns 1 to n
      which.rows <- which(RNV_mean_size[,(n+1)]==s&RNV_mean_size[,(n+2)]==1)
      points(int_var_sim[s,],RNV_mean_size[which.rows,1:n],col=mycol_sim[s],pch=16,...)
    }
    
    if (add.mean==TRUE) {
      #points of RNV mean for all selected simulation
      #points(int_var_mean,RNV_mean_by_enz,col='grey60',type="o",pch=15,lty=2,cex=1.5)
      points(RNV_and_var_mean_ord$int ,RNV_and_var_mean_ord$RNV,col=1,type="o",pch=15,lty=1,...)#cex=1.5)
    }
    #visualisation of 0 for B
    if (correl_fun=="RegNeg") {
      abline(v=0)
    }
    
    
    if (add.pred.e==TRUE) {
      #authorize superposition of plot
      par(new=T)
      #line predicted e fct A with points
      plot(RNV_and_var_mean_ord$int, RNV_and_var_mean_ord$e, pch=17, lty=2, xlab="", ylab="", ylim=c(0,1), axes=F, type="o", col="grey60",...)#,lwd=1.5, cex=1.5)
      #add right y-axis
      axis(4, ylim=c(0,1), col="grey60",col.axis="grey60",...)#,cex.axis=1.3)
      mtext("Relative concentrations",side=4,col="grey60",line=mar.lab,cex=cex.lab)
      
    }
  }
  
  
  
  return(invisible(list(RNV_all_sim=RNV_all_sim,RNV_mean_size=RNV_mean_size,rank_var_value=rank_var_sim,rank_var_name=rank_var_name,lm_RNV=lm_RNV)))
}

