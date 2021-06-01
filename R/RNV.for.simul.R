#' Compute RNV for simulation
#'
#' Computes different elements at RNV for evolution simulation
#' 
#'
#' @details 
#' This function is designed to computes RNV of \emph{one} simulation launched by \code{\link{simul.evol.enz.one}}.
#' Input \code{res_sim} is a result of \code{\link{simul.evol.enz.one}}.
#' If you used \code{\link{simul.evol.enz.multiple}}, input \code{tabR} for parameter \code{res_sim}.
#' If input is several simulations, remember that number of rows for one simulation is \code{npt}. 
#' \emph{See below.}
#'   
#' The \emph{Range of Neutral Variations} (RNV) are mutant concentration values such as coefficient selection is between \eqn{1/(2N)} and \eqn{-1/(2N)}.
#' 
#' Inferior (resp. superior) bound of RNV corresponds to selection coefficient equal to \eqn{-1/(2N)} (resp. \eqn{1/(2N)}).
#' 
#' Function \code{RNV.for.simul} computes the actual mutation effect \eqn{\delta_i} at RNV bounds (where \eqn{i} is the enzyme targeted by the mutation),
#' but also mutant concentration at RNV bounds and the RNV size.
#' 
#' 
#' 
#' Depending on applied constraint \code{correl_fun}, it exists 1 or 2 RNV.
#' In case of independence (\code{"SC"}) or positive regulation between all enzymes (\code{"RegPos"}), flux has no limit and there is only one RNV.
#' In other cases (competition and/or negative regulation), because flux can reach a maximum, there is two RNV: 
#' a "near" one, for small mutations, and a "far" one for big mutations that put mutants on the other side of flux dome.
#' 
#' 
#' 
#' \bold{RNV size}
#' 
#' The RNV size is the absolute value of \eqn{\delta_i^sup} minus \eqn{\delta_i^inf}.
#' If there is no superior bounds but there is two RNVs, RNV size is obtained by the difference of the two \eqn{\delta_i^inf}.
#'  
#' The mean of RNv size is the mean of every resident for each enzyme. If \code{end.mean=TRUE}, only last half of resident (the last half of \code{res_sim} rows) are used to compute RNV mean.
#'  
#'  
#' \bold{Use of \code{res_sim}}
#' 
#' \code{res_sim} is a numeric matrix of \code{(3*n+4)} columns and at least 2 rows.
#' Respective columns are: concentrations (\code{1:n}), kinetic parameters (\code{n+1:2n}), total concentration (\code{2n+1}), total kinetic (\code{2n+2}), flux/fitness (\code{2n+3}) and activities (\code{2n+4:3n+3}), corresponding to \emph{(E1 to En, kin_1 to kin_n, Etot, kin_tot, J, A1 to An)}.
#' \emph{See function \code{\link{simul.evol.enz.one}}.}
#' 
#' \code{res_sim} is normally output \code{$res_sim} of function \code{\link{simul.evol.enz.one}}.
#' Output \code{$tabR} of function \code{\link{simul.evol.enz.multiple}} is also possible, by selecting a simulation with \code{x$tabR[x$tabR$sim==i,]} where \code{i} is simulation number.
#' If input is several simulations, remember that number of rows for one simulation is \code{npt}. 
#' 
#' Other parameters (\code{n_fun, N_fun, correl_fun, beta_fun}) are available in output \code{$param} of simulation functions.
#' 
#' 
#' @usage RNV.for.simul(res_sim,n_fun,N_fun,correl_fun,beta_fun=NULL,
#' end.mean=TRUE)
#'
#' @inheritParams RNV.delta.all.enz
#' @param res_sim Dataframe corresponding to result of one simulation. \emph{See details}.
#' @param n_fun Integer number indicating the number of enzymes. 
#' @param end.mean Logical. If \code{FALSE}, compute RNV size mean for all rows of \code{res_sim}. If \code{TRUE}, compute RNV size mean for last half of \code{res_sim} rows.
#'
#'
#' @return Invisible list of 7 elements:
#' \itemize{
#'    \item \code{$RNV_delta}: list of \code{n} elements, one by enzyme.
#'    Each element contains a numeric matrix of \code{nb_resid} rows and two or four columns.
#'    For each enzyme \code{i} (a list element), each row corresponds to one resident and columns are actual mutation size \eqn{\delta_i} corresponding respectively to inferior bound (col 1) and superior bound (col 2) of RNV (x2 if there is a second RNV).
#'    Each row is in fact a result of \code{\link{RNV.delta.all.enz}}.
#'    \item \code{$RNV_enz}: same structure as \code{$RNV_delta}, but for mutant enzyme concentrations at RNV limits, \emph{i.e.} \eqn{E_i + \delta_i} for each target enzyme \code{i}.
#'    \item \code{$RNV_size}: list of one or two elements (depending RNV on RNV number).
#'    Each element contains a matrix of \code{n} columns (one by enzyme) and \code{nb_resid} rows.
#'    Each cell is the RNV size \emph{(see details)}.
#'    \item \code{$RNV_size_divEtot}: same structure as \code{$RNV_size}, but for RNV size divided by total concentration of corresponding resident.
#'    \item \code{$RNV_proxy}: numeric matrix of one or two rows (depending on RNV number) and \code{n} columns.
#'    Each cell is the mean of RNV size \emph{(see details)}.
#'    \item \code{$RNV_proxy_divEtot}: same structure as \code{$RNV_proxy}, and contains mean of RNV size divided by total concentration.
#'    \item \code{$RNV_flux}: numeric matrix of two columns (inferior and superior limits of neutral zone) and \code{nb_resid} rows.
#'    Each cell is the flux value at neutral zone limits.
#'    \item \code{$nb_RNV}: number of RNV. If constraint is \code{"SC"} or \code{"RegPos"}, there is one RNV, else two.
#'    \item \code{$limits_NZ}: numeric vector of the two limits of neutral zone
#' }
#' 
#' 
#' 
#' Note that \code{n} is the number of enzymes. \code{nb_resid} is the number of resident and is also the rows number of \code{res_sim}.
#' 
#' 
#' @seealso 
#' See function \code{\link{RNV.delta.all.enz}} to see how \eqn{\delta} is computed. 
#' 
#' Use function \code{\link{simul.evol.enz.one}} to launch a simulation, or \code{\link{simul.evol.enz.multiple}} for several simulations.
#'
#'
#' @examples
#' 
#'  #### Construction of false simulation
#' #for 2 resident genotypes and 3 enzymes
#' n <- 3
#' Er <- c(30,30,30)
#' kin <- c(1,10,30)
#' Keq <- c(1,1,1)
#' A <- activities(kin,Keq)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- compute.B.from.beta(beta)
#' correl <- "RegPos"
#' N <- 1000
#' 
#' #on one line
#' first_res <- cbind(t(Er),t(kin),sum(Er),sum(kin),flux(Er,A),t(A))
#' 
#' #second resident = theoretical equilibrium of first resident
#' Eq_th <- 100*predict_th(A,correl,B)$pred_e
#' second_res <- cbind(t(Eq_th),t(kin),sum(Eq_th),sum(kin),flux(Eq_th,A),t(A))
#' 
#' false_sim <- rbind(first_res,second_res)
#'
#' RNV_elements <- RNV.for.simul(false_sim,n,N,correl,beta)
#'
#' RNV_elements$RNV_delta #apparent mutation size at RNV for enzymes 1, 2 and 3
#' RNV_elements$RNV_enz #concentrations
#' RNV_elements$RNV_size #RNV size
#' RNV_elements$RNV_proxy #RNV size mean
#' RNV_elements$RNV_flux #flux at neutral zone
#' 
#' \donttest{
#' #With saved simulation
#' data(data_sim_RegPos)
#' RNV_elements <- RNV.for.simul(data_sim_RegPos$tabR,data_sim_RegPos$param$n,
#' data_sim_RegPos$param$N,data_sim_RegPos$param$correl,data_sim_RegPos$param$beta)
#' }
#'
#'
#' @export





RNV.for.simul <- function(res_sim,n_fun,N_fun,correl_fun,beta_fun=NULL,end.mean=TRUE) {
  

  ######### Parameters
  #Numnber of enzymes
  # n <- (ncol(res_sim)-3)/3
  n <- n_fun
  
  #Number of resident
  nb_resid <- nrow(res_sim)
  
  #limits of Neutral Zone in terms of selection coefficient
  limits_NZ <- c(-1/(2*N_fun),1/(2*N_fun))  #c("inf","sup")
  #limits_NZ <- limits_NZ[order(limits_NZ)] #for future improvment, if choice of neutral zone, add this railway line
  
  #how many RNV depending on constraints
  if (correl_fun=="SC"|correl_fun=="RegPos") {
    nb_RNV <- 1
  } else {
    nb_RNV <- 2
  }
  
  ########## Extract results
  #Take concentrations and total concentration from simulation
  #mat_E <- res_sim[,c(1:n,2*n+1)]
  #mat_E <- cbind(res_sim[,1:n],res_sim[,(2*n+1)])
  mat_E <- res_sim[,1:n]
  mat_Etot <- res_sim[,(2*n+1)]
  #tabE <- cbind(tabE,tabR[tabR$sim==s,(2*n+1)]) #total concentration
  #Take activities from simulation
  mat_A <- res_sim[,(2*n+4):(3*n+3)]
  #Take flux from simulation
  mat_J <- res_sim[,(2*n+3)]
  
  flag_length <- 0
  #if input a vector as res_sim or as matrix but with only 1 row
  #take only first value of class, because class(matrix) gives a vector containing c("matrix","array")
  if (class(res_sim)[1]=="numeric") { flag_length <- 1  }
  if (class(res_sim)[1]=="matrix"|class(res_sim)[1]=="data.frame") {if (nrow(res_sim)==1) {flag_length <- 1}}
  if (flag_length==1) {
    #because in this case, vector is consider as one column, and the following code is adapted to take elements by the different columns
    mat_E <- as.matrix(t(mat_E))
    mat_A <- as.matrix(t(mat_A))
  }
  
  
  #########" Computes RNV
  
  #list of n elements, one by enzyme. For each enzyme, rows = simulation and columns = delta corresponding to inf bounds and sup bounds of RNV (x2 if there is a 2nd RNV)
  RNV_dis_delta_graph <- vector("list",length=n)
  #idem, but for enzyme concentrations + delta
  RNV_dis_enz_graph <- vector("list",length=n)
  
  #for each result of current simulation s
  for (r in 1:nb_resid) {
    #extract resident concentrations for this line
    E_resid <- as.matrix(mat_E[r,1:n])
    A_resid <- as.matrix(mat_A[r,1:n])
    
    # compute delta at limits of neutral zone
    delta_RNV_all <- RNV.delta.all.enz(E_resid,A_resid,N_fun,correl_fun,beta_fun)
    #delta_all <- delta.RNV.all(E_resid,A_base,N,correl,beta,B,n,lim_calc) #E_fun,A_fun,N_fun,correl_fun,beta_fun,B_fun,n_fun,d_lim
    
    #for each enzyme
    for (i in 1:n) {
      #take corresponding delta at bounds of RNV
      RNV_dis_delta_graph[[i]] <- rbind(RNV_dis_delta_graph[[i]],delta_RNV_all[[i]])
      #add delta to concentrations to have corresponding "mutant" concentration at bounds of RNV
      RNV_dis_enz_graph[[i]] <- rbind(RNV_dis_enz_graph[[i]],E_resid[i] + delta_RNV_all[[i]])
      
    }
  }
  
  
  ####### RNV size
  ##list of one or two elemnts (number of RNV), and in each element, matrix of n columns and nb_resid rows (nb points for current simulation)
  #size of RNV (delta_sup - delta_inf)
  RNV_dis_size <- vector("list",length=nb_RNV)
  #RNV size divided by resident total concentration
  RNV_dis_size_divEtot <- vector("list",length=nb_RNV)
  
  ## matrix of one or two rows (number of RNV) and n columns
  #mean RNV size
  RNV_proxy <- matrix(0,nrow=nb_RNV,ncol=n)
  #mean RNV size divided by resident total concentration
  RNV_proxy_divEtot <- matrix(0,nrow=nb_RNV,ncol=n)
  
  #for each RNV
  for (nor in 1:nb_RNV) {
    #for each enzyme
    for (i in 1:n) {
      #take RNV size for all result 'r'
      RNV_size_vec_r <- NULL
      
      #for each result
      for (r in 1:nb_resid) {
        
        #superior and inferior bounds of current RNV 'nor' and current enzyme 'i' and current result 'r'
        RNV_sup_bounds_delta <- RNV_dis_delta_graph[[i]][r,2*nor]
        RNV_inf_bounds_delta <- RNV_dis_delta_graph[[i]][r,2*nor-1]
        
        #if there is two RNV and sup bounds is not accessible
        if (nb_RNV==2 & is.na(RNV_sup_bounds_delta)) {
          # RNV size = | inf bounds (RNV 2) - inf bounds (RNV 1) |
          RNV_size_num_r <- abs(RNV_dis_delta_graph[[i]][r,3] - RNV_dis_delta_graph[[i]][r,1])
          #RNV_dis_size[[nor]] <- cbind(RNV_dis_size[[nor]], abs(RNV_dis_delta_graph[[i]][,3] - RNV_dis_delta_graph[[i]][,1]))
        } else {
          #RNV size = |sup bounds - inf bounds| for RNV 'nor'
          RNV_size_num_r <- abs(RNV_sup_bounds_delta - RNV_inf_bounds_delta)
          #RNV_dis_size[[nor]] <- cbind(RNV_dis_size[[nor]], abs(RNV_sup_bounds_delta - RNV_inf_bounds_delta))
        }
        
        #add this RNV size for this result 'r' with precedent
        RNV_size_vec_r <- c(RNV_size_vec_r,RNV_size_num_r)
      }
      
      #add RNV size for all rows with precedent enzyme 
      RNV_dis_size[[nor]] <- cbind(RNV_dis_size[[nor]], RNV_size_vec_r)
    }
    
    #division by Etot : mat_E has same nrows than RNV_dis_size element, and R compute column by column
    RNV_dis_size_divEtot[[nor]] <- RNV_dis_size[[nor]]/mat_Etot #mat_E[,n+1]
    
    #mean of RNV size for half end of simul
    if (end.mean==TRUE) {
      RNV_proxy[nor,] <- apply(RNV_dis_size[[nor]][round(nb_resid/2):nb_resid,],2,mean,na.rm=TRUE)
      RNV_proxy_divEtot[nor,] <- apply(RNV_dis_size_divEtot[[nor]][round(nb_resid/2):nb_resid,],2,mean,na.rm=TRUE)
    } else { #mean for all resident
      RNV_proxy[nor,] <- apply(RNV_dis_size[[nor]],2,mean,na.rm=TRUE)
      RNV_proxy_divEtot[nor,] <- apply(RNV_dis_size_divEtot[[nor]],2,mean,na.rm=TRUE)
    }
  }
  #print(RNV_proxy)
  
  # delta_mean <- vector("list",length=n)
  # for (i in 1:n) {
  #   bfff <- apply(RNV_dis_delta_graph[[i]][round(nb_resid/2):nb_resid,],2,mean,na.rm=TRUE)
  #   delta_mean[[i]] <- c(delta_mean[[i]],bfff)
  # }
  
  
  ############## RNV flux
  #for flux, there is only neutral zone, and nb of RNV does not matter
  #matrx of 2 columns (inferior and superior bounds of neutral zone) and nb_resid rows for each resident value
  RNV_dis_flux <- matrix(0,ncol=2,nrow=nb_resid)
  for (r in 1:nb_resid) {
    #take flux value for this resident
    J_resid <- mat_J[r]
    
    #compute flux at NZ limits
    RNV_dis_flux[r,] <- J_resid*(1+limits_NZ)
    
    #if flux at sup RNV superior to max flux (if exists)
    # if (RNV_dis_flux[r,2]>max(mat_J)&correl_fun!="SC"&correl_fun!="RegPos") {
    #   #replaces by NA, because no sup limits for enzyme
    #   RNV_dis_flux[r,2] <- NA
    # }
  }
  
  
  return(invisible(list(RNV_delta=RNV_dis_delta_graph, RNV_enz=RNV_dis_enz_graph, RNV_size=RNV_dis_size,
              RNV_size_divEtot=RNV_dis_size_divEtot, RNV_proxy=RNV_proxy, RNV_proxy_divEtot=RNV_proxy_divEtot,
              RNV_flux=RNV_dis_flux,nb_RNV=nb_RNV,limits_NZ=limits_NZ)))
}

