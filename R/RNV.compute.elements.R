#' Compute RNV and associated elements for all enzymes \emph{(Deprecated)}
#'
#' \emph{(Deprecated)}. This function computes different elements at RNV, for different resident concentrations, based on function \code{\link{RNV.delta.all.enz}}
#' For simulations, preferably use \code{\link{RNV.for.simul}}.
#'
#' @details 
#' The \emph{Range of Neutral Variations} (RNV) are mutant concentration values such as coefficient selection is between \eqn{1/(2N)} and \emph{-1/(2N)}.
#' 
#' Inferior (resp. superior) bound of RNV corresponds to selection coefficient equal to \eqn{-1/(2N)} (resp. \eqn{1/(2N)}).
#'
#' 
#' Function \code{RNV.compute.elements} computes the \eqn{\delta_i} at RNV bounds (where \eqn{i} is the enzyme targeted by the mutation),
#' but also mutant concentrations at RNV bounds and the RNV size.
#' 
#' The RNV size is the absolute value of \eqn{\delta_i^sup} minus \eqn{\delta_i^inf}.
#'  
#' Depending on applied constraint \code{correl_fun}, it exists 1 or 2 RNV.
#' In case of independence (\code{"SC"}) or positive regulation between all enzymes (\code{"RegPos"}), flux has no limit and there is only one RNV.
#' In other cases (competition and/or negative regulation), because flux can reach a maximum, there is two RNVs:
#'  a "near" one, for small mutations, and a "far" one for big mutations that put mutant in the other side of flux dome.
#' 
#' This function \code{RNV.compute.elements} is designed to compute RNV of \emph{one} simulation launched by \code{\link{simul.evol.enz.multiple}}.
#' For example, for simulation 1, put \code{mat_E=tabR[tabR$sim==1,1:n]} and \code{mat_A=tabR[tabR$sim==1,(2*n+4):(3*n+3)]}. To reduce computation time, put \code{mat_J=tabR[tabR$sim==1,(2*n+3)]} for flux.
#' It also works for different resident values of concentrations and activities.
#' 
#' 
#' 
#' 
#' @usage RNV.compute.elements(mat_E,mat_A,N_fun,correl_fun,beta_fun=NULL,
#' end.mean=TRUE,add.RNV.J=FALSE,mat_J=NULL,X_fun=1)
#'
#' @inheritParams RNV.delta.all.enz
#' @param mat_E Numeric matrix of concentrations. Columns correspond to enzyme number, and rows correspond to different resident.
#' @param mat_A Numeric matrix of activities. Same dimensions as \code{mat_E}.
#' @param end.mean Logical. If \code{FALSE}, compute RNV size mean of all resident. If \code{TRUE}, compute RNV size mean for last half of resident only.
#' @param add.RNV.J Logical. Add value of flux at RNV ? If \code{TRUE}, possibility set also \code{mat_J}.
#' @inheritParams flux
#' @param mat_J Numeric matrix of flux. One column and same rows number as \code{mat_E}. Optional.
#'
#'
#' @return List of 7 elements:
#' \itemize{
#'    \item \code{$RNV_delta} : list of \code{n} elements, one by enzyme.
#'    Each element contains a numeric matrix of \code{nb_resid} rows and two or four columns.
#'    For each enzyme \eqn{i} (a list element), each row corresponds to one resident and columns are actual mutation effect \eqn{\delta_i} corresponding respectively to inferior bound (col 1) and superior bound (col 2) of RNV (x2 if there is a second RNV).
#'    Each row is a result of \code{\link{RNV.delta.all.enz}}.
#'    \item \code{$RNV_enz} : same structure as \code{$RNV_delta}, but for mutant enzyme concentrations at RNV limits, \emph{i.e.} \eqn{E_i + \delta_i} for each target enzyme \eqn{i} ;
#'    \item \code{$RNV_size} : list of one or two elements (depending on RNV number).
#'    Each element contains a matrix of \code{n} columns (one by enzyme) and \code{nb_resid} rows.
#'    Each cell is the RNV size for current enzyme (in column) and current resident (in row).
#'    The RNV size is absolute value of \eqn{\delta_i^sup} minus \eqn{\delta_i^inf}. If there is no superior bounds but there is two RNV, RNV size is obtained by the difference of the two \eqn{\delta_i^inf}.
#'    \item \code{$RNV_size_divEtot} : same structure as \code{$RNV_size}, but for RNV size divided by total concentration of corresponding resident.
#'    \item \code{$RNV_proxy} : numeric matrix of one or two rows (depending on RNV number) and \code{n} columns (one by enzyme).
#'    Each cell id the mean of RNV size. If \code{end.mean=TRUE}, the mean is computed from last half of resident.
#'    \item \code{$RNV_proxy_divEtot} : same structure as \code{$RNV_proxy}, but mean of RNV size is divided by mean of total concentration.
#'    \item \code{$RNV_flux} : numeric matrix of two columns (inferior and superior limits of neutral zone) and \code{nb_resid} rows.
#'    Each cell is the flux value at neutral zone limits.
#' }
#' 
#' 
#' 
#' Note that \code{n} is the number of enzymes. \code{nb_resid} is the number of resident and s also the row number of \code{mat_E} and \code{mat_A}.
#' 
#' 
#' @seealso 
#' See function \code{\link{RNV.delta.all.enz}} to see how actual mutation effect \eqn{\delta} is computed. 
#'
#' For simulations, preferably use \code{\link{RNV.for.simul}}. Code is almost the same for the two functions, but differs in input.
#'
#' @examples
#' \donttest{
#' #for 2 resident genotypes and 3 enzymes
#' Er <- matrix(30,ncol=3,nrow=2)
#' A <- matrix(c(1,10,30),byrow=TRUE,ncol=3,nrow=2)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- compute.B.from.beta(beta)
#' correl <- "CRPos"
#' N <- 1000
#' 
#' #second resident = theoretical equilibrium of first resident
#' Er[2,] <- 100*predict_th(A[1,],correl,B)$pred_e
#'
#' RNV.compute.elements(Er,A,N,correl,beta,add.RNV.J=TRUE)
#' }
#'
#'
#'
#' @export






RNV.compute.elements <- function(mat_E,mat_A,N_fun,correl_fun,beta_fun=NULL,end.mean=TRUE,add.RNV.J=FALSE,mat_J=NULL,X_fun=1) {
  
  ########### Tests
  if (nrow(mat_E)!=nrow(mat_A)) {
    stop("Resident is missing. Matrix 'mat_E' and 'mat_A' should have same number of row.")
  }
  if (ncol(mat_E)!=ncol(mat_A)) {
    stop("An enzyme is missing. Matrix 'mat_E' and 'mat_A' should have same number of column.")
  }
  if (add.RNV.J==TRUE&length(mat_J)!=0&length(mat_J)!=nrow(mat_E)) {
    stop("For neutral zone of flux, please put same number of resident for 'mat_E' and 'mat_J' or if flux is unknown, set 'mat_J' to NULL.")
  }

  
  
  ######### Parameters
  #Number of enzymes
  n <- ncol(mat_E)
  #Number of resident
  nb_resid <- nrow(mat_E)
  #limits of Neutral Zone in terms of selection coefficient
  limits_NZ <- c(-1/(2*N_fun),1/(2*N_fun))  #c("inf","sup")
  #total concentrations
  mat_Etot <- apply(mat_E,1,sum)
  
  #how many RNV depending on constraints
  if (correl_fun=="SC"|correl_fun=="RegPos") {
    nb_RNV <- 1
  } else {
    nb_RNV <- 2
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
  for (r in nb_resid) {
    
    #if flux is unknown
    if (length(mat_J)==0) {
      #extract resident concentrations for this line
      E_resid <- as.matrix(mat_E[r,1:n])
      A_resid <- as.matrix(mat_A[r,1:n])
      
      #compute resident flux
      J_resid <- flux(E_resid,A_resid,X_fun)
    } else { #if flux is already computed and input in parameters
      J_resid <- mat_J[r]
    }
    
    #compute flux at NZ limits
    RNV_dis_flux[r,] <- J_resid*(1+limits_NZ)
  }
  
  
  return(list(RNV_delta=RNV_dis_delta_graph, RNV_enz=RNV_dis_enz_graph, RNV_size=RNV_dis_size,
              RNV_size_divEtot=RNV_dis_size_divEtot, RNV_proxy=RNV_proxy, RNV_proxy_divEtot=RNV_proxy_divEtot,
              RNV_flux=RNV_dis_flux))
}

