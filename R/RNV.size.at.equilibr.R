#' Plot and linear model of RNV size at equilibrium
#'
#' Gives a plot and computes a linear model of the RNV size at equilibrium against the activities and the RNV-ranking-order factor
#' 
#' 
#'
#' @details 
#' Function \code{RNV.size.at.equilibr} gives a plot of RNV size in relation to activities, and then a plot of RNV size in relation to RNV-ranking-order-factor.
#' 
#' RNV size is computed at equilibrium with function \code{\link{RNV.for.simul}}. Only "near" RNV (small mutations) is used.
#' 
#' \bold{WARNING!} This function is not adapted or regulation groups (\code{1<sum(1/B)<n}).
#' 
#' \bold{Factors governing ranking order of RNV}
#' 
#' RNV-ranking-order factor depends on constraint. See function \code{\link{RNV.ranking.order.factor}}.
#' 
#' 
#' \bold{Random variables}
#' 
#' Input parameters \code{A_fun, B_fun, E0_fun} are chosen randomly if not specified (default \code{NA}).
#' 
#' Activities \code{A_fun} are taken in a uniform law between limits given by \code{A_lim}.
#' 
#' Initial concentrations \code{E0_fun} are taken in a uniform law between limits given by \code{E0_lim}, then leveled to have sum of \code{E0_fun} equal to \code{Etot_fun}.
#' 
#' Global co-regulation coefficients \code{B_fun} are chosen differently. The inverse of \code{B} are taken between limits given by \code{inv_B_lim}.
#' If \code{correl_fun} is equal to \code{"RegNeg"} or \code{"CRNeg"} (negative co-regulations), sign are randomly chosen.
#' Then these inverse values are leveled to have sum of \code{1/B} equal to 1.
#' Thus \code{B} are computed by the reverse operation, and therefore the matrix of \code{beta} by function \code{\link{compute.beta.from.B}}.
#' 
#' 
#' \bold{Equilibrium}
#' 
#' RNV is computed at equilibrium, theoretical one for cases \code{"SC"}, \code{"Comp"} and \code{"RegPos"}, and effective one for cases \code{"RegNeg"}, \code{"CRPos"} and \code{"CRNeg"}.
#' 
#' 
#' \bold{Graphics}
#' 
#' If \code{RNV.by.E=TRUE}, RNV size divided by its corresponding enzyme concentration is plotted, rather than RNV size.
#' 
#' 
#' @usage RNV.size.at.equilibr(n_fun,correl_fun,N_fun,A_fun=NA,B_fun=NA,E0_fun=NA,
#' Etot_0=sum(E0_fun),Etot_eq=100,A_lim=c(0.001,10000),inv_B_lim=c(0,100),E0_lim=c(1,100),
#' RNV.by.E=FALSE,show.plot=TRUE,enz.lab=FALSE,...)
#'
#' @inheritParams RNV.for.simul
#' @param A_fun,B_fun,E0_fun Numeric vectors (length \code{n_fun}) of activities, global co-regulation coefficients and initial concentrations respectively.
#' If \code{NA} (default), values are randomly chosen.
#' @param Etot_0 Numeric value of initial total concentration. \code{E0_fun} is rescaled to have \code{sum(E0_fun)=Etot_0}.
#' @param Etot_eq Numeric value of total concentration at equilibrium for cases \code{"SC"} and \code{"RegPos"}. Default is \code{100}.
#' @param A_lim Numeric vector of the two limits between which activities are chosen if \code{A_fun=NA}.
#' @param inv_B_lim Numeric vector of the two limits between which inverse of global co-regulations \code{B} are chosen if \code{B_fun=NA}.
#' @param E0_lim Numeric vector of the two limits between which initial concentrations \code{E0} are chosen if \code{E0_fun=NA}.
#' @param show.plot Logical. Show plot of RNV in relation to RNV-ranking-order factor?
#' @param ... Arguments to be passed to \code{plot} function, such \code{cex} or \code{cex.lab}
#' @param enz.lab Logical. Add enzyme name above points in graphics?
#' @param RNV.by.E Logical. Show RNV divided by Ei rather than RNV size?
#' 
#'
#' @importFrom stats lm
#'
#' @return Invisible list of 12 elements:
#' \itemize{
#'    \item \code{$RNV_size}: matrix of one or two rows and \code{n} columns indicating the RNV size of every enzymes (in columns) and current RNV (near or far, in rows). See output \code{$RNV_size} of function \code{\link{RNV.for.simul}};
#'    \item \code{$lm_RNV}: linear model of RNV size in relation to the ranking-order variable;
#'    \item \code{$E_eq}: numeric vector (length \code{n}) of concentrations at equilibrium;
#'    \item \code{$e_eq}: numeric vector (length \code{n}) of relative concentrations at equilibrium;
#'    \item \code{$A}: numeric vector (length \code{n}) of activities;
#'    \item \code{$B}: numeric vector (length \code{n}) of global co-regulation coefficients;
#'    \item \code{$beta}: numeric matrix of \code{n} rows and \code{n} columns indicating the co-regulation coefficients;
#'    \item \code{$E0}: numeric vector (length \code{n}) of initial concentrations;
#'    \item \code{$rank_var}: list of 2 elements: \itemize{
#'       \item \code{$value}: numeric vector (length \code{n}) of the RNV-ranking-order factor values for each enzyme \emph{(see details)};
#'       \item \code{$name}: character string indicating the name of the RNV-ranking-order factor.
#'    }
#'    \item \code{$N}, \code{$n}, \code{$correl}: numeric value of population size, number enzymes and applied constraints respectively (respectively input value of \code{N_fun}, \code{n_fun} and \code{correl_fun})
#' }
#' 
#' 
#' @seealso 
#' RNV is computed with the function \code{\link{RNV.for.simul}}. 
#' 
#' See function \code{\link{RNV.ranking.order.factor}} to have further details on RNV-ranking-order.
#'
#'
#' @examples
#' RNV.size.at.equilibr(20,"Comp",1000)
#' RNV.size.at.equilibr(100,"CRNeg",1000)
#' RNV.size.at.equilibr(3,"SC",1000,c(1,10,30),NA,c(30,30,30),100,200)
#'  #in this case, sum(E0)=100 and sum(E*)=200
#'
#'
#'
#' @export





RNV.size.at.equilibr <- function(n_fun,correl_fun,N_fun,A_fun=NA,B_fun=NA,E0_fun=NA,
                                           Etot_0=sum(E0_fun),Etot_eq=100,A_lim=c(0.001,10000),inv_B_lim=c(0,100),E0_lim=c(1,100),
                                           RNV.by.E=FALSE,show.plot=TRUE,enz.lab=FALSE,...) {
  
  ### Test
  is.correl.authorized(correl_fun)
  #for all specified parameters, length verification
  if (!all(is.na(A_fun))) {if (length(A_fun)!=n_fun) {stop("A_fun and n_fun need to have same length.")}}
  if (!all(is.na(E0_fun))) {if (length(E0_fun)!=n_fun) {stop("E0_fun and n_fun need to have same length.")}}
  if (!all(is.na(B_fun))) {is.B.accurate(B_fun,n_fun,correl_fun)}
  
  ### Random parameters if not specified
  
  #activities
  if (all(is.na(A_fun))) {
    A_fun <- runif(n_fun, min(A_lim), max(A_lim))
  }

  #initial concentrations
  if (all(is.na(E0_fun))) {
    E0_fun <- runif(n_fun,min(E0_lim),max(E0_lim))
  }

  
  #co-regulation coefficient
  if (all(is.na(B_fun))) {
    
    if (correl_fun=="SC"|correl_fun=="Comp") {
      B_fun <- rep(1,n_fun)
      beta_fun <- diag(1,n_fun)
    }
    if (correl_fun=="RegPos"|correl_fun=="CRPos"|correl_fun=="RegNeg"|correl_fun=="CRNeg") {
      #easier to choose randomly 1/B rather than B, because sum(1/B)=1
      
      #inv_B <- c(0.5,0.1,0.4)
      inv_B <- runif(n_fun,min(inv_B_lim),max(inv_B_lim))
      
      if (correl_fun=="RegNeg"|correl_fun=="CRNeg") {
        #inv_B <- c(0.5,-0.2,0.7)
        #sign randomly chosen
        has <- runif(n_fun,min=0,max=1)
        sign <- 1*(has<=0.5)-1*(has>0.5)
        inv_B <- sign*inv_B
      }
      
      #to have sum(1/B)=1
      inv_B <- inv_B/sum(inv_B)
      B_fun <- 1/inv_B

      
    } #end of co-regulation coefficient
  }
  
  ##### Test
  
  
  ########## Other parametrization
  #if Etot is specified
  if (!is.na(Etot_0)) {
    #cross product to have sum(E0)=Etot
    E0_fun <- Etot_0*E0_fun/sum(E0_fun)
  }
  
  #relative initial concentrations
  e0_fun <- E0_fun/sum(E0_fun)
  
  #compute beta
  beta_fun <- compute.beta.from.B(B_fun)
  #beta <- B%*%t(1/B) #multiplication matricielle
  
  
  ##### Compute equilibrium
  
  #theoretical
  eq_th <- predict_th(A_fun,correl_fun,B_fun)
  #effective
  if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
    eq_eff <- predict_eff(E0_fun,B_fun,A_fun,correl_fun)
  } 
  #choose equilibrium as resident to compute RNV at equilibrium
  if (correl_fun=="SC"|correl_fun=="Comp"|correl_fun=="RegPos") {
    #th eq
    e_eq <- eq_th$pred_e
    if (correl_fun=="Comp") {
      #fixed Etot
      E_eq <- eq_th$pred_e*sum(E0_fun)
    } else {
      #chosen Etot at eq
      E_eq <- eq_th$pred_e*sum(Etot_eq)
    }

  } else { #effective eq
    E_eq <- eq_eff$pred_E
    e_eq <- eq_eff$pred_e
  }
  
  ##### Compute RNV
  #delta_all <- RNV.delta.all.enz(E_eq,A_fun,N_fun,correl_fun,beta_fun)
  
  #format to correspond to input res_sim format
  false_sim <- cbind(t(E_eq),t(A_fun),sum(E_eq),sum(A_fun),flux(E_eq,A_fun),t(A_fun))
  #NB : flux will be not use after, so X_fun is not useful
  
  RNV_all <- RNV.for.simul(false_sim,n_fun,N_fun,correl_fun,beta_fun,end.mean=FALSE)
  
  #recup mean of RNV size
  #because mean for only one value is equal to this value, it is easier to use RNV_proxy (matrix) rather than RNV_size (list)
  RNV_size <- RNV_all$RNV_proxy
  
  
  
  ##### which RNV is plotted
  #Only nearest RNV is plotted
  RNV_plotted <- RNV_size[1,]
  RNV_lab <- "RNV size"
  if (RNV.by.E==TRUE) {
    #choice to plot RNV divided by Ei
    RNV_plotted <- RNV_size[1,]/E_eq
    RNV_lab <- "RNV size div by E"
  }
  
  
  
  
  ##### Compare RNV_size to RNV-ranking-order variable depending on constraints
  #plot and linear model of RNV size depending on RNV-ranking-order variable 
  
  rank_var <- RNV.ranking.order.factor(A_fun,correl_fun,E0_fun,B_fun)
  
  #### Enzyme name
  enz_name <- NULL
  for (j in 1:n_fun) {enz_name[j] <- paste0("E",j)}
  
  ### Group data
  compar_at_eq <- cbind.data.frame(A=A_fun,RNV=RNV_plotted,e=e_eq,B=B_fun,rank_var=rank_var$value,name=enz_name)
  
  ## Graphics
  if (show.plot==TRUE) {
    #par(mar=c(5,5,2,3))
    
    #RNV fct A
    plot(A_fun,RNV_plotted,xlab="Pseudo-activity A",ylab=RNV_lab,pch=16,...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3,lwd=1.5)
    if (enz.lab==TRUE) {
      # add name of enzyme above points
      text(A_fun,RNV_plotted,labels=enz_name,pos=3)
    }
    
    #RNV fct A well-presented
    if (correl_fun=="SC"|correl_fun=="Comp") {
      #ordered depending on A to easily rely points
      compar_at_eq_ord <- compar_at_eq[order(compar_at_eq$A),]
      plot(compar_at_eq_ord$A,compar_at_eq_ord$RNV,xlab="Pseudo-activity A",ylab=RNV_lab,pch=15,type="o",...)
      if (enz.lab==TRUE) {
        text(compar_at_eq_ord$A,compar_at_eq_ord$RNV,labels=compar_at_eq_ord$name,pos=3)
      }
    }
    
    #RNV fct B
    if (correl_fun=="RegPos"|correl_fun=="RegNeg") {
      #ordered depending on B to easily rely points
      compar_at_eq_ord <- compar_at_eq[order(compar_at_eq$B),]
      plot(compar_at_eq_ord$B,compar_at_eq_ord$RNV,xlab="Global co-regulation coefficient B",ylab=RNV_lab,pch=15,type="n",...)
      #separation between negative and positive B
      points(compar_at_eq_ord$B[which(compar_at_eq_ord$B<0)],compar_at_eq_ord$RNV[which(compar_at_eq_ord$B<0)],type="o",pch=15,...)
      points(compar_at_eq_ord$B[which(compar_at_eq_ord$B>=0)],compar_at_eq_ord$RNV[which(compar_at_eq_ord$B>=0)],type="o",pch=15,...)
      #show symmetry 
      abline(v=0,lwd=1.5)
      # add name of enzyme above points
      if (enz.lab==TRUE) {
        text(compar_at_eq_ord$B,compar_at_eq_ord$RNV,labels=compar_at_eq_ord$name,pos=3)
      }
    }
        
    #RNV fct |1/B| si interessant
    if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
      plot(abs(1/B_fun),RNV_plotted,xlab="Absolute value of 1/B",ylab=RNV_lab,pch=16,...)
      # add name of enzyme above points
      if (enz.lab==TRUE) {
        text(abs(1/B_fun),RNV_plotted,labels=enz_name,pos=3)
      }
    }
    
    #RNV fct ranking variable
    #ordered depending on ranking-order-factor to easily rely points
    compar_at_eq_ord <- compar_at_eq[order(compar_at_eq$rank_var),]
    plot(compar_at_eq_ord$rank_var,compar_at_eq_ord$RNV,xlab=rank_var$name,ylab=RNV_lab,pch=15,type="o",...)
    #plot(rank_var$value,RNV_plotted,xlab=rank_var$name,ylab=RNV_lab,pch=15,type="o",...)
    # add name of enzyme above points
    if (enz.lab==TRUE) {
      text(compar_at_eq_ord$rank_var,compar_at_eq_ord$RNV,labels=compar_at_eq_ord$name,pos=3)
    }
  }

  ## Linear model
  lm_RNV <- lm(as.vector(RNV_size[1,])~rank_var$value)
  summary(lm_RNV)
  
  
  
  
  # if (correl_fun=="SC") {
  #   #int_var <- A_fun^(1/3)
  #   plot(int_var,RNV_size[1,],xlab="Activities^(1/3)",ylab="RNV size")
  #   #my_lm <- lm(as.vector(RNV_size[1,])~int_var)
  #   #summary(lm(as.vector(RNV_size[1,])~int_var))
  #   # plot(-log(A,base=1/3),RNV_size[1,])
  #   # plot(-log(A,base=1/2),RNV_size[1,])
  # }
  # 
  # if (correl_fun=="Comp") {
  #   int_var <- A_fun^(-1/4)
  #   plot(int_var,RNV_size[1,],xlab="Activities^(-1/4)",ylab="RNV size")
  #   #summary(lm(as.vector(RNV_size[1,])~int_var))
  # }
  # 
  # if (correl_fun=="RegPos"|correl_fun=="RegNeg") {
  #   #plot(A_fun,RNV_size[1,],col=2)
  #   int_var <- abs(1/B_fun)
  #   plot(int_var,RNV_size[1,],xlab="Absolute of 1/B",ylab="RNV size")
  #   #summary(lm(as.vector(RNV_size[1,])~A_fun))
  #   #summary(lm(as.vector(RNV_size[1,])~int_var))
  # }
  # 
  # # if (correl_fun=="RegNeg") {
  # #   #plot(A_fun,RNV_size[1,],col=2)
  # #   plot(abs(1/B_fun),RNV_size[1,],col=3)
  # #   int_var <- abs(1/B_fun)
  # #   summary(lm(as.vector(RNV_size[1,])~A_fun))
  # #   summary(lm(as.vector(RNV_size[1,])~int_var))
  # # }
  # 
  # if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
  #   #plot(A_fun,RNV_size[1,],col=2)
  #   plot(abs(1/B_fun),RNV_size[1,],xlab="Absolute of 1/B",ylab="RNV size")
  #   
  #   int_var <- abs(1/B_fun-e0_fun)
  #   plot(int_var,RNV_size[1,],xlab="Absolute of (1/B - e0)",ylab="RNV size")
  #   
  #   #z <- abs(1/B_fun)
  #   # summary(lm(as.vector(RNV_size[1,])~A_fun))
  #   # summary(lm(as.vector(RNV_size[1,])~abs(1/B_fun)))
  #   # summary(lm(as.vector(RNV_size[1,])~y))
  # }
  # 
  # 
  # lm_RNV <- lm(as.vector(RNV_size[1,])~int_var)
  # summary(lm_RNV)
  
  
  return(invisible(list(E_eq=E_eq,e_eq=e_eq,A=A_fun,B=B_fun,beta=beta_fun,E0=E0_fun,RNV_size=RNV_size,N=N_fun,n=n_fun,correl=correl_fun,lm_RNV=lm_RNV,rank_var=rank_var)))
}



