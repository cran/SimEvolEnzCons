#' Flux shape computing for all points
#'
#' \code{flux.shape.for.all.points} computes flux and selection coefficient for various points (i.ei. vector of concentrations), giving the flux shape
#' 
#' 
#'
#' @details 
#' Every enzyme correspond to one dimension in a \emph{n}-dimensional graph.
#' 
#' For various concentration vectors, this function computes flux and selection coefficient.
#' 
#' Selection coefficient are computed from two different expressions (discrete by \code{\link{coef_sel.discrete}} and continuous by \code{\link{coef_sel.continue}}), for a mutation of canonical size \code{nu_fun}.
#'  
#' \code{E_ini_fun} is rescaled by a cross product to have sum of \code{E_ini_fun} equal to \code{Etot_fun}.
#' 
#' 
#' 	
#' 
#' @usage flux.shape.for.all.points(Etot_fun, A_fun, nu_fun, correl_fun, beta_fun=NULL,
#'  E_ini_fun=NULL, n_fun=3)
#'     
#' @inheritParams flux.shape.from.one.point
#' @param nu_fun Numeric value of \bold{canonical} mutation effect
#' @param n_fun Numeric. Number of enzymes. Necessarily equal to 3 for this function.
#'
#'
#' @return Invisible list of 3 elements:
#' \itemize{
#'    \item \code{$J} : Numeric vector of flux.
#'    \item \code{$sel_disc} : Numeric matrix of \code{n_fun} columns corresponding to discrete selection coefficient, and each row corresponds to one point (i.e. concentration vector).
#'    Each column correspond to the "mutated" enzyme.
#'    \item \code{$sel_cont} :Same as \code{$sel_disc}, but for continuous selection coefficient.
#'    Same properties.
#'    }
#' 
#'
#'
#' @seealso 
#' To study shape of flux only from a certain point and for any number of enzymes, see \code{\link{flux.shape.from.one.point}}
#' 
#' 
#'   
#' @export







flux.shape.for.all.points <- function(Etot_fun, A_fun, nu_fun, correl_fun, beta_fun=NULL, E_ini_fun=NULL, n_fun=3) {
  
  ###### Setting
  B_fun <- compute.B.from.beta(beta_fun)
  req <- 0.05 #bounds around equilibrium
  
  ##### Test
  if (n_fun!=3) {
    stop("This function is available only for three enzymes : n_fun=3.")
  }
  if (length(A_fun)!=n_fun) {
    stop("An enzyme is missing. Length of A_fun need to be equal to n_fun.")
  }
  is.correl.authorized(correl_fun)
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  #Regulation : need E0
  if (length(E_ini_fun)==0&correl_fun!="SC"&correl_fun!="Comp") {
    stop("In case of regulation, initial point 'E_ini_fun' is needed.")
  }

  
  #### Initialization
  
  # recalibrate E0 to have sum(E0) = Etot_fun
  E_ini_fun <- E_ini_fun*Etot_fun/sum(E_ini_fun)
  
  # Choix des concentrations
  e1<-rep(seq(0,1, by=0.01), 100) #nb of pts in x
  e2<-rep(seq(0,1, by=0.01), each=100) #nb of pts in y
  e3<-rep(seq(0,1, by=0.01), each=10100) #nb of pts in z
  e_sup<-cbind(e1,e2,e3)
  
  if (correl_fun=="SC"){
    E_test <- Etot_fun*e_sup
    e_fun <- matrix(0,nrow=nrow(E_test),ncol=n_fun)
    for (j in 1:nrow(E_test)) {
      e_fun[j,] <- E_test[j,]/sum(E_test[j,])
    }
  }
  
  if (correl_fun=="Comp"){
    #choose only if sum(e)=1 and each e !=0
    e_tri <- e_sup[apply(e_sup,1,sum)==1&(e_sup[,1]!=0&e_sup[,2]!=0&e_sup[,3]!=0),]
    e_fun <- e_tri
    E_test <- Etot_fun*e_tri
  }
  
  if (correl_fun=="RegPos"|correl_fun=="RegNeg") {
    E_test <- NULL
    #test sequence
    Ei <- seq(0,n_fun*Etot_fun)
    #for each enzyme in test sequence
    for (i in 1:n_fun) {
      E_reg <- matrix(0,nrow=length(Ei),ncol=n_fun)
      E_reg[,i] <- Ei
      for (j in 1:n_fun) {
        #computes value for other enzymes
        if (j != i){
          E_reg[,j] <- E_ini_fun[j]+beta[i,j]*(Ei-E_ini_fun[i])
        }
      }
      # keep only positive concentrations
      E_exist <- E_reg[which(E_reg[,1]>=0&E_reg[,2]>=0&E_reg[,3]>=0),]
      E_test <- rbind(E_test,E_exist)
    }
    
    e_fun <- matrix(0,nrow=nrow(E_test),ncol=n_fun)
    for (j in 1:nrow(E_test)) {
      e_fun[j,] <- E_test[j,]/sum(E_test[j,])
    }
  }
  
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    #Etot_fun <- sum(E_ini_fun)
    tau_CR_ens <- seq(-1,2,by=0.01)
    E_test_ens <- matrix(0,nrow=length(tau_CR_ens),ncol=n_fun)
    for (j in 1:length(tau_CR_ens)) {
      E_test_ens[j,] <- droite_E.CR(tau_CR_ens[j],E_ini_fun,B_fun)
    }
    keep_E <- which(E_test_ens[,1]>=0&E_test_ens[,2]>=0&E_test_ens[,3]>=0)
    tau <- tau_CR_ens[keep_E]
    E_test <- E_test_ens[keep_E,]
    e_fun <- E_test/sum(Etot_fun)
    
    # #mÃ©thode  de calcul, pour vÃ©rif si blocage par 1/B
    # E_test <- NULL
    # alpha <- alpha_ij(E0,beta,B,correl_fun,n)
    # Ei <- seq(0,n*Etot) #sÃ©quence Ã  tester
    # for (i in 1:n) { #pour chaque enzyme prise dan sla sÃ©quence Ã  tester
    #   E_CR <- matrix(0,nrow=length(Ei),ncol=n)
    #   E_CR[,i] <- Ei
    #   for (j in 1:n) {
    #     if (j != i){ #on calcule les valeurs des autres enzymes
    #       E_CR[,j] <- E0[j]+alpha[i,j]*(Ei-E0[i])
    #     }
    #   }
    # 
    #   E_exist <- E_CR[which(E_CR[,1]>=0&E_CR[,2]>=0&E_CR[,3]>=0),]
    #   E_test <- rbind(E_test,E_exist)
    # }
    # e <- E_test/sum(E0)
    # tau <- apply(E_test,1,droite_tau,E0,B)[1,]
    
  }
  
  
  J_test <- apply(E_test,1,flux,A_fun)
  J_limtop <- 0.9*max(J_test)
  #J0 <- flux(E0,A)
  
  eq_th <- predict_th(A_fun,correl_fun,B_fun)
  #J_th <- flux(eq_th$pred_e*sum(E0),A)
  if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
    eq_eff <- predict_eff(E_ini_fun,B_fun,A_fun,correl_fun)
  }
  
  if (correl_fun=="RegPos"|correl_fun=="RegNeg") {
    tau <- t(apply(E_test,1,droite_tau,E_ini_fun,B_fun))
  }
  
  coefrep_test <- matrix(0,nrow=nrow(E_test),ncol=n_fun)
  for (j in 1:nrow(E_test)) {
    coefrep_test[j,] <- coef_rep(E_test[j,],A_fun,correl_fun,beta_fun)
  }
  
  coefsel_test <- matrix(0,nrow=nrow(E_test),ncol=n_fun)
  for (j in 1:nrow(E_test)) { #pour chaque triplet
    delta_fun  <- compute.delta(nu_fun,E_test[j,],correl_fun,B_fun) 
    #coefsel_test[j,] <- coef_sel(E_test[j,],n,A,beta,B,correl_fun,delta)
    coefsel_test[j,] <- coef_sel.continue(0,E_test[j,],A_fun,delta_fun,correl_fun,beta_fun)
  }
  
  coefsel_dis <- matrix(0,nrow=nrow(E_test),ncol=n_fun)
  for (i in 1:n_fun) {
    E_mut <- t(apply(E_test,1,mut.E.direct,i,nu_fun,correl_fun,beta_fun))
    Jmut <- apply(E_mut,1,flux,A_fun)
    coefsel_dis[,i] <- (Jmut-J_test)/J_test
    #coefsel_dis[,i] <- apply(E_test,1,coef_sel.discrete,A_fun,E_mut)
  }
  
  #for flux > 0.9*J_max
  e_top <- e_fun[which(J_test>J_limtop),]
  E_top <- E_test[which(J_test>J_limtop),]
  coefrep_top <-  coefsel_test[which(J_test>J_limtop),]
  coefsel_top <-  coefsel_test[which(J_test>J_limtop),]
  J_top <- J_test[which(J_test>J_limtop)]
  if (correl_fun=="RegPos") {
    tau_top <- tau[which(J_test>J_limtop)]
  }
  
  #around theoretical equilibrium
  if (correl_fun!="RegNeg"&correl_fun!="CRNeg"){
    select_req <- which(e_fun[,1]>=(eq_th$pred_e[1]-req)&e_fun[,1]<=(eq_th$pred_e[1]+req)&e_fun[,2]>=(eq_th$pred_e[2]-req)&e_fun[,2]<=(eq_th$pred_e[2]+req)&e_fun[,3]>=(eq_th$pred_e[3]-req)&e_fun[,3]<=(eq_th$pred_e[3]+req))
    E_req <- E_test[select_req,]
    e_req <- e_fun[select_req,]
    coefrep_req <- coefsel_test[select_req,]
    coefsel_req <- coefsel_test[select_req,]
    J_req <- J_test[select_req]
    if (correl_fun=="RegPos") {
      tau_req <- tau[select_req]
    }
  }
  
  #around theoretical equilibrium & flux > 0.9*J_max
  if (correl_fun!="RegNeg"&correl_fun!="CRNeg"){
    select_topreq <- which(e_top[,1]>=(eq_th$pred_e[1]-req)&e_top[,1]<=(eq_th$pred_e[1]+req)&e_top[,2]>=(eq_th$pred_e[2]-req)&e_top[,2]<=(eq_th$pred_e[2]+req)&e_top[,3]>=(eq_th$pred_e[3]-req)&e_top[,3]<=(eq_th$pred_e[3]+req))
    E_topreq <- E_top[select_topreq,]
    e_topreq <- e_top[select_topreq,]
    coefrep_topreq <- coefsel_top[select_topreq,]
    coefsel_topreq <- coefsel_top[select_topreq,]
    J_topreq <- J_top[select_topreq]
  }
  
  #around effective equilibrium
  if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg"){
    select_reff <- which(e_fun[,1]>=(eq_eff$pred_e[1]-req)&e_fun[,1]<=(eq_eff$pred_e[1]+req)&e_fun[,2]>=(eq_eff$pred_e[2]-req)&e_fun[,2]<=(eq_eff$pred_e[2]+req)&e_fun[,3]>=(eq_eff$pred_e[3]-req)&e_fun[,3]<=(eq_eff$pred_e[3]+req))
    E_reff <- E_test[select_reff,]
    e_reff <- e_fun[select_reff,]
    coefrep_reff <- coefsel_test[select_reff,]
    coefsel_reff <- coefsel_test[select_reff,]
    J_reff <- J_test[select_reff]
    tau_reff <- tau[select_reff]
  }
  
  return(invisible(list(J=J_test,sel_cont=coefsel_test,sel_dis=coefsel_dis)))
}







###########" Graphics

# @rdname flux.shape.in.all.dimension
# @usage 
# @aliases flux.shape.in.all.dimension.graphics
# @return \code{flux.shape.in.all.dimension.graphics} returns a PDF document with graphics of \code{flux.shape.from.one.point} results
# @export







