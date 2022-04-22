#' Bayesian Illness-Death Model with Weibull Baselines, Semi-Markov structure and logit model for immediate event
#'
#' @param Formula formula
#' @param data data
#' @param na.action how NAs are treated. See \code{\link[stats]{model.frame}}.
#' @param subset a specification of the rows to be used: defaults to all rows. See \code{\link[stats]{model.frame}}.
#' @param hyperParams list of hyperparameters
#' @param mcmcParams list of hyperparameters
#' @param n_chains integer for number of chains to fit
#' @param start_mat a matrix of start values, with each column corresponding to a chain. If null, n_chains sets number of chains
#'
#' @return a list with outputs
#' @import Formula
#' @export
Bayes_SCR_logit <- function(Formula, data, na.action="na.fail", subset=NULL,
                      hyperParams, mcmcParams, n_chains, start_mat=NULL){
  browser()

  ##INITIALIZE DATA##
  ##*******************************##
  #rearrange input Formula object (which stores the different pieces of the input formula)
  #This line ensures that the formula is of type Formula, and not just formula
  #the below manipulations require special methods of Formula.
  form2 <- Formula::as.Formula(paste0(Formula[2], Formula[1], Formula[3]))
  #arrange data into correct form, extracting subset of variables
  #and observations needed in model
  data <- stats::model.frame(form2,data=data,na.action=na.action,subset=subset)
  #create matrices storing two outcomes, and then component vectors
  time1 <- Formula::model.part(form2, data=data, lhs=1)
  time2 <- Formula::model.part(form2, data=data, lhs=2)

  ### define flags characterizing three key groups
  ##
  #who has nonterminal event at all (aka at risk for immediate death)
  delta1 <- time1[[2]]
  #who has nonterminal event and immediate death
  delta1D <- Formula::model.part(form2, data=data, lhs=3)[[1]]
  #who has nonterminal event without immediate death (aka at risk for h3)
  delta1noD <- as.numeric(delta1==1 & delta1D==0)

#
#   #make vector of indices corresponding with those who experience nonterminal event
#   delta1_index = which(delta1==1)
#   #make n-length lookup vector connecting the index out of n with the index out of sum(delta1)
#   delta1_index_lookup <- delta1
#   delta1_index_lookup[delta1==1] <- delta1_index
#
#   #make vector of indices corresponding with those who experience nonterminal event without immediate death
#   delta1noD_index <- which(delta1noD==1)
#   delta1noD_index_lookup <- delta1noD
#   delta1noD_index_lookup[delta1noD==1] <- delta1noD_index

  #immediate death indicator among those who have non-terminal event
  delta1D_sub <- delta1D[delta1==1]

  #time of first event for everyone
  y1 <- time1[[1]]
  #indicator for whether first event was observed terminal event ("cr" = competing risk)
  delta_cr <- ifelse(time1[[1]] < time2[[1]], 0, time2[[2]])

  #time of death/censoring among those who have non-terminal event without immediate death
  y_sm <- (time2[[1]] - time1[[1]])[delta1noD==1]
  #death indicator among those who have non-terminal event without immediate death
  delta_sm <- time2[[2]][delta1noD==1]

  #ensure that number of obs in h3 is correctly accounted for
  stopifnot(length(y_sm) == sum(delta1noD))


  #Create covariate matrices for each of three transition hazards
  Xmat1 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=1),
                                        data=data))
  Xmat2 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=2),
                                        data=data))
  #subsetted to just those who experience the non-terminal event and don't immediately experience terminal event
  Xmat3 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=3),
                                        data=data))[delta1noD==1,,drop=FALSE]
  #subsetted to just those who experience the non-terminal event, INTERCEPT ADDED
  XmatD <- cbind(rep(1,sum(delta1)),
                 as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=4),
                                        data=data))[delta1==1,,drop=FALSE])
  p1 <- ncol(Xmat1)
  p2 <- ncol(Xmat2)
  p3 <- ncol(Xmat3)
  pD <- ncol(XmatD) #includes intercept
  n <- length(y1)

  # survreg(formula = Surv(y1,delta1) ~ Xmat1)
  # survreg(formula = Surv(y1,delta_cr) ~ Xmat2)
  # survreg(formula = Surv(y_sm,delta_sm) ~ Xmat3)


  ####SET HYPERPARAMETERS####
  hyper_vec <- as.vector(c(hyperParams$WB$WB.ab1, hyperParams$WB$WB.ab2, hyperParams$WB$WB.ab3,
                           hyperParams$WB$WB.cd1, hyperParams$WB$WB.cd2, hyperParams$WB$WB.cd3,
                           hyperParams$theta))
  tuning_vec <- as.vector(c(mhProp_alpha1_var=mcmcParams$tuning$mhProp_alphag_var[1],
                            mhProp_alpha2_var=mcmcParams$tuning$mhProp_alphag_var[2],
                            mhProp_alpha3_var=mcmcParams$tuning$mhProp_alphag_var[3],
                            mhProp_theta_var=mcmcParams$tuning$mhProp_theta_var))

  ####SET MCMC SETTINGS####
  numReps=mcmcParams$run$numReps
  thin=mcmcParams$run$thin
  burninPerc=mcmcParams$run$burninPerc

  n_burnin <- numReps * burninPerc
  n_sample <- numReps - n_burnin
  n_store  <- numReps/thin * (1 - burninPerc)
  stopifnot(n_burnin %% 1 == 0)
  stopifnot(n_sample > 0)
  if(n_store %% 1 != 0){ stop("numReps * burninPerc  must be divisible by thin")}


  ####ASSIGN START VALUES####

  if(is.null(start_mat)){
    if(is.null(n_chains)){stop("either supply non-null start_list, or specify n_chains.")}
    theta_start_temp <- stats::runif(n_chains, 0.1, 1.1) #theta
    start_mat <- rbind(
      stats::runif(n_chains,0.1,1.1), #kappa1
      stats::runif(n_chains,0.8,1.2), #alpha1
      stats::runif(n_chains,0.1,1.1), #kappa2
      stats::runif(n_chains,0.8,1.2), #alpha2
      stats::runif(n_chains,0.1,1.1), #kappa3
      stats::runif(n_chains,0.8,1.2), #alpha3
      theta_start_temp, #theta
      if(p1>0) matrix(stats::runif(n_chains*p1, -0.1, 0.1),ncol=n_chains), #beta1
      if(p2>0) matrix(stats::runif(n_chains*p2, -0.1, 0.1),ncol=n_chains), #beta2
      if(p3>0) matrix(stats::runif(n_chains*p3, -0.1, 0.1),ncol=n_chains), #beta3
      if(pD>0) matrix(stats::runif(n_chains*pD, -0.1, 0.1),ncol=n_chains), #betaD
      matrix(stats::rgamma(n*n_chains, 1/theta_start_temp, 1/theta_start_temp),
             ncol = n_chains,byrow=TRUE) #frailties
    )
    # theta_start_temp <- stats::runif(1, 0.1, 1.1)
    # start_vec <- c(kappa1=0.1,alpha1=1, kappa2=0.1, alpha2=1,
    #                0.1,1, theta_start_temp,
    #                if(p1>0) stats::runif(p1, -0.1, 0.1),
    #                if(p2>0) stats::runif(p2, -0.1, 0.1),
    #                if(p3>0) stats::runif(p3, -0.1, 0.1),
    #                stats::rgamma(n, 1/theta_start_temp, 1/theta_start_temp))
    # start_mat <- matrix(data=start_vec,ncol=n_chains,nrow=length(start_vec))
    rownames(start_mat) <- c("kappa1","alpha1","kappa2","alpha2","kappa3","alpha3",
                             "theta",
                             if(p1>0) paste0("beta1_",1:p1),
                             if(p2>0) paste0("beta2_",1:p2),
                             if(p3>0) paste0("beta3_",1:p3),
                             if(pD>0) paste0("betaD_",1:pD),
                             paste0("gamma",1:n))
  }

  # #TODO: parallelize this loop (I know it can be done!)
  out_list <- list()
  # #generate an array to store the resulting samples
  out_list[["samples"]] <- array(dim = c(n_store, n_chains, nrow(start_mat)),
                                 dimnames = list(as.character(1:n_store),
                                                 paste0("chain:",1:n_chains),
                                                 rownames(start_mat)))
  # out_list[["accept"]] <- list()
  mcmcRet <- list()
  for(i in 1:n_chains){
    print(paste0("Chain: ", i))
    mcmcRet <- WeibSCRlogitmcmc(y1, y_sm,
                                delta1, #indicator for nonterminal event followed by immediate death
                                delta1noD, #indicator for nonterminal event not followed by immediate death
                                delta_cr, #terminal first event indicator
                                delta_sm, #terminal event indicator in group with nonterminal first event and no immediate terminal
                                delta1D_sub, #immediate terminal event indicator in group with nonterminal first event
                                Xmat1, Xmat2, Xmat3, XmatD,
                                hyper_vec, tuning_vec, start_mat[,i],
                                n_burnin, n_sample, thin)

    out_list[["samples"]][,i,1] <- mcmcRet[["samples"]][["kappa1"]]
    out_list[["samples"]][,i,2] <- mcmcRet[["samples"]][["alpha1"]]
    out_list[["samples"]][,i,3] <- mcmcRet[["samples"]][["kappa2"]]
    out_list[["samples"]][,i,4] <- mcmcRet[["samples"]][["alpha2"]]
    out_list[["samples"]][,i,5] <- mcmcRet[["samples"]][["kappa3"]]
    out_list[["samples"]][,i,6] <- mcmcRet[["samples"]][["alpha3"]]
    out_list[["samples"]][,i,7] <- mcmcRet[["samples"]][["theta"]]
    if(p1>0) out_list[["samples"]][,i,8:(7+p1)] <- mcmcRet[["samples"]][["beta1"]]
    if(p2>0) out_list[["samples"]][,i,(8+p1):(7+p1+p2)] <- mcmcRet[["samples"]][["beta2"]]
    if(p3>0) out_list[["samples"]][,i,(8+p1+p2):(7+p1+p2+p3)] <- mcmcRet[["samples"]][["beta3"]]
    if(pD>0) out_list[["samples"]][,i,(8+p1+p2+p3):(7+p1+p2+p3+pD)] <- mcmcRet[["samples"]][["betaD"]]
    out_list[["samples"]][,i,(8+p1+p2+p3+pD):(7+p1+p2+p3+pD+n)] <- mcmcRet[["samples"]][["gamma"]]

    out_list[["accept"]][[paste0("chain",i)]] <- mcmcRet$accept

  }

  #for now, my plan is going to be to leverage the bayesplot package to
  #make visuals

  out_list[["setup"]]	<-
    # mcmcRet[["setup"]]	<-
    list(start_mat = start_mat,
         mcmc_para = c(n_burnin=n_burnin,n_sample=n_sample,thin=thin))

  # return(mcmcRet)
  return(out_list)

}




#this is simple experimental version with a covariate in the logit model.



#' Bayesian Illness-Death Model with Weibull Baselines, Semi-Markov structure and logit model for immediate event
#'
#' @param Formula formula
#' @param data data
#' @param na.action how NAs are treated. See \code{\link[stats]{model.frame}}.
#' @param subset a specification of the rows to be used: defaults to all rows. See \code{\link[stats]{model.frame}}.
#' @param hyperParams list of hyperparameters
#' @param mcmcParams list of hyperparameters
#' @param n_chains integer for number of chains to fit
#' @param start_mat a matrix of start values, with each column corresponding to a chain. If null, n_chains sets number of chains
#'
#' @return a list with outputs
#' @import Formula
#' @export
Bayes_SCR_logit2 <- function(Formula, data, na.action="na.fail", subset=NULL,
                            hyperParams, mcmcParams, n_chains, start_mat=NULL){
  browser()

  ##INITIALIZE DATA##
  ##*******************************##
  #rearrange input Formula object (which stores the different pieces of the input formula)
  #This line ensures that the formula is of type Formula, and not just formula
  #the below manipulations require special methods of Formula.
  form2 <- Formula::as.Formula(paste0(Formula[2], Formula[1], Formula[3]))
  #arrange data into correct form, extracting subset of variables
  #and observations needed in model
  data <- stats::model.frame(form2,data=data,na.action=na.action,subset=subset)
  #create matrices storing two outcomes, and then component vectors
  time1 <- Formula::model.part(form2, data=data, lhs=1)
  time2 <- Formula::model.part(form2, data=data, lhs=2)

  ### define flags characterizing three key groups
  ##
  #who has nonterminal event at all (aka at risk for immediate death)
  delta1 <- time1[[2]]
  #who has nonterminal event and immediate death
  delta1D <- Formula::model.part(form2, data=data, lhs=3)[[1]]
  #who has nonterminal event without immediate death (aka at risk for h3)
  delta1noD <- as.numeric(delta1==1 & delta1D==0)

  #
  #   #make vector of indices corresponding with those who experience nonterminal event
  #   delta1_index = which(delta1==1)
  #   #make n-length lookup vector connecting the index out of n with the index out of sum(delta1)
  #   delta1_index_lookup <- delta1
  #   delta1_index_lookup[delta1==1] <- delta1_index
  #
  #   #make vector of indices corresponding with those who experience nonterminal event without immediate death
  #   delta1noD_index <- which(delta1noD==1)
  #   delta1noD_index_lookup <- delta1noD
  #   delta1noD_index_lookup[delta1noD==1] <- delta1noD_index

  #immediate death indicator among those who have non-terminal event
  delta1D_sub <- delta1D[delta1==1]

  #time of first event for everyone
  y1 <- time1[[1]]
  #indicator for whether first event was observed terminal event ("cr" = competing risk)
  delta_cr <- ifelse(time1[[1]] < time2[[1]], 0, time2[[2]])

  #time of death/censoring among those who have non-terminal event without immediate death
  y_sm <- (time2[[1]] - time1[[1]])[delta1noD==1]
  #death indicator among those who have non-terminal event without immediate death
  delta_sm <- time2[[2]][delta1noD==1]

  #ensure that number of obs in h3 is correctly accounted for
  stopifnot(length(y_sm) == sum(delta1noD))


  #Create covariate matrices for each of three transition hazards
  Xmat1 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=1),
                                        data=data))
  Xmat2 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=2),
                                        data=data))
  #subsetted to just those who experience the non-terminal event and don't immediately experience terminal event
  Xmat3 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=3),
                                        data=data))[delta1noD==1,,drop=FALSE]
  #subsetted to just those who experience the non-terminal event, INTERCEPT ADDED
  XmatD <- cbind(numeric(sum(delta1)), #placeholder column for frailties
                 rep(1,sum(delta1)),
                 as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=4),
                                              data=data))[delta1==1,,drop=FALSE])
  p1 <- ncol(Xmat1)
  p2 <- ncol(Xmat2)
  p3 <- ncol(Xmat3)
  pD <- ncol(XmatD) #includes intercept and column of frailties
  n <- length(y1)

  # survreg(formula = Surv(y1,delta1) ~ Xmat1)
  # survreg(formula = Surv(y1,delta_cr) ~ Xmat2)
  # survreg(formula = Surv(y_sm,delta_sm) ~ Xmat3)


  ####SET HYPERPARAMETERS####
  hyper_vec <- as.vector(c(hyperParams$WB$WB.ab1, hyperParams$WB$WB.ab2, hyperParams$WB$WB.ab3,
                           hyperParams$WB$WB.cd1, hyperParams$WB$WB.cd2, hyperParams$WB$WB.cd3,
                           hyperParams$theta))
  tuning_vec <- as.vector(c(mhProp_alpha1_var=mcmcParams$tuning$mhProp_alphag_var[1],
                            mhProp_alpha2_var=mcmcParams$tuning$mhProp_alphag_var[2],
                            mhProp_alpha3_var=mcmcParams$tuning$mhProp_alphag_var[3],
                            mhProp_theta_var=mcmcParams$tuning$mhProp_theta_var))

  ####SET MCMC SETTINGS####
  numReps=mcmcParams$run$numReps
  thin=mcmcParams$run$thin
  burninPerc=mcmcParams$run$burninPerc

  n_burnin <- numReps * burninPerc
  n_sample <- numReps - n_burnin
  n_store  <- numReps/thin * (1 - burninPerc)
  stopifnot(n_burnin %% 1 == 0)
  stopifnot(n_sample > 0)
  if(n_store %% 1 != 0){ stop("numReps * burninPerc  must be divisible by thin")}


  ####ASSIGN START VALUES####

  if(is.null(start_mat)){
    if(is.null(n_chains)){stop("either supply non-null start_list, or specify n_chains.")}
    theta_start_temp <- stats::runif(n_chains, 0.1, 1.1) #theta
    start_mat <- rbind(
      stats::runif(n_chains,0.1,1.1), #kappa1
      stats::runif(n_chains,0.8,1.2), #alpha1
      stats::runif(n_chains,0.1,1.1), #kappa2
      stats::runif(n_chains,0.8,1.2), #alpha2
      stats::runif(n_chains,0.1,1.1), #kappa3
      stats::runif(n_chains,0.8,1.2), #alpha3
      theta_start_temp, #theta
      if(p1>0) matrix(stats::runif(n_chains*p1, -0.1, 0.1),ncol=n_chains), #beta1
      if(p2>0) matrix(stats::runif(n_chains*p2, -0.1, 0.1),ncol=n_chains), #beta2
      if(p3>0) matrix(stats::runif(n_chains*p3, -0.1, 0.1),ncol=n_chains), #beta3
      if(pD>0) matrix(stats::runif(n_chains*pD, -0.1, 0.1),ncol=n_chains), #betaD
      matrix(stats::rgamma(n*n_chains, 1/theta_start_temp, 1/theta_start_temp),
             ncol = n_chains,byrow=TRUE) #frailties
    )
    # theta_start_temp <- stats::runif(1, 0.1, 1.1)
    # start_vec <- c(kappa1=0.1,alpha1=1, kappa2=0.1, alpha2=1,
    #                0.1,1, theta_start_temp,
    #                if(p1>0) stats::runif(p1, -0.1, 0.1),
    #                if(p2>0) stats::runif(p2, -0.1, 0.1),
    #                if(p3>0) stats::runif(p3, -0.1, 0.1),
    #                stats::rgamma(n, 1/theta_start_temp, 1/theta_start_temp))
    # start_mat <- matrix(data=start_vec,ncol=n_chains,nrow=length(start_vec))
    rownames(start_mat) <- c("kappa1","alpha1","kappa2","alpha2","kappa3","alpha3",
                             "theta",
                             if(p1>0) paste0("beta1_",1:p1),
                             if(p2>0) paste0("beta2_",1:p2),
                             if(p3>0) paste0("beta3_",1:p3),
                             if(pD>0) paste0("betaD_",1:pD),
                             paste0("gamma",1:n))
  }

  # #TODO: parallelize this loop (I know it can be done!)
  out_list <- list()
  # #generate an array to store the resulting samples
  out_list[["samples"]] <- array(dim = c(n_store, n_chains, nrow(start_mat)),
                                 dimnames = list(as.character(1:n_store),
                                                 paste0("chain:",1:n_chains),
                                                 rownames(start_mat)))
  # out_list[["accept"]] <- list()
  mcmcRet <- list()
  for(i in 1:n_chains){
    print(paste0("Chain: ", i))
    mcmcRet <- WeibSCRlogitmcmc2(y1, y_sm,
                                delta1, #indicator for nonterminal event followed by immediate death
                                delta1noD, #indicator for nonterminal event not followed by immediate death
                                delta_cr, #terminal first event indicator
                                delta_sm, #terminal event indicator in group with nonterminal first event and no immediate terminal
                                delta1D_sub, #immediate terminal event indicator in group with nonterminal first event
                                Xmat1, Xmat2, Xmat3, XmatD,
                                hyper_vec, tuning_vec, start_mat[,i],
                                n_burnin, n_sample, thin)

    out_list[["samples"]][,i,1] <- mcmcRet[["samples"]][["kappa1"]]
    out_list[["samples"]][,i,2] <- mcmcRet[["samples"]][["alpha1"]]
    out_list[["samples"]][,i,3] <- mcmcRet[["samples"]][["kappa2"]]
    out_list[["samples"]][,i,4] <- mcmcRet[["samples"]][["alpha2"]]
    out_list[["samples"]][,i,5] <- mcmcRet[["samples"]][["kappa3"]]
    out_list[["samples"]][,i,6] <- mcmcRet[["samples"]][["alpha3"]]
    out_list[["samples"]][,i,7] <- mcmcRet[["samples"]][["theta"]]
    if(p1>0) out_list[["samples"]][,i,8:(7+p1)] <- mcmcRet[["samples"]][["beta1"]]
    if(p2>0) out_list[["samples"]][,i,(8+p1):(7+p1+p2)] <- mcmcRet[["samples"]][["beta2"]]
    if(p3>0) out_list[["samples"]][,i,(8+p1+p2):(7+p1+p2+p3)] <- mcmcRet[["samples"]][["beta3"]]
    if(pD>0) out_list[["samples"]][,i,(8+p1+p2+p3):(7+p1+p2+p3+pD)] <- mcmcRet[["samples"]][["betaD"]]
    out_list[["samples"]][,i,(8+p1+p2+p3+pD):(7+p1+p2+p3+pD+n)] <- mcmcRet[["samples"]][["gamma"]]

    out_list[["accept"]][[paste0("chain",i)]] <- mcmcRet$accept

  }

  #for now, my plan is going to be to leverage the bayesplot package to
  #make visuals

  out_list[["setup"]]	<-
    # mcmcRet[["setup"]]	<-
    list(start_mat = start_mat,
         mcmc_para = c(n_burnin=n_burnin,n_sample=n_sample,thin=thin))

  # return(mcmcRet)
  return(out_list)

}
