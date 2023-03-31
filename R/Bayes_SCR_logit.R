#version where there is an explicit window of width epsilon, and a landmark-style hazard after
#to simplify things, I'm just going to directly accept inputs for times:
  #y1 (first event time),
  #y_sm (sojourn time post landmark of y1+eps),
#and relevant indicators:
  #delta1 (non-terminal event observed or not)
  #delta2 (terminal event observed or not)
  #deltaD (terminal event observed within "immediate)

#in fact, I think doing it this way will kind of subsume the earlier point mass approach
#upon setting epsilon = 0.

#' Bayesian Illness-Death Model with Weibull Baselines, Semi-Markov structure and logit model for immediate event
#'
#' @param y1 time to first event (nonterminal, terminal only, censored)
#' @param y_sm sojourn time from end of terminal window (y1+eps) (not length n! Just as long as sum(delta1D))
#' @param delta1 indicator for nonterminal event
#' @param delta1D indicator for nonterminal event followed by immediate death
#' @param delta2 indicator for terminal event
#' @param eps width of "immediate" event window
#' @param data data
#' @param vars1,vars2,vars3,varsD string vectors with variable names for each submodel
#' @param na.action how NAs are treated. See \code{\link[stats]{model.frame}}.
#' @param subset a specification of the rows to be used: defaults to all rows. See \code{\link[stats]{model.frame}}.
#' @param hyperParams list of hyperparameters
#' @param mcmcParams list of mcmc control parameters
#' @param n_chains integer for number of chains to fit
#' @param start_mat a matrix of start values, with each column corresponding to a chain. If null, n_chains sets number of chains
#' @param frailty Boolean indicating whether a gamma distributed subject-specific frailty should
#'   be included.
#' @param frail_path filepath for where to write frailty output
#' @param logLHi_path filepath for where to write individual likelihood contribution output
#'
#'
#' @return a list with outputs
#' @import Formula
#' @export
Bayes_SCR_logit_eps <- function(y1, y_sm, delta1, delta1D, delta2, eps,
                            vars1, vars2, vars3, varsD,
                            data, na.action="na.fail", subset=NULL,
                            hyperParams, mcmcParams, n_chains, start_mat=NULL, frailty=TRUE,
                            frail_path=NULL, logLHi_path=NULL){
  # browser()

  ##INITIALIZE DATA##
  ##*******************************##

  #provided values
  y1=y1 #time to first event (nonterminal, terminal only, censored)
  y_sm=y_sm #sojourn time from end of terminal window (y1+eps) [not length n! Just as long as sum(delta1D)]
  delta1=delta1 #indicator for nonterminal event
  delta1D=delta1D #indicator for nonterminal event followed by immediate death
  delta2=delta2 #indicator for terminal event

  #derived values of length n
  delta1noD <- as.numeric(delta1==1 & delta1D==0) #indicator for nonterminal event not followed by immediate death
  delta_cr= as.numeric(delta1==0 & delta2==1) #terminal first event indicator

  #derived values on subsets
  delta1D_sub <- delta1D[delta1==1] #immediate terminal event indicator in group with nonterminal first event
  delta_sm <- delta2[delta1noD==1] #terminal event indicator in group with nonterminal first event and no immediate terminal

  #note this is not necessarily the all time max!! just the max "first event"
  ymax <- max(y1)

  y_sm_sub <- y_sm[delta1noD==1]
  # #ensure that number of obs in h3 is correctly accounted for
  # stopifnot(length(y_sm) == sum(delta1noD))

  #Create covariate matrices for each of three transition hazards

  Xmat1 <- if(length(vars1)>0) as.matrix(data[,vars1,drop=FALSE]) else matrix(nrow = n,ncol=0)
  Xmat2 <- if(length(vars2)>0) as.matrix(data[,vars2,drop=FALSE]) else matrix(nrow = n,ncol=0)
  #subsetted to just those who experience the non-terminal event and don't immediately experience terminal event
  Xmat3 <- if(length(vars3)>0) as.matrix(data[delta1noD==1,vars3,drop=FALSE]) else matrix(nrow = sum(delta1noD==1),ncol=0)
  #subsetted to just those who experience the non-terminal event, INTERCEPT ADDED
  if(length(varsD)>0){
    XmatD <- cbind(logit_intercept=rep(1,sum(delta1)),
                   as.matrix(data[delta1==1,varsD,drop=FALSE]))
  } else {
    XmatD <- cbind(logit_intercept=rep(1,sum(delta1)),
                   matrix(nrow = sum(delta1noD==1),ncol=0))
  }
  p1 <- ncol(Xmat1)
  p2 <- ncol(Xmat2)
  p3 <- ncol(Xmat3)
  pD <- ncol(XmatD) #includes intercept
  n <- length(y1)

  ####SET HYPERPARAMETERS####
  hyper_vec <- as.vector(c(hyperParams$WB$WB.ab1, hyperParams$WB$WB.ab2, hyperParams$WB$WB.ab3,
                           hyperParams$WB$WB.cd1, hyperParams$WB$WB.cd2, hyperParams$WB$WB.cd3,
                           if(frailty) hyperParams$theta else c(0.7,0.7), #placeholders
                           if(is.null(hyperParams$logit$logit.m0)) numeric(pD) else hyperParams$logit$logit.m0,
                           if(is.null(hyperParams$logit$logit.P0diag)) numeric(pD) else hyperParams$logit$logit.P0diag
  ))
  tuning_vec <- as.vector(c(mhProp_alpha1_var=mcmcParams$tuning$mhProp_alphag_var[1],
                            mhProp_alpha2_var=mcmcParams$tuning$mhProp_alphag_var[2],
                            mhProp_alpha3_var=mcmcParams$tuning$mhProp_alphag_var[3],
                            mhProp_theta_var= if(frailty) mcmcParams$tuning$mhProp_theta_var else 0.1)) #placeholder

  ####SET MCMC SETTINGS####
  numReps=mcmcParams$run$numReps
  thin=mcmcParams$run$thin
  burninPerc=mcmcParams$run$burninPerc

  #set up where if it's less than 0, set to n, otherwise take what it is.
  nGam_save=if(mcmcParams$storage$nGam_save <0) n else mcmcParams$storage$nGam_save
  nlogLHi_save=if(mcmcParams$storage$nlogLHi_save <0) n else mcmcParams$storage$nlogLHi_save


  n_burnin <- numReps * burninPerc
  n_sample <- numReps - n_burnin
  n_store  <- numReps/thin * (1 - burninPerc)
  stopifnot(n_burnin %% 1 == 0)
  stopifnot(n_sample > 0)
  if(n_store %% 1 != 0){ stop("numReps * burninPerc  must be divisible by thin")}

  if(!is.null(frail_path)){
    dir.create(paste(frail_path), recursive = TRUE, showWarnings = FALSE)
  }
  if(!is.null(logLHi_path)){
    dir.create(paste(logLHi_path), recursive = TRUE, showWarnings = FALSE)
  }

  ####ASSIGN START VALUES####
  if(is.null(start_mat)){
    if(is.null(n_chains)){stop("either supply non-null start_list, or specify n_chains.")}
    theta_start_temp <- if(frailty) stats::runif(n_chains, 0.1, 1.1) else numeric(n_chains) #theta
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
      #if frailties included, draw them randomly, otherwise set them to 1.
      if(frailty) matrix(stats::rgamma(n*n_chains, 1/theta_start_temp, 1/theta_start_temp),ncol = n_chains,byrow=TRUE) else matrix(1,nrow=n,ncol=n_chains) #frailties (initialized by gamma, even if drawn from log-normal)
    )
    rownames(start_mat) <- c("kappa1","alpha1","kappa2","alpha2","kappa3","alpha3",
                             "theta", if(p1>0) paste0("beta1_",1:p1),
                             if(p2>0) paste0("beta2_",1:p2),
                             if(p3>0) paste0("beta3_",1:p3),
                             if(pD>0) paste0("betaD_",1:pD),
                             paste0("gamma",1:n))
  }

  gh_weights <- get_ghquad_pointsweights(n_quad = 15)$weights
  gh_nodes <- get_ghquad_pointsweights(n_quad = 15)$points

  #### PREALLOCATE OUTPUT LIST ####
  out_list <- list(
    #generate an array to store the resulting samples
    samples = array(dim = c(n_store, n_chains, 7 + p1 + p2 + p3 + pD), #even if frailty==FALSE, leave empty row for it.
                    dimnames = list(as.character(1:n_store),
                                    paste0("chain:",1:n_chains),
                                    rownames(start_mat)[1:(7 + p1 + p2 + p3 + pD)])),
    covnames = list(covNames1= if(p1 > 0) colnames(Xmat1) else NULL,
                    covNames2= if(p2 > 0) colnames(Xmat2) else NULL,
                    covNames3= if(p3 > 0) colnames(Xmat3) else NULL,
                    covNamesD= if(pD > 1) colnames(XmatD[,-1,drop=FALSE]) else NULL),
    #list final characteristics useful to export
    setup = list(Formula=Formula,
                 frailty=frailty,
                 nCov0 = c(2,2,2),
                 nCov = c(p1,p2,p3,pD),
                 hyper_vec = hyper_vec, start_mat = start_mat,
                 tuning_vec = tuning_vec,
                 nGam_save = nGam_save, nlogLHi_save = nlogLHi_save,
                 numReps = numReps, thin = thin,
                 burninPerc = burninPerc, n_store=n_store,
                 ymax=ymax, n=n,
                 frail_path = frail_path, logLHi_path = logLHi_path,
                 hz.type = "Weibull", model = "semi-Markov", nChain = n_chains,
                 mcmc_para = c(n_burnin=n_burnin,n_sample=n_sample,thin=thin)), #slight duplication but that's ok
    #empty list to keep "accept" values in.
    accept = vector(mode = "list", length = n_chains),
    move = vector(mode = "list", length = n_chains),
    diagnostics=list(dev=NA, DIC = NA, LPML = NA, dev_marg=NA, DIC_marg = NA, LPML_marg = NA,
                     #the mean log likelihood at each stored iteration
                     logLH_mat = matrix(data = NA,nrow = n_store, ncol = n_chains),
                     logLH_marg_mat = if(frailty) matrix(data = NA,nrow = n_store, ncol = n_chains) else NULL,
                     #the running mean of the likelihood and inverse likelihood for each subject across the samples
                     LH_mean_mat = matrix(data = NA,nrow = n, ncol = n_chains),
                     invLH_mean_mat = matrix(data = NA,nrow = n, ncol = n_chains),
                     LH_marg_mean_mat = if(frailty) matrix(data = NA,nrow = n, ncol = n_chains) else NULL,
                     invLH_marg_mean_mat = if(frailty) matrix(data = NA,nrow = n, ncol = n_chains) else NULL),
    class = c("Bayes_HReg2", "IDlogit", "Ind", "WB")
  )
  names(out_list$accept) <- paste0("chain",1:n_chains)

  #### RUN SAMPLER ####
  # #TODO: parallelize this loop (I know it can be done!)
  for(i in 1:n_chains){
    print(paste0("Chain: ", i))

    #initialize all of the vectors of individually sampled parameters
    sample_alpha1 <- numeric(n_store)
    sample_alpha2 <- numeric(n_store)
    sample_alpha3 <- numeric(n_store)
    sample_kappa1 <- numeric(n_store)
    sample_kappa2 <- numeric(n_store)
    sample_kappa3 <- numeric(n_store)
    sample_theta <- numeric(n_store)

    #initialize all of the matrices of regression parameters
    sample_beta1 <- matrix(data=0,nrow=p1,ncol=n_store)
    sample_beta2 <- matrix(data=0,nrow=p2,ncol=n_store)
    sample_beta3 <- matrix(data=0,nrow=p3,ncol=n_store)
    sample_betaD <- matrix(data=0,nrow=pD,ncol=n_store)

    #initialize matrices to store frailties and log-likelihood, as needed
    sample_frail <- matrix(data=0,nrow=nGam_save,ncol=n_store)
    sample_logLHi <- matrix(data=0,nrow=nlogLHi_save,ncol=n_store)
    sample_logLHi_marg <- if(frailty) matrix(data=0,nrow=nlogLHi_save,ncol=n_store) else matrix(nrow=0,ncol=0)

    #initialize integers counting acceptance of individually sampled parameters
    accept_base <- numeric(8)

    #initialize vectors counting acceptance of sampled vectors
    accept_frail <- numeric(n)
    accept_beta1 <- numeric(p1)
    accept_beta2 <- numeric(p2)
    accept_beta3 <- numeric(p3)

    #initialize vectors for diagnostics
    LH_marg_mean_vec <- if(frailty) numeric(n) else numeric(0)
    invLH_marg_mean_vec <- if(frailty) numeric(n) else numeric(0)
    sample_logLH_marg <- if(frailty) numeric(n_store) else numeric(0)

    LH_mean_vec <- numeric(n)
    invLH_mean_vec <- numeric(n)
    sample_logLH <- numeric(n_store)

    #initialize move vector
    move_vec <- numeric(numReps)

    mcmcRet <- WeibSCRlogitmcmc(y1=y1, y_sm=y_sm_sub,
                                delta1=delta1, #indicator for nonterminal event followed by immediate death
                                delta1noD=delta1noD, #indicator for nonterminal event not followed by immediate death
                                delta_cr=delta_cr, #terminal first event indicator
                                delta_sm=delta_sm, #terminal event indicator in group with nonterminal first event and no immediate terminal
                                delta1D_sub=delta1D_sub, #immediate terminal event indicator in group with nonterminal first event
                                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, XmatD=XmatD, eps=eps,
                                hyper_vec=hyper_vec, tuning_vec=tuning_vec, start_vec=start_mat[,i],
                                sample_alpha1 = sample_alpha1,
                                sample_alpha2 = sample_alpha2,
                                sample_alpha3 = sample_alpha3,
                                sample_kappa1 = sample_kappa1,
                                sample_kappa2 = sample_kappa2,
                                sample_kappa3 = sample_kappa3,
                                sample_beta1 = sample_beta1,
                                sample_beta2 = sample_beta2,
                                sample_beta3 = sample_beta3,
                                sample_betaD = sample_betaD,
                                sample_frail = sample_frail,
                                sample_theta = sample_theta,
                                accept_base = accept_base,
                                accept_frail = accept_frail,
                                accept_beta1 = accept_beta1,
                                accept_beta2 = accept_beta2,
                                accept_beta3 = accept_beta3,
                                LH_marg_mean_vec = LH_marg_mean_vec,
                                invLH_marg_mean_vec = invLH_marg_mean_vec,
                                sample_logLH_marg = sample_logLH_marg,
                                sample_logLHi_marg = sample_logLHi_marg,
                                LH_mean_vec = LH_mean_vec,
                                invLH_mean_vec = invLH_mean_vec,
                                sample_logLH = sample_logLH,
                                sample_logLHi = sample_logLHi,
                                move_vec=move_vec,
                                n_burnin=n_burnin, n_sample=n_sample, thin=thin,
                                frail_ind = as.integer(frailty),
                                nGam_save=nGam_save, nlogLHi_save=nlogLHi_save,
                                gh_nodes = gh_nodes, gh_weights = gh_weights)

    out_list[["samples"]][,i,1] <- sample_kappa1
    out_list[["samples"]][,i,2] <- sample_alpha1
    out_list[["samples"]][,i,3] <- sample_kappa2
    out_list[["samples"]][,i,4] <- sample_alpha2
    out_list[["samples"]][,i,5] <- sample_kappa3
    out_list[["samples"]][,i,6] <- sample_alpha3
    out_list[["samples"]][,i,7] <- sample_theta
    if(p1>0) out_list[["samples"]][,i,8:(7+p1)] <- t(sample_beta1)
    if(p2>0) out_list[["samples"]][,i,(8+p1):(7+p1+p2)] <- t(sample_beta2)
    if(p3>0) out_list[["samples"]][,i,(8+p1+p2):(7+p1+p2+p3)] <- t(sample_beta3)
    if(pD>0) out_list[["samples"]][,i,(8+p1+p2+p3):(7+p1+p2+p3+pD)] <- t(sample_betaD)
    out_list$diagnostics$logLH_mat[,i] <- sample_logLH
    out_list$diagnostics$LH_mean_mat[,i] <- LH_mean_vec
    out_list$diagnostics$invLH_mean_mat[,i] <- invLH_mean_vec
    if(frailty){
      out_list$diagnostics$logLH_marg_mat[,i] <- sample_logLH_marg
      out_list$diagnostics$LH_marg_mean_mat[,i] <- LH_marg_mean_vec
      out_list$diagnostics$invLH_marg_mean_mat[,i] <- invLH_marg_mean_vec
    }

    #save the gammas and the log-likelihood matrix as RDS objects.
    if(!is.null(frail_path) & nGam_save > 0){
      saveRDS(sample_frail, file=paste0(frail_path,"/frail_chain", i, ".RDS"))
    }
    #save log-likelihood contribution matrix
    if(!is.null(logLHi_path) & nlogLHi_save > 0){
      saveRDS(sample_logLHi, file=paste0(logLHi_path,"/logLHi_chain", i, ".RDS"))
      if(frailty){
        saveRDS(sample_logLHi_marg, file=paste0(logLHi_path,"/logLHi_marg_chain", i, ".RDS"))
      }
    }

    #reassign all the acceptance things later.
    out_list[["accept"]][[paste0("chain",i)]] <-
      list(accept_alpha1 = accept_base[1], accept_alpha2 = accept_base[2],
           accept_alpha3 = accept_base[3], accept_kappa1 = accept_base[4],
           accept_kappa2 = accept_base[5], accept_kappa3 = accept_base[6],
           accept_theta = accept_base[7], accept_frail = accept_frail,
           accept_beta1 = accept_beta1, accept_beta2 = accept_beta2,
           accept_beta3 = accept_beta3, accept_betaD = accept_base[8])
    out_list[["move"]][[paste0("chain",i)]] <- move_vec
  }

  #compute the final model diagnostics based on the sampled outputs
  out_list$diagnostics$dev <- -2*mean(out_list$diagnostics$logLH_mat)
  out_list$diagnostics$DIC = 2*out_list$diagnostics$dev + 2*sum(log(apply(out_list$diagnostics$LH_mean_mat, 1, mean)))
  out_list$diagnostics$LPML = -sum(log(apply(out_list$diagnostics$invLH_mean_mat, 1, mean)))
  if(frailty){
    out_list$diagnostics$dev_marg <- -2*mean(out_list$diagnostics$logLH_marg_mat)
    out_list$diagnostics$DIC_marg = 2*out_list$diagnostics$dev_marg + 2*sum(log(apply(out_list$diagnostics$LH_marg_mean_mat, 1, mean)))
    out_list$diagnostics$LPML_marg = -sum(log(apply(out_list$diagnostics$invLH_marg_mean_mat, 1, mean)))
  }

  class(out_list) <- "Bayes_HReg2"
  #for now, my plan is going to be to leverage the bayesplot package to make visuals
  return(out_list)

}




####EARLIER VERSIONS BELOW####


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
#' @param frailty Boolean indicating whether a gamma distributed subject-specific frailty should
#'   be included.
#' @param frail_path filepath for where to write frailty output
#' @param logLHi_path filepath for where to write individual likelihood contribution output
#'
#'
#' @return a list with outputs
#' @import Formula
#' @export
Bayes_SCR_logit <- function(Formula, data, na.action="na.fail", subset=NULL,
                            hyperParams, mcmcParams, n_chains, start_mat=NULL, frailty=TRUE,
                            frail_path=NULL, logLHi_path=NULL){
  # browser()

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

  ymax <- max(time2[[1]])

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
  XmatD <- cbind(logit_intercept=rep(1,sum(delta1)),
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
                           if(frailty) hyperParams$theta else c(0.7,0.7), #placeholders
                           if(is.null(hyperParams$logit$logit.m0)) numeric(pD) else hyperParams$logit$logit.m0,
                           if(is.null(hyperParams$logit$logit.P0diag)) numeric(pD) else hyperParams$logit$logit.P0diag
  ))
  tuning_vec <- as.vector(c(mhProp_alpha1_var=mcmcParams$tuning$mhProp_alphag_var[1],
                            mhProp_alpha2_var=mcmcParams$tuning$mhProp_alphag_var[2],
                            mhProp_alpha3_var=mcmcParams$tuning$mhProp_alphag_var[3],
                            mhProp_theta_var= if(frailty) mcmcParams$tuning$mhProp_theta_var else 0.1)) #placeholder

  ####SET MCMC SETTINGS####
  numReps=mcmcParams$run$numReps
  thin=mcmcParams$run$thin
  burninPerc=mcmcParams$run$burninPerc

  #set up where if it's less than 0, set to n, otherwise take what it is.
  nGam_save=if(mcmcParams$storage$nGam_save <0) n else mcmcParams$storage$nGam_save
  nlogLHi_save=if(mcmcParams$storage$nlogLHi_save <0) n else mcmcParams$storage$nlogLHi_save


  n_burnin <- numReps * burninPerc
  n_sample <- numReps - n_burnin
  n_store  <- numReps/thin * (1 - burninPerc)
  stopifnot(n_burnin %% 1 == 0)
  stopifnot(n_sample > 0)
  if(n_store %% 1 != 0){ stop("numReps * burninPerc  must be divisible by thin")}

  if(!is.null(frail_path)){
    dir.create(paste(frail_path), recursive = TRUE, showWarnings = FALSE)
  }
  if(!is.null(logLHi_path)){
    dir.create(paste(logLHi_path), recursive = TRUE, showWarnings = FALSE)
  }

  ####ASSIGN START VALUES####
  if(is.null(start_mat)){
    if(is.null(n_chains)){stop("either supply non-null start_list, or specify n_chains.")}
    theta_start_temp <- if(frailty) stats::runif(n_chains, 0.1, 1.1) else numeric(n_chains) #theta
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
      #if frailties included, draw them randomly, otherwise set them to 1.
      if(frailty) matrix(stats::rgamma(n*n_chains, 1/theta_start_temp, 1/theta_start_temp),ncol = n_chains,byrow=TRUE) else matrix(1,nrow=n,ncol=n_chains) #frailties (initialized by gamma, even if drawn from log-normal)
    )
    rownames(start_mat) <- c("kappa1","alpha1","kappa2","alpha2","kappa3","alpha3",
                             "theta", if(p1>0) paste0("beta1_",1:p1),
                             if(p2>0) paste0("beta2_",1:p2),
                             if(p3>0) paste0("beta3_",1:p3),
                             if(pD>0) paste0("betaD_",1:pD),
                             paste0("gamma",1:n))
  }

  gh_weights <- get_ghquad_pointsweights(n_quad = 15)$weights
  gh_nodes <- get_ghquad_pointsweights(n_quad = 15)$points

  #### PREALLOCATE OUTPUT LIST ####
  out_list <- list(
    #generate an array to store the resulting samples
    samples = array(dim = c(n_store, n_chains, 7 + p1 + p2 + p3 + pD), #even if frailty==FALSE, leave empty row for it.
                    dimnames = list(as.character(1:n_store),
                                    paste0("chain:",1:n_chains),
                                    rownames(start_mat)[1:(7 + p1 + p2 + p3 + pD)])),
    covnames = list(covNames1= if(p1 > 0) colnames(Xmat1) else NULL,
                    covNames2= if(p2 > 0) colnames(Xmat2) else NULL,
                    covNames3= if(p3 > 0) colnames(Xmat3) else NULL,
                    covNamesD= if(pD > 1) colnames(XmatD[,-1,drop=FALSE]) else NULL),
    #list final characteristics useful to export
    setup = list(Formula=Formula,
                 frailty=frailty,
                 nCov0 = c(2,2,2),
                 nCov = c(p1,p2,p3,pD),
                 hyper_vec = hyper_vec, start_mat = start_mat,
                 tuning_vec = tuning_vec,
                 nGam_save = nGam_save, nlogLHi_save = nlogLHi_save,
                 numReps = numReps, thin = thin,
                 burninPerc = burninPerc, n_store=n_store,
                 ymax=ymax, n=n,
                 frail_path = frail_path, logLHi_path = logLHi_path,
                 hz.type = "Weibull", model = "semi-Markov", nChain = n_chains,
                 mcmc_para = c(n_burnin=n_burnin,n_sample=n_sample,thin=thin)), #slight duplication but that's ok
    #empty list to keep "accept" values in.
    accept = vector(mode = "list", length = n_chains),
    move = vector(mode = "list", length = n_chains),
    diagnostics=list(dev=NA, DIC = NA, LPML = NA, dev_marg=NA, DIC_marg = NA, LPML_marg = NA,
                     #the mean log likelihood at each stored iteration
                     logLH_mat = matrix(data = NA,nrow = n_store, ncol = n_chains),
                     logLH_marg_mat = if(frailty) matrix(data = NA,nrow = n_store, ncol = n_chains) else NULL,
                     #the running mean of the likelihood and inverse likelihood for each subject across the samples
                     LH_mean_mat = matrix(data = NA,nrow = n, ncol = n_chains),
                     invLH_mean_mat = matrix(data = NA,nrow = n, ncol = n_chains),
                     LH_marg_mean_mat = if(frailty) matrix(data = NA,nrow = n, ncol = n_chains) else NULL,
                     invLH_marg_mean_mat = if(frailty) matrix(data = NA,nrow = n, ncol = n_chains) else NULL),
    class = c("Bayes_HReg2", "IDlogit", "Ind", "WB")
  )
  names(out_list$accept) <- paste0("chain",1:n_chains)

  #### RUN SAMPLER ####
  # #TODO: parallelize this loop (I know it can be done!)
  for(i in 1:n_chains){
    print(paste0("Chain: ", i))

    #initialize all of the vectors of individually sampled parameters
    sample_alpha1 <- numeric(n_store)
    sample_alpha2 <- numeric(n_store)
    sample_alpha3 <- numeric(n_store)
    sample_kappa1 <- numeric(n_store)
    sample_kappa2 <- numeric(n_store)
    sample_kappa3 <- numeric(n_store)
    sample_theta <- numeric(n_store)

    #initialize all of the matrices of regression parameters
    sample_beta1 <- matrix(data=0,nrow=p1,ncol=n_store)
    sample_beta2 <- matrix(data=0,nrow=p2,ncol=n_store)
    sample_beta3 <- matrix(data=0,nrow=p3,ncol=n_store)
    sample_betaD <- matrix(data=0,nrow=pD,ncol=n_store)

    #initialize matrices to store frailties and log-likelihood, as needed
    sample_frail <- matrix(data=0,nrow=nGam_save,ncol=n_store)
    sample_logLHi <- matrix(data=0,nrow=nlogLHi_save,ncol=n_store)
    sample_logLHi_marg <- if(frailty) matrix(data=0,nrow=nlogLHi_save,ncol=n_store) else matrix(nrow=0,ncol=0)

    #initialize integers counting acceptance of individually sampled parameters
    accept_base <- numeric(8)

    #initialize vectors counting acceptance of sampled vectors
    accept_frail <- numeric(n)
    accept_beta1 <- numeric(p1)
    accept_beta2 <- numeric(p2)
    accept_beta3 <- numeric(p3)

    #initialize vectors for diagnostics
    LH_marg_mean_vec <- if(frailty) numeric(n) else numeric(0)
    invLH_marg_mean_vec <- if(frailty) numeric(n) else numeric(0)
    sample_logLH_marg <- if(frailty) numeric(n_store) else numeric(0)

    LH_mean_vec <- numeric(n)
    invLH_mean_vec <- numeric(n)
    sample_logLH <- numeric(n_store)

    #initialize move vector
    move_vec <- numeric(numReps)

    mcmcRet <- WeibSCRlogitmcmc(y1=y1, y_sm=y_sm,
                                delta1=delta1, #indicator for nonterminal event followed by immediate death
                                delta1noD=delta1noD, #indicator for nonterminal event not followed by immediate death
                                delta_cr=delta_cr, #terminal first event indicator
                                delta_sm=delta_sm, #terminal event indicator in group with nonterminal first event and no immediate terminal
                                delta1D_sub=delta1D_sub, #immediate terminal event indicator in group with nonterminal first event
                                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, XmatD=XmatD, eps=0,
                                hyper_vec=hyper_vec, tuning_vec=tuning_vec, start_vec=start_mat[,i],
                                sample_alpha1 = sample_alpha1,
                                sample_alpha2 = sample_alpha2,
                                sample_alpha3 = sample_alpha3,
                                sample_kappa1 = sample_kappa1,
                                sample_kappa2 = sample_kappa2,
                                sample_kappa3 = sample_kappa3,
                                sample_beta1 = sample_beta1,
                                sample_beta2 = sample_beta2,
                                sample_beta3 = sample_beta3,
                                sample_betaD = sample_betaD,
                                sample_frail = sample_frail,
                                sample_theta = sample_theta,
                                accept_base = accept_base,
                                accept_frail = accept_frail,
                                accept_beta1 = accept_beta1,
                                accept_beta2 = accept_beta2,
                                accept_beta3 = accept_beta3,
                                LH_marg_mean_vec = LH_marg_mean_vec,
                                invLH_marg_mean_vec = invLH_marg_mean_vec,
                                sample_logLH_marg = sample_logLH_marg,
                                sample_logLHi_marg = sample_logLHi_marg,
                                LH_mean_vec = LH_mean_vec,
                                invLH_mean_vec = invLH_mean_vec,
                                sample_logLH = sample_logLH,
                                sample_logLHi = sample_logLHi,
                                move_vec=move_vec,
                                n_burnin=n_burnin, n_sample=n_sample, thin=thin,
                                frail_ind = as.integer(frailty),
                                nGam_save=nGam_save, nlogLHi_save=nlogLHi_save,
                                gh_nodes = gh_nodes, gh_weights = gh_weights)

    out_list[["samples"]][,i,1] <- sample_kappa1
    out_list[["samples"]][,i,2] <- sample_alpha1
    out_list[["samples"]][,i,3] <- sample_kappa2
    out_list[["samples"]][,i,4] <- sample_alpha2
    out_list[["samples"]][,i,5] <- sample_kappa3
    out_list[["samples"]][,i,6] <- sample_alpha3
    out_list[["samples"]][,i,7] <- sample_theta
    if(p1>0) out_list[["samples"]][,i,8:(7+p1)] <- t(sample_beta1)
    if(p2>0) out_list[["samples"]][,i,(8+p1):(7+p1+p2)] <- t(sample_beta2)
    if(p3>0) out_list[["samples"]][,i,(8+p1+p2):(7+p1+p2+p3)] <- t(sample_beta3)
    if(pD>0) out_list[["samples"]][,i,(8+p1+p2+p3):(7+p1+p2+p3+pD)] <- t(sample_betaD)
    out_list$diagnostics$logLH_mat[,i] <- sample_logLH
    out_list$diagnostics$LH_mean_mat[,i] <- LH_mean_vec
    out_list$diagnostics$invLH_mean_mat[,i] <- invLH_mean_vec
    if(frailty){
      out_list$diagnostics$logLH_marg_mat[,i] <- sample_logLH_marg
      out_list$diagnostics$LH_marg_mean_mat[,i] <- LH_marg_mean_vec
      out_list$diagnostics$invLH_marg_mean_mat[,i] <- invLH_marg_mean_vec
    }

    #save the gammas and the log-likelihood matrix as RDS objects.
    if(!is.null(frail_path) & nGam_save > 0){
      saveRDS(sample_frail, file=paste0(frail_path,"/frail_chain", i, ".RDS"))
    }
    #save log-likelihood contribution matrix
    if(!is.null(logLHi_path) & nlogLHi_save > 0){
      saveRDS(sample_logLHi, file=paste0(logLHi_path,"/logLHi_chain", i, ".RDS"))
      if(frailty){
        saveRDS(sample_logLHi_marg, file=paste0(logLHi_path,"/logLHi_marg_chain", i, ".RDS"))
      }
    }

    #reassign all the acceptance things later.
    out_list[["accept"]][[paste0("chain",i)]] <-
      list(accept_alpha1 = accept_base[1], accept_alpha2 = accept_base[2],
           accept_alpha3 = accept_base[3], accept_kappa1 = accept_base[4],
           accept_kappa2 = accept_base[5], accept_kappa3 = accept_base[6],
           accept_theta = accept_base[7], accept_frail = accept_frail,
           accept_beta1 = accept_beta1, accept_beta2 = accept_beta2,
           accept_beta3 = accept_beta3, accept_betaD = accept_base[8])
    out_list[["move"]][[paste0("chain",i)]] <- move_vec
  }

  #compute the final model diagnostics based on the sampled outputs
  out_list$diagnostics$dev <- -2*mean(out_list$diagnostics$logLH_mat)
  out_list$diagnostics$DIC = 2*out_list$diagnostics$dev + 2*sum(log(apply(out_list$diagnostics$LH_mean_mat, 1, mean)))
  out_list$diagnostics$LPML = -sum(log(apply(out_list$diagnostics$invLH_mean_mat, 1, mean)))
  if(frailty){
    out_list$diagnostics$dev_marg <- -2*mean(out_list$diagnostics$logLH_marg_mat)
    out_list$diagnostics$DIC_marg = 2*out_list$diagnostics$dev_marg + 2*sum(log(apply(out_list$diagnostics$LH_marg_mean_mat, 1, mean)))
    out_list$diagnostics$LPML_marg = -sum(log(apply(out_list$diagnostics$invLH_marg_mean_mat, 1, mean)))
  }

  class(out_list) <- "Bayes_HReg2"
  #for now, my plan is going to be to leverage the bayesplot package to make visuals
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
