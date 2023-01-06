#' Bayesian Illness-Death Model with Weibull Baselines and Semi-Markov structure
#'
#' @inheritParams Bayes_SCR_logit
#'
#' @return a list with outputs
#' @import Formula
#' @export
Bayes_SCR <- function(Formula, data, na.action="na.fail", subset=NULL,
                      hyperParams, mcmcParams, n_chains, start_mat=NULL, frailty=TRUE, frail_dist="gamma",
                      frail_path=NULL, logLHi_path=NULL){
  # browser()

  frail_dist <- tolower(frail_dist)

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

  #indicator for first event being non-terminal event
  delta1 <- time1[[2]]
  #subtly different indicator for first event being non-terminal event with nonzero sojourn time after
  delta1noD <- as.numeric( delta1 & time2[[1]] - time1[[1]] > 0)

  #time of first event/censoring
  y1 <- time1[[1]]
  #indicator for first event being terminal event (nonterminal is competing risk)
  delta_cr <- ifelse(time1[[1]] < time2[[1]], 0, time2[[2]])

  #time of terminal event after non-terminal event excl. obs with immediate event (sojourn time 0)
  y_sm <- (time2[[1]] - time1[[1]])[delta1noD==1]
  #indicator of terminal event after non-terminal event excl. obs with immediate event (sojourn time 0)
  delta_sm <- time2[[2]][delta1noD==1]

  ymax <- max(time2[[1]])

  #ensure that number of obs in h3 is correctly accounted for
  stopifnot(length(y_sm) == sum(delta1noD))


  #Create covariate matrices for each of three transition hazards
  Xmat1 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=1),
                                        data=data))
  Xmat2 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=2),
                                        data=data))
  #subsetted to just those who experience the non-terminal event
  Xmat3 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=3),
                                        data=data))[delta1noD==1,,drop=FALSE]
  p1 <- ncol(Xmat1)
  p2 <- ncol(Xmat2)
  p3 <- ncol(Xmat3)
  n <- length(y1)

  ####SET HYPERPARAMETERS####
  hyper_vec <- as.vector(c(hyperParams$WB$WB.ab1, hyperParams$WB$WB.ab2, hyperParams$WB$WB.ab3,
                             hyperParams$WB$WB.cd1, hyperParams$WB$WB.cd2, hyperParams$WB$WB.cd3,
                           if(frailty) hyperParams$theta else c(0.7,0.7))) #placeholders))
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

  #we're gonna write these values separately because there's so many of them
  #so make sure there's a directory to write to.
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
      #if frailties included, draw them randomly, otherwise set them to 1.
      if(frailty) matrix(stats::rgamma(n*n_chains, 1/theta_start_temp, 1/theta_start_temp),ncol = n_chains,byrow=TRUE) else matrix(1,nrow=n,ncol=n_chains) #frailties (initialized by gamma, even if drawn from log-normal)
    )
    rownames(start_mat) <- c("kappa1","alpha1","kappa2","alpha2","kappa3","alpha3",
                             "theta", if(p1>0) paste0("beta1_",1:p1),
                             if(p2>0) paste0("beta2_",1:p2),
                             if(p3>0) paste0("beta3_",1:p3),
                             paste0("gamma",1:n))
  }

  gh_weights <- get_ghquad_pointsweights(n_quad = 15)$weights
  gh_nodes <- get_ghquad_pointsweights(n_quad = 15)$points

  #### PREALLOCATE OUTPUT LIST ####
  out_list <- list(
    #generate an array to store the resulting samples
    samples = array(dim = c(n_store, n_chains, 7 + p1 + p2 + p3), #even if frailty==FALSE, leave empty row for it.
                    dimnames = list(as.character(1:n_store),
                                    paste0("chain:",1:n_chains),
                                    rownames(start_mat)[1:(7 + p1 + p2 + p3)])),
    covnames = list(covNames1= if(p1 > 0) colnames(Xmat1) else NULL,
                    covNames2= if(p2 > 0) colnames(Xmat2) else NULL,
                    covNames3= if(p3 > 0) colnames(Xmat3) else NULL),
    #list final characteristics useful to export
    setup = list(Formula=Formula,
                 frailty=frailty,
                 nCov0 = c(2,2,2),
                 nCov = c(p1,p2,p3),
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
    class = c("Bayes_HReg2", "ID", "Ind", "WB")
  )
  names(out_list$accept) <- paste0("chain",1:n_chains)
  class(out_list) <- "Bayes_HReg2"

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

    #initialize matrices to store frailties and log-likelihood, as needed
    sample_frail <- matrix(data=0,nrow=nGam_save,ncol=n_store)
    sample_logLHi <- matrix(data=0,nrow=nlogLHi_save,ncol=n_store)
    sample_logLHi_marg <- if(frailty) matrix(data=0,nrow=nlogLHi_save,ncol=n_store) else matrix(nrow=0,ncol=0)

    #initialize integers counting acceptance of individually sampled parameters
    accept_base <- numeric(7)

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

    mcmcRet <- WeibSCRmcmc(y1=y1, y_sm=y_sm,
                                delta1=delta1, #indicator for nonterminal event followed by immediate death
                                delta1noD=delta1noD, #indicator for nonterminal event not followed by immediate death
                                delta_cr=delta_cr, #terminal first event indicator
                                delta_sm=delta_sm, #terminal event indicator in group with nonterminal first event and no immediate terminal
                                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
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
                                frail_dist = "gamma",
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
           accept_beta3 = accept_beta3)
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

  #for now, my plan is going to be to leverage the bayesplot package to make visuals
  return(out_list)

}
