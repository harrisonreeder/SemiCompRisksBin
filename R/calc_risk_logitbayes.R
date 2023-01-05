#' Calculate absolute risk profiles
#'
#' This function calculates absolute risk profiles under weibull baseline hazard specifications.
#'
#' @inheritParams calc_risk
#' @param XmatD Numeric matrices with \eqn{n} rows and \eqn{q_1,q_2,q_3} columns containing covariates.
#' @param logit_ind Boolean for whether to include logit submodel.
#' @param logit_tv String indicating whether there is an effect of t1 on the logit model.
#' @param logit_tv_knots for piecewise effect of t1 in h3, these are the knots at which the effect jumps
#' @param beta2frail,beta3frail,betaDfrail coefficients on the log-frailty in each hazard, only used for conditional prediction
#' @param h3tv_knots for piecewise effect of t1 in h3, these are the knots at which the effect jumps
#' @param n_quad Scalar for number of Gaussian quadrature points used to evaluate numerical integral of B-spline.
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 4 element row matrix
#'   of probabilities. If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix of probabilities.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
calc_risk_logitbayes <- function(para, Xmat1, Xmat2, Xmat3, XmatD, hazard,
                            t_cutoff, tol=1e-3, frailty=TRUE,
                            beta2frail=1, beta3frail=1, betaDfrail=1, #temporary, just get them in the mix
                            type="marginal", gamma=1,model="semi-markov",
                            h3_tv="none",h3tv_knots=NULL,
                            logit_ind=TRUE,logit_tv="none",logit_tv_knots=NULL,
                            n_quad = 15){
  # browser()
  ##TO START, EXTRACT POINT ESTIMATES OF ALL PARAMETERS FROM MODEL OBJECT##
  ##*********************************************************************##

  n <- max(1,nrow(Xmat1),nrow(Xmat2),nrow(Xmat3))
  t_length <- length(t_cutoff)
  #standardize namings for the use of "switch" below
  #for now, restrict to weibull
  stopifnot(tolower(hazard) %in% c("wb","weibull"))
  hazard <- switch(tolower(hazard),
                   wb="weibull",weibull="weibull")
  stopifnot(tolower(type) %in% c("c","conditional","m","marginal"))
  type <- switch(tolower(type),
                 c="conditional",conditional="conditional",m="marginal",marginal="marginal")
  #for now, restrict to semi-markov
  stopifnot(tolower(model) %in% c("sm","semi-markov"))
  model <- switch(tolower(model),
                  sm="semi-markov","semi-markov"="semi-markov",m="markov",markov="markov")

  #right now, we assume that h3 and logit have same knots and denote that by leaving logit_tv_knots blank,
  #so throw an error if not
  if(!is.null(logit_tv_knots)){warning("logit_tv_knots isn't actually used, we assume it's the same as h3_tv_knots")}

  #set number of baseline parameters to be found in para vector
  nP01 <- nP02 <- nP03 <- 2
  #always assume theta parameter is present even for non-frailty model
  nP0 <- nP01 + nP02 + nP03 + 1
  if(frailty){
    theta <- para[nP0] #for bayesian model, theta is not logged, so just store directly
    if(type=="conditional" & length(gamma)==1){
      gamma <- rep(gamma,n)
    }
  } else{ #if no frailty, can't be "marginalized" over frailties
    theta <- 0
    type <- "conditional"
    gamma <- rep(1,n)
  }

  #set linear predictors corresponding to each cause-specific hazard
  if(!is.null(Xmat1) && !(ncol(Xmat1)==0)){
    nP1 <- ncol(Xmat1)
    beta1 <- para[(1+nP0):(nP0+nP1)]
    eta1 <- as.vector(Xmat1 %*% beta1)
  } else{
    nP1 <- 0
    eta1 <- 0
  }
  if(!is.null(Xmat2) && !(ncol(Xmat2)==0)){
    nP2 <- ncol(Xmat2)
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    eta2 <- as.vector(Xmat2 %*% beta2)
  } else{
    nP2 <- 0
    eta2 <- 0
  }
  if(!is.null(Xmat3) && !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }

  #specify different forms by which t1 can be incorporated into h3
  #assumes that effect of t1 is at the end of the set of beta3 coefficints
  if(tolower(h3_tv) == "linear"){
    if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
    n_h3_tv <- 1
    beta3_tv_linear <- para[(1+nP0+nP1+nP2+nP3):(nP0+nP1+nP2+nP3+n_h3_tv)]
  } else if(tolower(h3_tv) %in% c("pw","piecewise")){
    if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
    stopifnot(!is.null(h3tv_knots))
    if(h3tv_knots[1] != 0){h3tv_knots <- c(0,h3tv_knots)}
    if(utils::tail(h3tv_knots, n=1) != Inf){h3tv_knots <- c(h3tv_knots,Inf)}
    n_h3_tv <- length(h3tv_knots) - 2
    beta3_tv <- c(0,para[(1+nP0+nP1+nP2+nP3):(nP0+nP1+nP2+nP3+n_h3_tv)])
    beta3_tv_linear <- 0
  } else{
    n_h3_tv <- 0
    beta3_tv_linear <- 0
  }

  #coefficients related to logit model
  if(logit_ind){
    beta0D <- para[(1+nP0+nP1+nP2+nP3+n_h3_tv):(nP0+nP1+nP2+nP3+n_h3_tv+1)]
    #now, add in the component corresponding with the logit model
    if(!is.null(XmatD) && !(ncol(XmatD)==0)){
      nPD <- ncol(XmatD)
      betaD <- para[(1+nP0+nP1+nP2+nP3+n_h3_tv+1):(nP0+nP1+nP2+nP3+n_h3_tv+1+nPD)]
      etaD <- as.vector(XmatD %*% betaD)
    } else{
      nPD <- 0
      etaD <- 0
    }

    #specify different forms by which t1 can be incorporated into the logit model
    #assumes that effect of t1 is at the end of the set of betaD coefficients
    if(tolower(logit_tv) == "linear"){
      n_logit_tv <- 1
      betaD_tv_linear <- para[(1+nP0+nP1+nP2+nP3+n_h3_tv+1+nPD):(nP0+nP1+nP2+nP3+n_h3_tv+1+nPD+n_logit_tv)]
    } else if(tolower(logit_tv) %in% c("pw","piecewise")){
      if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
      stopifnot(!is.null(logit_tv_knots))
      if(logit_tv_knots[1] != 0){logit_tv_knots <- c(0,logit_tv_knots)}
      if(utils::tail(logit_tv_knots, n=1) != Inf){logit_tv_knots <- c(logit_tv_knots,Inf)}
      n_logit_tv <- length(logit_tv_knots) - 2
      betaD_tv <- c(0,para[(1+nP0+nP1+nP2+nP3+n_h3_tv+1+nPD):(nP0+nP1+nP2+nP3+n_h3_tv+1+nPD+n_logit_tv)])
      betaD_tv_linear <- 0
    } else{
      n_logit_tv <- 0
      betaD_tv_linear <- 0
    }
    nP_logit <- 1+nPD+n_logit_tv
  } else{
    nP_logit <- 0
  }

  #if the size of the parameter vector doesn't match the expected size, throw a fuss
  if(!(length(para) == nP0 + nP1 + nP2 + nP3 + n_h3_tv + logit_ind * nP_logit)){
    warning("provided parameter vector is different length than expected...")
  }

  ##Set up the hazard functions##
  ##***************************##
  alpha1=para[2]
  alpha2=para[4]
  alpha3=para[6]
  kappa1=para[1]
  kappa2=para[3]
  kappa3=para[5]

  #first, compute some helper quantities
  h1_const=alpha1 * kappa1 * exp(eta1)
  h2_const=alpha2 * kappa2 * exp(eta2)
  h3_const=alpha3 * kappa3 * exp(eta3)
  H1_const=kappa1 * exp(as.vector(eta1))
  H2_const=kappa2 * exp(as.vector(eta2))
  H3_const=kappa3 * exp(as.vector(eta3))
  alpha1_m1=alpha1 - 1
  alpha2_m1=alpha2 - 1
  alpha3_m1=alpha3 - 1

  ##set up the components for numerically integrating over the frailties
  gh_adj_nodes <- get_ghquad_pointsweights(n_quad=n_quad)$points * sqrt(2 * theta)
  gh_weights <- get_ghquad_pointsweights(n_quad=n_quad)$weights

  #now, in the below, define an additional version of each function we integrate over
  #to compute `marginally`

  if(logit_ind){
    pi_vec <- function(t1,index,gamma,betaD_tv_const,betaD_tv_lin){
      if(length(gamma) > 1){
        return(stats::plogis(q = beta0D + etaD[index] + betaD_tv_const + betaD_tv_lin * t1 + betaDfrail * log(gamma[index])))
      } else{
        return(stats::plogis(q = beta0D + etaD[index] + betaD_tv_const + betaD_tv_lin * t1 + betaDfrail * log(gamma)))
      }
    }
  } else{
    pi_vec <- function(t1,index,gamma,betaD_tv_const,betaD_tv_lin){ return(0)}
  }


  ##*************************************##
  ## Calculating predicted risk profiles ##
  ##*************************************##

  ## First, write functions that compute the integrand,
  ## which we feed into integration function
  ##***********************************************##

  #the univariate function when T1=infinity
  #aka, the integrand needed to generate the probability of just the terminal event occurring.
  f_t2 <- switch(type,
           marginal=function(t2, index){
             temp <- 0
             for(x in 1:n_quad){
               temp2 <- exp(gh_adj_nodes[x] * beta2frail) *
                 h2_const[index] * (t2)^alpha2_m1 *
                 exp(-exp(gh_adj_nodes[x]) * H1_const[index] * (t2)^alpha1 -
                       exp(gh_adj_nodes[x] * beta2frail) * H2_const[index] * (t2)^alpha2)
               temp <- temp + temp2 * gh_weights[x]
             }
             temp / sqrt(pi)
           },
           conditional=function(t2, index){
             gamma[index]^beta2frail * h2_const[index] * (t2)^alpha2_m1 *
               exp(-gamma[index]*H1_const[index] * (t2)^alpha1 -
                     gamma[index]^beta2frail * H2_const[index] * (t2)^alpha2)
           })

  #next, the different regions of the joint density on the upper triangle
  # here is the generic (semi-markov) formula if we pre-integrate t2 from u to v
  # \int \lambda_1(t_1)\exp(-[\Lambda_1(t_1)+\Lambda_2(t_1)]) \left[\exp(-\Lambda_3(u-t_1)) - \exp(-\Lambda_3(v-t_1))\right] dt_1
  # here is the generic (markov) formula if we pre-integrate t2 from u to v
  # \int \lambda_1(t_1)\exp(-[\Lambda_1(t_1)+\Lambda_2(t_1)-\Lambda_3(t_1)]) \left[\exp(-\Lambda_3(u)) - \exp(-\Lambda_3(v))\right] dt_1

  #function of t1 if we pre-integrate t2 from t_cutoff to infinity
  #aka, integrand needed to generate probability of just the non-terminal event occurring by t_cutoff
  f_joint_t1_nonTerm <- switch(type,
    marginal=function(t1,t_cutoff,index,beta3_tv_const=0,beta3_tv_lin=0, betaD_tv_const=0,betaD_tv_lin=0){
      # browser()
      temp <- 0
      for(x in 1:n_quad){
        temp2 <- exp(gh_adj_nodes[x]) * h1_const[index] * t1^alpha1_m1 *
          (1 - pi_vec(t1 = t1,index = index,gamma = exp(gh_adj_nodes[x]),
                      betaD_tv_const = betaD_tv_const,betaD_tv_lin = betaD_tv_lin)) *
          exp(-exp(gh_adj_nodes[x]) * H1_const[index] * t1^alpha1 -
                exp(gh_adj_nodes[x] * beta2frail) * H2_const[index] * t1^alpha2 -
                exp(gh_adj_nodes[x] * beta3frail + beta3_tv_const + beta3_tv_lin * t1) *
                H3_const[index] * (t_cutoff - t1)^alpha3)
        temp <- temp + temp2 * gh_weights[x]
      }
      temp / sqrt(pi)
    },
    conditional=function(t1,t_cutoff,index,beta3_tv_const=0,beta3_tv_lin=0, betaD_tv_const=0,betaD_tv_lin=0){
      # browser()
      H3_temp <- H3_const[index] * (t_cutoff - t1)^alpha3 * exp(beta3_tv_const + beta3_tv_lin * t1)
      #return the right value corresponding to the hazard and type
      #ASSUME SEMI-MARKOV, CONDITIONAL, AND WEIBULL
      gamma[index] * h1_const[index] * t1^alpha1_m1 *
        (1 - pi_vec(t1 = t1,index = index,gamma = gamma,
                    betaD_tv_const = betaD_tv_const,betaD_tv_lin = betaD_tv_lin)) *
        exp(-gamma[index] * H1_const[index] * t1^alpha1 -
              gamma[index]^beta2frail * H2_const[index] * t1^alpha2 -
              gamma[index]^beta3frail * H3_temp)
    }
  )

  #function of t1 if we pre-integrate t2 from t1 to t_cutoff
  #aka, integrand needed to generate probability of both events occurring by t_cutoff
  f_joint_t1_both <- switch(type,
    marginal=function(t1,t_cutoff,index,beta3_tv_const=0,beta3_tv_lin=0,betaD_tv_const=0,betaD_tv_lin=0){
      temp <- 0
      for(x in 1:n_quad){
        temp2 <- exp(gh_adj_nodes[x]) * h1_const[index] * t1^alpha1_m1 *
          (1 - pi_vec(t1 = t1,index = index, gamma = exp(gh_adj_nodes[x]),
                      betaD_tv_const = betaD_tv_const,betaD_tv_lin = betaD_tv_lin)) *
          exp( -exp(gh_adj_nodes[x]) * H1_const[index] * t1^alpha1 -
                 exp(gh_adj_nodes[x] * beta2frail) * H2_const[index] * t1^alpha2 ) *
          ( 1 - exp(-exp(gh_adj_nodes[x] * beta3frail + beta3_tv_const + beta3_tv_lin * t1) *
                      H3_const[index] * (t_cutoff - t1)^alpha3))
        temp <- temp + temp2 * gh_weights[x]
      }
      temp / sqrt(pi)
    },
    conditional=function(t1,t_cutoff,index,beta3_tv_const=0,beta3_tv_lin=0,betaD_tv_const=0,betaD_tv_lin=0){
    H3_temp <- H3_const[index] * (t_cutoff - t1)^alpha3 * exp(beta3_tv_const + beta3_tv_lin * t1)
    #return the right value corresponding to the hazard and type
    #ASSUME SEMI-MARKOV, CONDITIONAL, AND WEIBULL
    gamma[index] * h1_const[index] * t1^alpha1_m1 *
      (1 - pi_vec(t1 = t1,index = index,gamma = gamma,betaD_tv_const = betaD_tv_const,betaD_tv_lin = betaD_tv_lin)) *
      exp( -gamma[index]*H1_const[index] * t1^alpha1 -
             gamma[index]^beta2frail * H2_const[index] * t1^alpha2 ) *
      ( 1 - exp( -gamma[index]^beta3frail * H3_temp))
    }
  )

  #function of t1 if we pre-integrate t2 from t1 to t_cutoff
  #aka, integrand needed to generate probability of both events occurring by t_cutoff
  f_joint_t1_both_inst <- switch(type,
    marginal=function(t1, index, beta3_tv_const=0, beta3_tv_lin=0, betaD_tv_const=0, betaD_tv_lin=0){
      temp <- 0
      for(x in 1:n_quad){
        temp2 <- exp(gh_adj_nodes[x]) * h1_const[index] * t1^alpha1_m1 *
          pi_vec(t1 = t1,index = index,gamma = exp(gh_adj_nodes[x]),
                 betaD_tv_const = betaD_tv_const,betaD_tv_lin = betaD_tv_lin) *
          exp( -exp(gh_adj_nodes[x]) * H1_const[index] * t1^alpha1 -
                 exp(gh_adj_nodes[x] * beta2frail) * H2_const[index] * t1^alpha2 )
        temp <- temp + temp2 * gh_weights[x]
      }
      temp / sqrt(pi)
    },
    conditional=function(t1, index, beta3_tv_const=0, beta3_tv_lin=0, betaD_tv_const=0, betaD_tv_lin=0){
    #return the right value corresponding to the hazard and type
    #ASSUME SEMI-MARKOV, CONDITIONAL, AND WEIBULL
    gamma[index] * h1_const[index] * t1^alpha1_m1 *
      pi_vec(t1 = t1,index = index,gamma = gamma, betaD_tv_const = betaD_tv_const,betaD_tv_lin = betaD_tv_lin) *
      exp( -gamma[index] * H1_const[index] * t1^alpha1 -
             gamma[index]^beta2frail * H2_const[index] * t1^alpha2 )
    }
  )


  ##finally, p_neither has a closed form, so we can write a function for it directly##
  #this derivation is actually identical to the "no event" likelihood contribution
  p_neither_func <- switch(type,
    marginal=function(t2){
      temp <- 0
      for(x in 1:n_quad){
        temp2 <- exp(-exp(gh_adj_nodes[x]) * H1_const * (t2)^alpha1 -
                       exp(gh_adj_nodes[x]*beta2frail) * H2_const * (t2)^alpha2)
        temp <- temp + temp2 * gh_weights[x]
      }
      temp / sqrt(pi)
    },
    conditional=function(t2){
    exp(-gamma * H1_const * (t2)^alpha1 -
          gamma^beta2frail * H2_const * (t2)^alpha2)
    }
  )

  ##ACTUALLY COMPUTING THE PROBABILITIES##
  ##************************************##

  n_cat <- if(logit_ind) 5 else 4
  cat_labs <- c("p_ntonly",
                if(logit_ind) c("p_both_noinst","p_both_inst") else "p_both",
                "p_tonly","p_neither")

  #this function allows inputs with multiple subjects, multiple time points, or both
  #therefore, we need to create the right data structure to contain the output.
  if(n > 1){
    if(t_length > 1){
      out_mat <- array(dim=c(t_length,n_cat,n),dimnames = list(paste0("t",t_cutoff),cat_labs,paste0("i",1:n)))
    } else{
      out_mat <- matrix(nrow=n,ncol=n_cat,dimnames = list(paste0("i",1:n),cat_labs))
    }
  } else{
    out_mat <- matrix(nrow=t_length,ncol=n_cat,dimnames = list(paste0("t",t_cutoff),cat_labs))
  }

  #loop through each time point, and compute the predicted probability at that time point for all subjects
  #each probability is an n-length vector

  #TO-DO: MAKE THIS GO INCREMENTALLY, WHICH I THINK WILL BE FASTER / MORE STABLE?

  for(t_ind in 1:t_length){
    t_temp <- t_cutoff[t_ind]
    p_tonly <- sapply(1:n,function(x){tryCatch(stats::integrate(f_t2, lower=0, upper=t_temp, index=x)$value,
                                               error=function(cnd){return(NA)}) })
    p_neither <- p_neither_func(t2=t_temp)

    #now, we have to be careful computing the probabilities that include h3 because it may depend on t1
    # p_both_start <- p_ntonly_start <- rep(0,n)

    if(h3_tv %in% c("piecewise")){
      #a piecewise effect cannot be integrated in one go, because it is discontinuous
      #instead it must be divided into constant regions, integrated one region at a time and summed up
      curr_interval <- findInterval(x = t_temp,h3tv_knots,left.open = TRUE)
      #this is now a vector starting at 0, with elements at each change point up to the target time t_temp
      tv_knots_temp <- c(h3tv_knots[1:curr_interval],t_temp)
      if(t_temp == 0){ tv_knots_temp <- c(0,0)}

      #loop through the regions of constant effect, and add them together
      p_both <- p_ntonly <- rep(0,n)
      if(logit_ind) p_both_inst <- rep(0,n)
      for(i in 1:(length(tv_knots_temp)-1)){
        p_both <- p_both + sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both, lower=tv_knots_temp[i], upper=tv_knots_temp[i+1],
                                                                            t_cutoff=t_temp, index=x,
                                                                            beta3_tv_const=beta3_tv[i], beta3_tv_lin=0,
                                                                            betaD_tv_const=betaD_tv[i],betaD_tv_lin=0)$value,
                                                           error=function(cnd){return(NA)}) })

        p_ntonly <- p_ntonly + sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_nonTerm, lower=tv_knots_temp[i], upper=tv_knots_temp[i+1],
                                                                                t_cutoff=t_temp, index=x,
                                                                                beta3_tv_const=beta3_tv[i], beta3_tv_lin=0,
                                                                                betaD_tv_const=betaD_tv[i],betaD_tv_lin=0)$value,
                                                               error=function(cnd){return(NA)}) })
        if(logit_ind){
          p_both_inst <- p_both_inst + sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both_inst, lower=tv_knots_temp[i], upper=tv_knots_temp[i+1],
                                                                                        index=x, beta3_tv_const=beta3_tv[i], beta3_tv_lin=0,
                                                                                        betaD_tv_const=betaD_tv[i],betaD_tv_lin=0)$value,
                                                                       error=function(cnd){return(NA)}) })
        }
      }
    } else{ #just normal, time-invariant prediction or a linear effect of time
      p_both <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both, lower=0, upper=t_temp,
                                                                 t_cutoff=t_temp, index=x,
                                                                 beta3_tv_const=0,beta3_tv_lin=beta3_tv_linear,
                                                                 betaD_tv_const=0,betaD_tv_lin=betaD_tv_linear)$value,
                                                error=function(cnd){return(NA)}) })
      p_ntonly <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_nonTerm, lower=0, upper=t_temp,
                                                                   t_cutoff=t_temp, index=x,
                                                                   beta3_tv_const=0,beta3_tv_lin=beta3_tv_linear,
                                                                   betaD_tv_const=0,betaD_tv_lin=betaD_tv_linear)$value,
                                                  error=function(cnd){return(NA)}) })
      if(logit_ind){
        p_both_inst <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both_inst, lower=0, upper=t_temp,
                                                                        index=x, beta3_tv_const=0,beta3_tv_lin=beta3_tv_linear,
                                                                        betaD_tv_const=0,betaD_tv_lin=betaD_tv_linear)$value,
                                                       error=function(cnd){return(NA)}) })
      }
    }

    out_temp <- cbind(p_ntonly=p_ntonly,
                      if(logit_ind) cbind(p_both_noinst=p_both,p_both_inst=p_both_inst) else cbind(p_both=p_both),
                      p_tonly=p_tonly,
                      p_neither=p_neither)

    #I noticed that sometimes, if exactly one category has an NA, then we could back out the value
    #from the other categories. However, we don't want any to be negative.
    #I will also supply a warning if any of the rows are way off from 1.
    # out_temp <- t(apply(out_temp,1,
    #                     function(x){
    #                       if(sum(is.na(x))==1){
    #                         x[is.na(x)]<- max(1-sum(x,na.rm=TRUE),0)
    #                       }
    #                       return(x)}))
    if(any(is.na(out_temp)) | any(abs(1-rowSums(out_temp))>tol)){
      warning(paste0("some predicted probabilities at time ", t_temp," do not sum to within",tol,"of 1."))
    }

    if(n > 1){
      if(t_length > 1){
        out_mat[t_ind,,] <- t(out_temp)
      } else{
        out_mat <- out_temp
      }
    } else{
      out_mat[t_ind,] <- out_temp
    }
  }

  return(out_mat)
}







#' Calculate absolute risk profiles after non-terminal event
#'
#' This function calculates absolute risk profiles conditional on non-terminal
#' event already occurring.
#'
#' @inheritParams calc_risk_logitbayes
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 4 element row matrix
#'   of probabilities. If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix of probabilities.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
# calc_risk_term <- function(para, Xmat1, Xmat2, Xmat3, XmatD, hazard,
#                            t_cutoff, tol=1e-3, frailty=TRUE,
#                            beta2frail=1, beta3frail=1, betaDfrail=1, #temporary, just get them in the mix
#                            type="marginal", gamma=1,model="semi-markov",
#                            h3_tv="none",h3tv_knots=NULL,
#                            logit_ind=TRUE,logit_tv="none",logit_tv_knots=NULL,
#                            n_quad = 15){
#   #notice reduced default tolerance
#   # browser()
#   ##TO START, EXTRACT POINT ESTIMATES OF ALL PARAMETERS FROM MODEL OBJECT##
#   ##*********************************************************************##
#
#   n <- max(1,nrow(Xmat1),nrow(Xmat2),nrow(Xmat3))
#   t_length <- length(t_cutoff)
#   #standardize namings for the use of "switch" below
#   #for now, restrict to weibull
#   stopifnot(tolower(hazard) %in% c("wb","weibull"))
#   hazard <- switch(tolower(hazard),
#                    wb="weibull",weibull="weibull")
#   stopifnot(tolower(type) %in% c("c","conditional","m","marginal"))
#   type <- switch(tolower(type),
#                  c="conditional",conditional="conditional",m="marginal",marginal="marginal")
#   #for now, restrict to semi-markov
#   stopifnot(tolower(model) %in% c("sm","semi-markov"))
#   model <- switch(tolower(model),
#                   sm="semi-markov","semi-markov"="semi-markov",m="markov",markov="markov")
#
#   #right now, we assume that h3 and logit have same knots and denote that by leaving logit_tv_knots blank,
#   #so throw an error if not
#   stopifnot(is.null(logit_tv_knots))
#
#
#   #set number of baseline parameters to be found in para vector
#   nP01 <- nP02 <- nP03 <- 2
#   #always assume theta parameter is present even for non-frailty model
#   nP0 <- nP01 + nP02 + nP03 + 1
#   if(frailty){
#     theta <- para[nP0] #for bayesian model, theta is not logged, so just store directly
#     if(type=="conditional" & length(gamma)==1){
#       gamma <- rep(gamma,n)
#     }
#   } else{ #if no frailty, can't be "marginalized" over frailties
#     theta <- 0
#     type <- "conditional"
#     gamma <- rep(1,n)
#   }
#
#   #now, work backwards to dig out all of the relevant parameters
#   #coefficients related to logit model
#   if(logit_ind){
#     #specify different forms by which t1 can be incorporated into the logit model
#     #assumes that effect of t1 is at the end of the set of betaD coefficients
#     if(tolower(logit_tv) == "linear"){
#       n_logit_tv <- 1
#       betaD_tv_linear <- utils::tail(para,n=n_logit_tv)
#     } else if(tolower(h3_tv) %in% c("pw","piecewise")){
#       if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
#       stopifnot(!is.null(logit_tv_knots))
#       if(logit_tv_knots[1] != 0){logit_tv_knots <- c(0,logit_tv_knots)}
#       if(utils::tail(logit_tv_knots, n=1) != Inf){logit_tv_knots <- c(logit_tv_knots,Inf)}
#       n_logit_tv <- length(logit_tv_knots) - 2
#       betaD_tv <- c(0,utils::tail(para,n=n_logit_tv))
#       betaD_tv_linear <- 0
#     } else{
#       n_logit_tv <- 0
#       betaD_tv_linear <- 0
#     }
#
#     #now, add in the component corresponding with the logit model
#     if(!is.null(XmatD) && !(ncol(XmatD)==0)){
#       nPD <- ncol(XmatD)
#       betaD <- util::tail(para,n=n_logit_tv + nPD)[1:nPD]
#       etaD <- as.vector(XmatD %*% betaD)
#     } else{
#       nPD <- 0
#       etaD <- 0
#     }
#
#     nP_logit <- 1+nPD+n_logit_tv
#     beta0D <- util::tail(para,n=nP_logit)[1]
#   } else{
#     nP_logit <- 0
#     etaD <- 0
#   }
#
#   #specify different forms by which t1 can be incorporated into h3
#   if(tolower(h3_tv) == "linear"){
#     if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
#     n_h3_tv <- 1
#     beta3_tv_linear <- utils::tail(para,n = n_h3_tv + nP_logit)[1:n_h3_tv]
#   } else if(tolower(h3_tv) %in% c("pw","piecewise")){
#     if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
#     stopifnot(!is.null(tv_knots))
#     if(tv_knots[1] != 0){tv_knots <- c(0,tv_knots)}
#     if(utils::tail(tv_knots, n=1) != Inf){tv_knots <- c(tv_knots,Inf)}
#     n_h3_tv <- length(tv_knots) - 2
#     beta3_tv_pw <- c(0,utils::tail(para,n=n_h3_tv + nP_logit)[1:n_h3_tv])
#     beta3_tv_linear <- 0
#   } else{
#     n_h3_tv <- 0
#     beta3_tv_linear <- 0
#   }
#
#   if(!is.null(Xmat3) && !(ncol(Xmat3)==0)){
#     nP3 <- ncol(Xmat3)
#     beta3 <- utils::tail(para,n = nP3 + n_h3_tv + nP_logit)[1:nP3]
#     eta3 <- as.vector(Xmat3 %*% beta3)
#   } else{
#     nP3 <- 0
#     eta3 <- 0
#   }
#
#   ##Set up the hazard functions##
#   ##***************************##
#   alpha3=para[6]
#   kappa3=para[5]
#
#   ##set up the components for numerically integrating over the frailties
#   gh_adj_nodes <- get_ghquad_pointsweights(n_quad=n_quad)$points * sqrt(2 * theta)
#   gh_weights <- get_ghquad_pointsweights(n_quad=n_quad)$weights
#
#   #now, in the below, define an additional version of each function we integrate over
#   #to compute `marginally`
#
#   if(logit_ind){
#     pi_vec <- function(t1,index,gamma,betaD_tv_const,betaD_tv_lin){
#       if(length(gamma) > 1){
#         return(stats::plogis(q = beta0D + etaD[index] + betaD_tv_const + betaD_tv_lin * t1 + betaDfrail * log(gamma[index])))
#       } else{
#         return(stats::plogis(q = beta0D + etaD[index] + betaD_tv_const + betaD_tv_lin * t1 + betaDfrail * log(gamma)))
#       }
#     }
#   } else{
#     pi_vec <- function(t1,index,gamma,betaD_tv_const,betaD_tv_lin){ return(0)}
#   }
#
#
#   ##******************************************##
#   ## Calculating posterior predictive density ##
#   ##******************************************##
#
#   #probability of surviving from to time t2 given that non-terminal event happened at t1
#   #naming convention comes from Putter (2007)
#   S_2r <- function(t2,t1,index){
#     if(tolower(h3_tv) %in% c("pw","piecewise")){
#       curr_interval <- findInterval(x = t1,vec = tv_knots,left.open = TRUE)
#       beta3_tv <- beta3_tv[curr_interval]
#     } else if(tolower(h3_tv) == "linear"){
#       beta3_tv <- beta3_tv_linear * t1
#     } else{
#       beta3_tv <- 0
#     }
#     if(logit_ind){
#       #assume they use the same intervals for now.
#       betaD_tv_const <- if(tolower(logit_tv) %in% c("pw","piecewise")) betaD_tv[curr_interval] else 0
#       betaD_tv_lin <- betaD_tv_linear * t1
#     } else{
#       betaD_tv_const <- 0
#       betaD_lin <- 0
#     }
#     if(type == "marginal"){
#       temp <- 0
#       for(x in 1:n_quad){
#         temp2 <- (1-pi_vec(t1 = t1,index = index,gamma = exp(gh_adj_nodes[x]),
#                  betaD_tv_const = betaD_tv_const,betaD_tv_lin = betaD_tv_lin)) *
#           exp( -exp(gh_adj_nodes[x] * beta3frail + as.vector(eta3[index]) + beta3_tv) *
#                  kappa3 * (t2 - t1)^alpha3)
#         temp <- temp + temp2 * gh_weights[x]
#       }
#       temp / sqrt(pi)
#     } else{
#       (1-pi_vec(t1 = t1,index = index,gamma = gamma,
#           betaD_tv_const = betaD_tv_const,betaD_tv_lin = betaD_tv_lin)) *
#         exp(-gamma[index] * exp(as.vector(eta3[index]) + beta3_tv) *
#             kappa3 * (t2 - t1)^alpha3)
#     }
#   }
#
#   ##ACTUALLY COMPUTING THE PROBABILITIES##
#   ##************************************##
#
#   if(n > 1){
#     if(t_length > 1){
#       if(t_start_length > 1){
#         out_mat <- array(dim=c(t_length,t_start_length,n),dimnames = list(paste0("t",t_cutoff),paste0("t",t_start,"_1"),paste0("i",1:n)))
#         for(i in 1:n){out_mat[,,i] <- outer(t_cutoff,t_start,function(x,y){S_2r(x,y,i)})}
#       } else{
#         out_mat <- outer(ids,t_cutoff,function(x,y){S_2r(y,t_start,x)})
#       }
#     } else{
#       if(t_start_length > 1){
#         out_mat <- outer(ids,t_start,function(x,y){S_2r(t_cutoff,y,x)})
#       }
#     }
#   } else{
#     out_mat <- outer(t_cutoff,t_start,function(x,y){S_2r(x,y,1)})
#   }
#
#
#   return(out_mat)
# }


