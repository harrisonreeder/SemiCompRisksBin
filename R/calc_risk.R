#' Calculate absolute risk profiles
#'
#' This function calculates absolute risk profiles under weibull baseline hazard specifications.
#'
#' @param para A numeric vector of parameters, arranged as follows:
#'   the first \eqn{k_1+k_2+k_3} elements correspond to the baseline hazard parameters,
#'   then the \eqn{k_1+k_2+k_3+1} element corresponds to the gamma frailty log-variance parameter,
#'   then the last\eqn{q_1+q_2+q_3} elements correspond with the regression parameters.
#' @param Xmat1,Xmat2,Xmat3 Numeric matrices with \eqn{n} rows and \eqn{q_1,q_2,q_3} columns containing covariates.
#' @param t_cutoff Numeric vector indicating the time(s) to compute the risk profile.
#' @param t_start Numeric scalar indicating the dynamic start time to compute the risk profile. Set to 0 by default.
#' @param tol Numeric value for the tolerance of the numerical integration procedure.
#' @param type String either indicating 'marginal' for population-averaged probabilities,
#'   or 'conditional' for probabilities computed at the specified gamma
#' @param gamma Numeric value indicating the fixed level of the frailty assumed for predicted probabilities,
#'   if 'type' is set to 'conditional'
#' @param h3_tv String indicating whether there is an effect of t1 on hazard 3.
#' @param tv_knots for piecewise effect of t1 in h3, these are the knots at which the effect jumps
#' @param hazard String specifying the form of the baseline hazard.
#' @param frailty Boolean indicating whether a gamma distributed subject-specific frailty should
#'   be included. Currently this must be set to TRUE.
#' @param model String specifying the transition assumption
#' @param knots_list Used for hazard specifications besides Weibull, a
#'   list of three increasing sequences of integers, each corresponding to
#'   the knots for the flexible model on the corresponding transition baseline hazard.
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 4 element row matrix
#'   of probabilities. If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix of probabilities.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
calc_risk <- function(para, Xmat1, Xmat2, Xmat3,hazard,knots_list=NULL,
                         t_cutoff, t_start=0, tol=1e-3, frailty=TRUE,
                         type="marginal", gamma=1,model="semi-markov",
                         h3_tv="none",tv_knots=NULL){
  #notice reduced default tolerance
  # browser()
  ##TO START, EXTRACT POINT ESTIMATES OF ALL PARAMETERS FROM MODEL OBJECT##
  ##*********************************************************************##

  n <- max(1,nrow(Xmat1),nrow(Xmat2),nrow(Xmat3))
  t_length <- length(t_cutoff)
  #standardize namings for the use of "switch" below
  stopifnot(tolower(hazard) %in% c("wb","weibull","pw","piecewise"))
  hazard <- switch(tolower(hazard),
    wb="weibull",weibull="weibull",pw="piecewise",piecewise="piecewise")
  stopifnot(tolower(type) %in% c("c","conditional","m","marginal"))
  type <- switch(tolower(type),
    c="conditional",conditional="conditional",m="marginal",marginal="marginal")
  stopifnot(tolower(model) %in% c("sm","semi-markov","m","markov"))
  model <- switch(tolower(model),
    sm="semi-markov","semi-markov"="semi-markov",m="markov",markov="markov")

  if(hazard == "weibull"){
    nP01 <- nP02 <- nP03 <- 2
  } else{
    stopifnot(!is.null(knots_list))
    #left pad with a zero if it is not already present
    if(knots_list[[1]][1] != 0){knots_list[[1]] <- c(0,knots_list[[1]])}
    if(knots_list[[2]][1] != 0){knots_list[[2]] <- c(0,knots_list[[2]])}
    if(knots_list[[3]][1] != 0){knots_list[[3]] <- c(0,knots_list[[3]])}
    nP01 <- length(knots_list[[1]])
    nP02 <- length(knots_list[[2]])
    nP03 <- length(knots_list[[3]])
  }

  if(frailty){
    nP0 <- nP01 + nP02 + nP03 + 1
    theta <- exp(para[nP0])
    if(type=="conditional" & length(gamma)==1){
      gamma <- rep(gamma,n)
    }
  } else{
    nP0 <- nP01 + nP02 +nP03
    type <- "conditional"
    gamma <- rep(1,n)
  }

  if(!is.null(Xmat1) & !(ncol(Xmat1)==0)){
    nP1 <- ncol(Xmat1)
    beta1 <- para[(1+nP0):(nP0+nP1)]
    eta1 <- as.vector(Xmat1 %*% beta1)
  } else{
    nP1 <- 0; eta1 <- 0
  }
  if(!is.null(Xmat2) & !(ncol(Xmat2)==0)){
    nP2 <- ncol(Xmat2)
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    eta2 <- as.vector(Xmat2 %*% beta2)
  } else{
    nP2 <- 0; eta2 <- 0
  }
  if(!is.null(Xmat3) & !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0; eta3 <- 0
  }

  #specify different forms by which t1 can be incorporated into h3
  if(tolower(h3_tv) == "linear"){
    if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
    n_tv <- 1
    beta3_tv_linear <- utils::tail(para,n = n_tv)
  } else if(tolower(h3_tv) %in% c("pw","piecewise")){
    if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
    stopifnot(!is.null(tv_knots))
    if(tv_knots[1] != 0){tv_knots <- c(0,tv_knots)}
    if(utils::tail(tv_knots, n=1) != Inf){tv_knots <- c(tv_knots,Inf)}
    n_tv <- length(tv_knots) - 2
    beta3_tv <- c(0,utils::tail(para,n=n_tv))
    beta3_tv_linear <- 0
  } else{
    n_tv <- 0
    beta3_tv_linear <- 0
  }

  #if the size of the parameter vector doesn't match the expected size, throw a fuss
  stopifnot(length(para) == nP0 + nP1 + nP2 + nP3 + n_tv)


  ##Set up the hazard functions##
  ##***************************##
  if(hazard == "weibull"){
    alpha1=exp(para[2])
    alpha2=exp(para[4])
    alpha3=exp(para[6])
    kappa1=exp(para[1])
    kappa2=exp(para[3])
    kappa3=exp(para[5])

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
  } else{
    phi1 <- as.numeric(para[(1):(nP01)])
    phi2 <- as.numeric(para[(1+nP01):(nP01+nP02)])
    phi3 <- as.numeric(para[(1+nP01+nP02):(nP01+nP02+nP03)])
    haz <- function(t,phi,knots){
      exp(phi)[findInterval(x=t, vec=knots, left.open=TRUE)]
    }
    Haz <- function(t,phi,knots){
      rowSums(sweep(x=pw_cum_mat(t,knots),MARGIN=2,STATS=exp(phi),FUN ="*"))
    }
  }

  ##******************************************##
  ## Calculating posterior predictive density ##
  ##******************************************##

  ##First, write functions that compute the integrand,
  ## which we feed into integration function
  ##***********************************************##

  #the univariate function when T1=infinity
  f_t2 <- switch(hazard,
    weibull=switch(type,
                   marginal=function(t2, index){
                     h2_const[index] * (t2)^alpha2_m1 *
                       (1 + theta*(H1_const[index] * (t2)^alpha1 +
                                     H2_const[index] * (t2)^alpha2) )^(-theta^(-1) - 1)
                     },
                   conditional=function(t2, index){
                     gamma[index] * h2_const[index] * (t2)^alpha2_m1 *
                       exp(-gamma[index]*(H1_const[index] * (t2)^alpha1 +
                                            H2_const[index] * (t2)^alpha2))
                     }),
    piecewise=switch(type,
                     marginal=function(t2, index){
                       haz(t=t2,phi=phi2,knots=knots_list[[2]]) * exp(eta2[index]) *
                         (1 + theta * (Haz(t=t2,phi=phi1,knots=knots_list[[1]])* exp(eta1[index]) +
                                         Haz(t=t2,phi=phi2,knots=knots_list[[2]])* exp(eta2[index])))^(-theta^(-1)-1)
                     },
                     conditional=function(t2, index){
                       gamma[index] * haz(t=t2,phi=phi2,knots=knots_list[[2]]) * exp(eta2[index]) *
                         exp(-gamma[index]*(Haz(t=t2,phi=phi1,knots=knots_list[[1]])* exp(eta1[index]) +
                                              Haz(t=t2,phi=phi2,knots=knots_list[[2]])* exp(eta2[index])))
                     })
  )

  #next, the different regions of the joint density on the upper triangle
  # here is the generic (semi-markov) formula if we pre-integrate t2 from u to v
  # \int \lambda_1(t_1)\exp(-[\Lambda_1(t_1)+\Lambda_2(t_1)]) \left[\exp(-\Lambda_3(u-t_1)) - \exp(-\Lambda_3(v-t_1))\right] dt_1
  # here is the generic (markov) formula if we pre-integrate t2 from u to v
  # \int \lambda_1(t_1)\exp(-[\Lambda_1(t_1)+\Lambda_2(t_1)-\Lambda_3(t_1)]) \left[\exp(-\Lambda_3(u)) - \exp(-\Lambda_3(v))\right] dt_1

  #function of t1 if we pre-integrate t2 from t1 to t_cutoff
  f_joint_t1_both <- function(time_pt1,t_cutoff,index,beta3_tv_const=0,beta3_tv_lin=0){
    H3_temp <- switch(hazard,
                      weibull=switch(model,
                        "semi-markov"= H3_const[index] * (t_cutoff - time_pt1)^alpha3 *
                          exp(beta3_tv_const + beta3_tv_lin * time_pt1),
                        "markov"=H3_const[index] * (t_cutoff^alpha3 - time_pt1^alpha3)),
                      piecewise=switch(model,
                        "semi-markov"=Haz(t=t_cutoff-time_pt1,phi=phi3,knots=knots_list[[3]]) *
                          exp(eta3[index] + beta3_tv_const + beta3_tv_lin*time_pt1),
                        "markov"=(Haz(t=t_cutoff,phi=phi3,knots=knots_list[[3]])-
                                    Haz(t=time_pt1,phi=phi3,knots=knots_list[[3]])) * exp(eta3[index])))
    #return the right value corresponding to the hazard and type
    switch(hazard,
           weibull=switch(type,
                          marginal={
                            h1_const[index] * time_pt1^alpha1_m1 * (
                              (1 + theta*(H1_const[index] * time_pt1^alpha1 +
                                            H2_const[index] * time_pt1^alpha2))^(-theta^(-1)-1) -
                                (1 + theta*(H1_const[index] * time_pt1^alpha1 +
                                              H2_const[index] * time_pt1^alpha2 +
                                              H3_temp))^(-theta^(-1)-1))
                          },
                          conditional={
                            gamma[index] * h1_const[index] * time_pt1^alpha1_m1 *
                              exp( -gamma[index]*(H1_const[index] * time_pt1^alpha1 + H2_const[index] * time_pt1^alpha2 ) ) *
                              ( 1 - exp( -gamma[index] * H3_temp))
                          }),
          piecewise=switch(type,
                           marginal={
                             haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) * (
                               (1 + theta*(Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                             Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index])))^(-theta^(-1)-1) -
                                 (1 + theta*(Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                               Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) +
                                               H3_temp))^(-theta^(-1)-1))
                           },
                           conditional={
                             gamma[index] * haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) *
                               exp(-gamma[index] * (Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                                      Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) )) *
                               ( 1 - exp(-gamma[index] * H3_temp))
                           }))
  }

  #function of t1 if we pre-integrate t2 from t_cutoff to infinity
  f_joint_t1_nonTerm <- function(time_pt1,t_cutoff,index,beta3_tv_const=0,beta3_tv_lin=0){
    H3_temp <- switch(hazard,
                      weibull=switch(model,
                       "semi-markov"= H3_const[index] * (t_cutoff - time_pt1)^alpha3 *
                         exp(beta3_tv_const + beta3_tv_lin * time_pt1),
                       "markov"=H3_const[index] * (t_cutoff^alpha3 - time_pt1^alpha3)),
                      piecewise=switch(model,
                       "semi-markov"=Haz(t=t_cutoff-time_pt1,phi=phi3,knots=knots_list[[3]]) *
                         exp(eta3[index] + beta3_tv_const + beta3_tv_lin*time_pt1),
                       "markov"=(Haz(t=t_cutoff,phi=phi3,knots=knots_list[[3]])-
                                   Haz(t=time_pt1,phi=phi3,knots=knots_list[[3]])) * exp(eta3[index])))
    #return the right value corresponding to the hazard and type
    switch(hazard,
           weibull=switch(type,
                          marginal={
                            h1_const[index] * time_pt1^alpha1_m1 *
                              (1 + theta * (H1_const[index] * time_pt1^alpha1 +
                                              H2_const[index] * time_pt1^alpha2 +
                                              H3_temp))^(-theta^(-1)-1)
                          },
                          conditional={
                            gamma[index] * h1_const[index] * time_pt1^alpha1_m1 *
                              exp(-gamma[index] * (H1_const[index] * time_pt1^alpha1 +
                                                     H2_const[index] * time_pt1^alpha2 +
                                                     H3_temp))
                          }),
           piecewise=switch(type,
                            marginal={
                              haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) *
                                (1 + theta * (Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                                Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) +
                                                H3_temp))^(-theta^(-1)-1)
                            },
                            conditional={
                              gamma[index] * haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) *
                                exp(-gamma[index] * (Haz(t=time_pt1,phi=phi1,knots=knots_list[[1]]) * exp(eta1[index]) +
                                                       Haz(t=time_pt1,phi=phi2,knots=knots_list[[2]])* exp(eta2[index]) +
                                                       H3_temp))
                            }))
  }


  ##finally, p_neither has a closed form, so we can write a function for it directly##
  #this derivation is actually identical to the "no event" likelihood contribution
  p_neither_func <- switch(hazard,
                           weibull=switch(type,
                              marginal=function(t2){
                                (1 + theta*(H1_const * (t2)^alpha1 +
                                              H2_const * (t2)^alpha2) )^(-theta^(-1))
                              },
                              conditional=function(t2){
                                exp(-gamma*(H1_const * (t2)^alpha1 +
                                                     H2_const * (t2)^alpha2))
                              }),
                           piecewise=switch(type,
                              marginal=function(t2){
                                (1 + theta * (Haz(t=t2,phi=phi1,knots=knots_list[[1]])* exp(eta1) +
                                                Haz(t=t2,phi=phi2,knots=knots_list[[2]])* exp(eta2)))^(-theta^(-1))
                              },
                              conditional=function(index,t2){
                                exp(-gamma*(Haz(t=t2,phi=phi1,knots=knots_list[[1]])* exp(eta1) +
                                                     Haz(t=t2,phi=phi2,knots=knots_list[[2]])* exp(eta2)))
                              }))

  ##ACTUALLY COMPUTING THE PROBABILITIES##
  ##************************************##

  #If we are computing 'dynamic probabilities' updated to some later timepoint t_start,
  #then we need to compute the probability of having experienced the event by t_start.
  #This is easier for 'neither' and 'tonly' outcomes because the inner integral does not depend on t_start.
  if(t_start > 0){
    p_neither_start <- p_neither_func(t2=t_start)
    p_tonly_start <- sapply(1:n,function(x){tryCatch(stats::integrate(f_t2, lower=0, upper=t_start, index=x)$value,
                                               error=function(cnd){return(NA)}) })
  } else{
    p_tonly_start <- 0
    p_neither_start <- 1
  }

  #this function allows inputs with multiple subjects, multiple time points, or both
  #therefore, we need to create the right data structure to contain the output.
  if(n > 1){
    if(t_length > 1){
      out_mat <- array(dim=c(t_length,4,n),dimnames = list(paste0("t",t_cutoff),c("p_ntonly","p_both","p_tonly","p_neither"),paste0("i",1:n)))
    } else{
      out_mat <- matrix(nrow=n,ncol=4,dimnames = list(paste0("i",1:n),c("p_ntonly","p_both","p_tonly","p_neither")))
    }
  } else{
      out_mat <- matrix(nrow=t_length,ncol=4,dimnames = list(paste0("t",t_cutoff),c("p_ntonly","p_both","p_tonly","p_neither")))
  }

  #loop through each time point, and compute the predicted probability at that time point for all subjects
  #each probability is an n-length vector
  for(t_ind in 1:t_length){
    t_temp <- t_cutoff[t_ind]
    p_tonly <- sapply(1:n,function(x){tryCatch(stats::integrate(f_t2, lower=0, upper=t_temp, index=x)$value,
                                               error=function(cnd){return(NA)}) })
    p_neither <- p_neither_func(t2=t_temp)

    #now, we have to be careful computing the probabilities that include h3 because it may depend on t1
    p_both_start <- p_ntonly_start <- rep(0,n)
    if(h3_tv %in% c("piecewise")){
      #a piecewise effect cannot be integrated in one go, because it is discontinuous
      #instead it must be divided into constant regions, integrated one region at a time and summed up
      curr_interval <- findInterval(x = t_temp,tv_knots,left.open = TRUE)
      #this is now a vector starting at 0, with elements at each change point up to the target time t_temp
      tv_knots_temp <- c(tv_knots[1:curr_interval],t_temp)
      if(t_temp == 0){ tv_knots_temp <- c(0,0)}

      #loop through the regions of constant effect, and add them together
      p_both <- p_ntonly <- rep(0,n)
      for(i in 1:(length(tv_knots_temp)-1)){
        p_both <- p_both + sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both, lower=tv_knots_temp[i], upper=tv_knots_temp[i+1],
                                                                   t_cutoff=t_temp, index=x,
                                                                   beta3_tv_const=beta3_tv[i], beta3_tv_lin=0)$value,
                                                  error=function(cnd){return(NA)}) })

        p_ntonly <- p_ntonly + sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_nonTerm, lower=tv_knots_temp[i], upper=tv_knots_temp[i+1],
                                                                     t_cutoff=t_temp, index=x,
                                                                     beta3_tv_const=beta3_tv[i], beta3_tv_lin=0)$value,
                                                    error=function(cnd){return(NA)}) })
      }
      if(t_start > 0){
        #need to similarly loop through the constant regions, now up to t_start instead of t_temp
        curr_interval <- findInterval(x = t_start,tv_knots,left.open = TRUE)
        tv_knots_temp <- c(tv_knots[1:curr_interval],t_start)
        for(i in 1:(length(tv_knots_temp)-1)){
          p_both_start <- p_both_start + sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both, lower=tv_knots_temp[i], upper=tv_knots_temp[i+1],
                                                                                          t_cutoff=t_temp, index=x,
                                                                                          beta3_tv_const=beta3_tv[i], beta3_tv_lin=0)$value,
                                                                         error=function(cnd){return(NA)}) })

          p_ntonly_start <- p_ntonly_start + sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_nonTerm, lower=tv_knots_temp[i], upper=tv_knots_temp[i+1],
                                                                                              t_cutoff=t_temp, index=x,
                                                                                              beta3_tv_const=beta3_tv[i], beta3_tv_lin=0)$value,
                                                                             error=function(cnd){return(NA)}) })
        }
      }
    } else{ #just normal, time-invariant prediction
      p_both <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both, lower=0, upper=t_temp,
                                                                          t_cutoff=t_temp, index=x,
                                                                   beta3_tv_const=0,beta3_tv_lin=beta3_tv_linear)$value,
                                                         error=function(cnd){return(NA)}) })
      p_ntonly <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_nonTerm, lower=0, upper=t_temp,
                                                                              t_cutoff=t_temp, index=x,
                                                                   beta3_tv_const=0,beta3_tv_lin=beta3_tv_linear)$value,
                                                             error=function(cnd){return(NA)}) })
      if(t_start > 0){
        p_both_start <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_both, lower=0, upper=t_start,
                                                                                        t_cutoff=t_temp, index=x,
                                                                   beta3_tv_const=0,beta3_tv_lin=beta3_tv_linear)$value,
                                                                       error=function(cnd){return(NA)}) })
        p_ntonly_start <- sapply(1:n,function(x){tryCatch(stats::integrate(f_joint_t1_nonTerm, lower=0, upper=t_start,
                                                                                            t_cutoff=t_temp, index=x,
                                                                   beta3_tv_const=0,beta3_tv_lin=beta3_tv_linear)$value,
                                                                           error=function(cnd){return(NA)}) })
      }
    }


    out_temp <- cbind(p_ntonly=(p_ntonly-p_ntonly_start)/p_neither_start,
                      p_both=(p_both-p_both_start)/p_neither_start,
                      p_tonly=(p_tonly-p_tonly_start)/p_neither_start,
                      p_neither=(p_neither)/p_neither_start)

    #I noticed that sometimes, if exactly one category has an NA, then we could back out the value
    #from the other categories. However, we don't want any to be negative.
    #I will also supply a warning if any of the rows are way off from 1.
    out_temp <- t(apply(out_temp,1,
                        function(x){
                          if(sum(is.na(x))==1){
                            x[is.na(x)]<- max(1-sum(x,na.rm=TRUE),0)
                          }
                          return(x)}))
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
#' @inheritParams calc_risk
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 4 element row matrix
#'   of probabilities. If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix of probabilities.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
calc_risk_term <- function(para, Xmat3,hazard,knots_list=NULL,
                      t_cutoff, t_start, tol=1e-3, frailty=TRUE,
                      type="marginal", gamma=1,model="semi-markov",
                      h3_tv="none",tv_knots=NULL){
  #notice reduced default tolerance
  # browser()
  ##TO START, EXTRACT POINT ESTIMATES OF ALL PARAMETERS FROM MODEL OBJECT##
  ##*********************************************************************##

  n <- nrow(Xmat3)
  t_length <- length(t_cutoff)
  t_start_length <- length(t_start)
  names(t_cutoff) <- paste0("t",t_cutoff)
  names(t_start) <- paste0("t",t_start,"_1")
  ids <- 1:n; names(ids) <- paste0("i",ids)
  #for now, assume that start is also t1
  # if(is.null(t1)){t_1 <- t_start}
  #standardize namings for the use of "switch" below
  stopifnot(tolower(hazard) %in% c("wb","weibull","pw","piecewise"))
  hazard <- switch(tolower(hazard),
                   wb="weibull",weibull="weibull",pw="piecewise",piecewise="piecewise")
  stopifnot(tolower(type) %in% c("c","conditional","m","marginal"))
  type <- switch(tolower(type),
                 c="conditional",conditional="conditional",m="marginal",marginal="marginal")
  stopifnot(tolower(model) %in% c("sm","semi-markov","m","markov"))
  model <- switch(tolower(model),
                  sm="semi-markov","semi-markov"="semi-markov",m="markov",markov="markov")

  if(hazard == "weibull"){
    nP01 <- nP02 <- nP03 <- 2
  } else{
    stopifnot(!is.null(knots_list))
    #left pad with a zero if it is not already present
    if(knots_list[[1]][1] != 0){knots_list[[1]] <- c(0,knots_list[[1]])}
    if(knots_list[[2]][1] != 0){knots_list[[2]] <- c(0,knots_list[[2]])}
    if(knots_list[[3]][1] != 0){knots_list[[3]] <- c(0,knots_list[[3]])}
    nP01 <- length(knots_list[[1]])
    nP02 <- length(knots_list[[2]])
    nP03 <- length(knots_list[[3]])
  }

  if(frailty){
    nP0 <- nP01 + nP02 + nP03 + 1
    theta <- exp(para[nP0])
    if(type=="conditional" & length(gamma)==1){
      gamma <- rep(gamma,n)
    }
  } else{
    nP0 <- nP01 + nP02 +nP03
    type <- "conditional"
    gamma <- rep(1,n)
  }

  #specify different forms by which t1 can be incorporated into h3
  if(tolower(h3_tv) == "linear"){
    if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
    n_tv <- 1
    beta3_tv_linear <- utils::tail(para,n = n_tv)
  } else if(tolower(h3_tv) %in% c("pw","piecewise")){
    if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
    stopifnot(!is.null(tv_knots))
    if(tv_knots[1] != 0){tv_knots <- c(0,tv_knots)}
    if(utils::tail(tv_knots, n=1) != Inf){tv_knots <- c(tv_knots,Inf)}
    n_tv <- length(tv_knots) - 2
    beta3_tv_pw <- c(0,utils::tail(para,n=n_tv))
    beta3_tv_linear <- 0
  } else{
    n_tv <- 0
    beta3_tv_linear <- 0
  }

  if(!is.null(Xmat3) & !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- utils::tail(para,n = nP3+n_tv)[1:nP3]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0
    eta3 <- 0
  }

  ##Set up the hazard functions##
  ##***************************##
  if(hazard == "weibull"){
    alpha3=exp(para[6])
    kappa3=exp(para[5])
  } else{
    phi3 <- as.numeric(para[(1+nP01+nP02):(nP01+nP02+nP03)])
    Haz <- function(t,phi,knots){
      rowSums(sweep(x=pw_cum_mat(t,knots),MARGIN=2,STATS=exp(phi),FUN ="*"))
    }
  }

  ##******************************************##
  ## Calculating posterior predictive density ##
  ##******************************************##

  #probability of surviving from to time t2 given that non-terminal event happened at t1
  #naming convention comes from Putter (2007)
  S_2r <- function(t2,t1,index){
    if(tolower(h3_tv) %in% c("pw","piecewise")){
      curr_interval <- findInterval(x = t1,vec = tv_knots,left.open = TRUE)
      beta3_tv <- beta3_tv_pw[curr_interval]
    } else if(tolower(h3_tv) == "linear"){
      beta3_tv <- beta3_tv_linear * t1
    } else{
      beta3_tv <- 0
    }
    H3_temp <- switch(hazard,
                      weibull=switch(model,
                                     "semi-markov"= kappa3 * exp(as.vector(eta3[index]) + beta3_tv) * (t2 - t1)^alpha3,
                                     "markov"=kappa3 * exp(as.vector(eta3[index])) * (t2^alpha3 - t1^alpha3)),
                      piecewise=switch(model,
                                       "semi-markov"=Haz(t=t2 - t1,phi=phi3,knots=knots_list[[3]]) *
                                         exp(eta3[index] + beta3_tv),
                                       "markov"=(Haz(t=t2,phi=phi3,knots=knots_list[[3]])-
                                                   Haz(t=t1,phi=phi3,knots=knots_list[[3]])) * exp(eta3[index])))
    #return the right value corresponding to the type
    switch(type,
           marginal=(1+theta*H3_temp)^(-theta^(-1)),
           conditional=exp(-gamma[index]*H3_temp))
  }

  ##ACTUALLY COMPUTING THE PROBABILITIES##
  ##************************************##

  if(n > 1){
    if(t_length > 1){
      if(t_start_length > 1){
        out_mat <- array(dim=c(t_length,t_start_length,n),dimnames = list(paste0("t",t_cutoff),paste0("t",t_start,"_1"),paste0("i",1:n)))
        for(i in 1:n){out_mat[,,i] <- outer(t_cutoff,t_start,function(x,y){S_2r(x,y,i)})}
      } else{
        out_mat <- outer(ids,t_cutoff,function(x,y){S_2r(y,t_start,x)})
      }
    } else{
      if(t_start_length > 1){
        out_mat <- outer(ids,t_start,function(x,y){S_2r(t_cutoff,y,x)})
      }
    }
  } else{
    out_mat <- outer(t_cutoff,t_start,function(x,y){S_2r(x,y,1)})
  }


  return(out_mat)
}


