#' The function that simulates independent/cluster-correlated semi-competing
#'   risks data under semi-Markov Weibull models with logit submodel for immediate terminal event.
#'
#' @param x1,x2,x3,xD Covariate matrices with \code{n} rows.
#' @param beta1.true,beta2.true,beta3.true,betaD.true Vectors of true regression parameter values.
#'   The length of each vector should equal the number of columns in the corresponding covariate matrix.
#' @param beta2frail.true,beta3frail.true,betaDfrail.true scalar for the coefficient of the shared frailty in the logit submodel.
#' @param beta3tv.true,betaDtv.true Vectors of true regression parameter values for effects of T1 on the logit and h3 submodels.
#' @param alpha1.true,alpha2.true,alpha3.true,kappa1.true,kappa2.true,kappa3.true Vectors of true baseline parameter values.
#' @param theta.true True value for \eqn{\theta}.
#' @param anyD Boolean for whether to allow any "immediate" terminal events
#' @param h3tv_degree,Dtv_degree either the string "cs" indicating restricted cubic spline, or an integer for degree of time-varying hazard/odds ratio B-spline basis. (0 is piecewise constant)
#' @param frailty_type string denoting "gamma" for gamma-distributed frailty with variance \code{theta}, or "lognormal" for lognormal distributed frailty with log-frailty variance \code{theta}
#' @param cens A numeric vector of two elements. The right censoring times are generated from Uniform(\eqn{cens[1]}, \eqn{cens[2]}).
#'
#' @return returns a data.frame containing semi-competing risks outcomes from \code{n} subjects.
#'   It is of dimension \eqn{n\times 4}: the columns correspond to \eqn{y_1}, \eqn{\delta_1}, \eqn{y_2}, \eqn{\delta_2}. \cr
#'   \itemize{
#'   \item{y1}{a vector of \code{n} times to the non-terminal event}
#'   \item{y2}{a vector of \code{n} times to the terminal event}
#'   \item{delta1}{a vector of \code{n} censoring indicators for the non-terminal event time (1=event occurred, 0=censored)}
#'   \item{delta2}{a vector of \code{n} censoring indicators for the terminal event time (1=event occurred, 0=censored)}
#'   \item{deltaD}{a vector of \code{n} indicators whether the terminal event occurred immediately}
#'   }
#'
#' @export
simID_logit <- function(x1, x2, x3, xD=NULL,
                  beta1.true, beta2.true, beta3.true,
                  alpha1.true, alpha2.true, alpha3.true,
                  kappa1.true, kappa2.true, kappa3.true,
                  theta.true, frailty_type="gamma",
                  beta2frail.true=1, beta3frail.true=1, betaDfrail.true=1,
                  beta3tv.true=NULL, h3tv_degree=3,
                  anyD=TRUE, D_eps=0, D_dist="uniform",
                  betaD.true=NULL, betaDtv.true=NULL, Dtv_degree=3,
                  cens = c(0,0)) {
  browser()
  n <- dim(x1)[1]
  p1 <- dim(x1)[2]
  p2 <- dim(x2)[2]
  p3 <- dim(x3)[2]
  if(anyD & !is.null(xD)){
    pD <- dim(xD)[2]
    stopifnot(pD == length(betaD.true))
  } else{
    pD <- 0
  }

  if(theta.true >0) {
    if(tolower(frailty_type)=="gamma"){
      gamma.true <- stats::rgamma(n, 1/theta.true, 1/theta.true)
    } else {
      gamma.true <- exp(stats::rnorm(n,0,sqrt(theta.true)))
    }
  } else if(theta.true == 0){
    gamma.true <- rep(1, n)
  }

  LP1	<- if(p1>0) as.vector(beta1.true %*% t(x1)) else 0
  LP2	<- if(p2>0) as.vector(beta2.true %*% t(x2)) else 0
  LP3	<- if(p3>0) as.vector(beta3.true %*% t(x3)) else numeric(n) #made a vector bc it is subset automatically by delta1 below
  LPD	<- if(pD>0) as.vector(betaD.true %*% t(xD)) else 0

  T1 <- stats::rweibull(n, shape = alpha1.true, scale = exp(-(log(kappa1.true) +
                                                        LP1 + log(gamma.true))/alpha1.true))
  T2_temp <- stats::rweibull(n, shape = alpha2.true, scale = exp(-(log(kappa2.true) +
                                                        LP2 + beta2frail.true * log(gamma.true))/alpha2.true))
  yesT1 <- T1 < T2_temp
  if(cens[2] == 0){
    Cen <- rep(Inf,n)
  } else{
    Cen <- stats::runif(n, cens[1], cens[2])
  }

  #now, incorporate a possibly time-varying component into LP3
  if(!is.null(beta3tv.true)){ #if we set total number of parameters to 0, then we have no time-varying component.
    p3tv <- length(beta3tv.true)
    if(h3tv_degree == "linear"){ #linear
      stopifnot(p3tv==1)
      x3tv <- as.matrix(pmin(T1,T2_temp,Cen))
      colnames(x3tv) <- paste0("h3tv",1)
      LP3 <- LP3 + x3tv %*% beta3tv.true
      h3_knots <- c(0,Inf)

    } else if(h3tv_degree == "log1p") {
      stopifnot(p3tv==1)
      x3tv <- as.matrix(log1p(pmin(T1,T2_temp,Cen)))
      colnames(x3tv) <- paste0("h3tv",1)
      LP3 <- LP3 + x3tv %*% beta3tv.true
      h3_knots <- c(0,Inf)

    } else if(h3tv_degree == "cs"){ #cubic spline model
      #in cubic spline model, boundary knots are set directly at min/max endpoints,
      #so no need to fix at 0
      h3_quantile_seq <- seq(from = 0,to = 1, length.out = p3tv+1)
      #pmin(T1,Cen) gives vector of y1 times, from which we use just the observed non-terminal events
      #to define quantiles for the knots
      h3_knots <- stats::quantile(pmin(T1,Cen)[yesT1 & T1<Cen], h3_quantile_seq)
      x3tv <- splines::ns(x = pmin(T1,T2_temp,Cen), knots = h3_knots[-c(1,length(h3_knots))],
                          Boundary.knots = h3_knots[c(1,length(h3_knots))],
                          intercept = FALSE)
    } else { #if we don't use restricted cubic, then we are using a regular b-spline with specified degree
      h3tv_degree <- as.numeric(h3tv_degree)
      stopifnot(p3tv>=h3tv_degree)
      h3_quantile_seq <- seq(from = 0,to = 1, length.out = p3tv+2-h3tv_degree)[-c(1,p3tv+2-h3tv_degree)]
      #fixing piecewise endpoint at maximum is ok, because splines2 prediction will extrapolate beyond it
      h3_knots <- c(0,stats::quantile(pmin(T1,Cen)[yesT1 & T1<Cen],
                                      h3_quantile_seq),max(pmin(T1,T2_temp,Cen)))
      x3tv <- splines2::bSpline(x = pmin(T1,T2_temp,Cen), intercept = FALSE, degree = h3tv_degree,
                                knots = h3_knots[-c(1,length(h3_knots))],
                                Boundary.knots = h3_knots[c(1,length(h3_knots))])
    }
    colnames(x3tv) <- paste0("h3tv",1:p3tv)
    LP3 <- LP3 + x3tv %*% beta3tv.true
  } else{
    p3tv <- 0
  }

  #now, generate sojourn times for those with non-terminal event
  Sojourn <- rep(NA, n)
  Sojourn[yesT1] <- stats::rweibull(sum(yesT1), shape = alpha3.true,
                  scale = exp(-(log(kappa3.true) + LP3[yesT1] +
                                  beta3frail.true * log(gamma.true[yesT1]))/alpha3.true))

  T2 <- T2_temp
  T2[yesT1] <- T1[yesT1] + Sojourn[yesT1]

  y1 <- T1
  y2 <- T2
  delta1 <- delta2 <- rep(NA, n)

  #cases where terminal occurs before non-terminal and censoring
  ind01 <- which(T2 < T1 & T2 < Cen)
  y1[ind01] <- T2[ind01]
  delta1[ind01] <- 0
  delta2[ind01] <- 1

  #cases where nonterminal occurs, then censoring before terminal
  ind10 <- which(T1 < T2 & T1 < Cen & T2 >= Cen)
  y2[ind10] <- Cen[ind10]
  delta1[ind10] <- 1
  delta2[ind10] <- 0

  #cases where censoring occurs first
  ind00 <- which(T1 >= Cen & T2 >= Cen)
  y1[ind00] <- Cen[ind00]
  y2[ind00] <- Cen[ind00]
  delta1[ind00] <- 0
  delta2[ind00] <- 0

  #cases where nonterminal occurs, then terminal, then censoring
  ind11 <- which(T1 < Cen & T2 < Cen & T1 < T2)
  delta1[ind11] <- 1
  delta2[ind11] <- 1

  #for censored observations, replace the sojourn time with the censored time
  Sojourn[delta1==1 & delta2==0] <- (y2-y1)[delta1==1 & delta2==0]

  #this is a gut-check that the values I build the basis for h3tv on above
  #match the values of y1 for all observations that matter (e.g., those with delta1==1)
  stopifnot(all((y1==pmin(T1,Cen))[delta1==1]))

  if(anyD){
    LPD <- LPD + betaDfrail.true * log(gamma.true)
    if(!is.null(betaDtv.true)){ #if we set total number of parameters to 0, then we have no time-varying component.
      pDtv <- length(betaDtv.true)
      if(Dtv_degree == "linear"){ #linear
        stopifnot(pDtv==1)
        xDtv <- as.matrix(y1)
        colnames(xDtv) <- paste0("xDtv",1)
        LPD <- LPD + (xDtv %*% betaDtv.true)
        Dtv_knots <- c(0,Inf)

      } else if(Dtv_degree == "log1p") {
        stopifnot(pDtv==1)
        xDtv <- as.matrix(log1p(y1))
        colnames(xDtv) <- paste0("xDtv",1)
        LPD <- LPD + (xDtv %*% betaDtv.true)
        Dtv_knots <- c(0,Inf)

      } else if(Dtv_degree == "cs"){ #cubic spline model
        #in cubic spline model, boundary knots are set directly at min/max endpoints,
        #so no need to fix at 0
        Dtv_quantile_seq <- seq(from = 0,to = 1, length.out = pDtv+1)
        Dtv_knots <- stats::quantile(y1[delta1==1],Dtv_quantile_seq)
        xDtv <- splines::ns(x = y1, knots = Dtv_knots[-c(1,length(Dtv_knots))],
                            Boundary.knots = Dtv_knots[c(1,length(Dtv_knots))],
                            intercept = FALSE)
      } else { #if we don't use restricted cubic, then we are using a regular b-spline with specified degree
        Dtv_degree <- as.numeric(Dtv_degree)
        stopifnot(pDtv>=Dtv_degree)
        Dtv_quantile_seq <- seq(from = 0,to = 1, length.out = pDtv+2-Dtv_degree)[-c(1,pDtv+2-Dtv_degree)]
        #fixing piecewise endpoint at maximum is ok, because splines2 prediction will extrapolate beyond it
        Dtv_knots <- c(0,stats::quantile(y1[delta1==1],Dtv_quantile_seq),max(y1))
        xDtv <- splines2::bSpline(x = y1, intercept = FALSE, degree = Dtv_degree,
                                 knots = Dtv_knots[-c(1,length(Dtv_knots))],
                                 Boundary.knots = Dtv_knots[c(1,length(Dtv_knots))])
      }
      colnames(xDtv) <- paste0("xDtv",1:pDtv)
      LPD <- LPD + (xDtv %*% betaDtv.true)
    } else{
      pDtv <- 0
    }

    #this generates for all n obs, not just those with delta1==1, but that's ok.
    if(pD > 0 | pDtv > 0){
      probD.true <- stats::plogis(q = LPD)
      Dimmediate <- stats::rbinom(n,size = 1,prob = probD.true)
    } else{
      Dimmediate <- ProbD <- numeric(n)
    }

    delta1D <- as.numeric(delta1 & Dimmediate)
    delta2[delta1D==1] <- 0

    #for now, D_dist is only "uniform"
    #generate the "hypothetical" immediate event time uniformly on the interval (0,eps)
    D_eps.time <- stats::runif(n,min=0,max=D_eps)

    #if you had a non-immediate event, y2 is y1 + y_sm + D_eps (width of window)
    #note this shifts censored times as well as event times...
    y2[delta1==1 & Dimmediate==0] <- y2[delta1==1 & Dimmediate==0] + D_eps
    #if you had the immediate event, y2 is y1 + D_eps.time (i.e., a random little bit)
    y2[delta1D==1] <- y1[delta1D==1] + D_eps.time[delta1D==1]

    #this should correctly account for censoring because we've adjusted Sojourn above
    y_sm <- rep(0,n)
    y_sm[delta1==1 & Dimmediate==0] <- Sojourn[delta1==1 & Dimmediate==0]

  } else{
    pDtv <- 0
    y_sm <- rep(0,n)
    y_sm[delta1==1] <- Sojourn[delta1==1]
  }

  # #there are no differences between the spline used, and that used
  # all.equal(splines2:::predict.bSpline2(object = x3tv, newx = y1)[delta1==1,],
  #           x3tv[delta1==1,])
  # all.equal(splines2:::predict.bSpline2(object = x3tv, newx = y1), x3tv)
  # all.equal(splines:::predict.ns(object = x3tv, newx = y1)[delta1==1,],
  #           x3tv[delta1==1,])
  # all.equal(splines:::predict.ns(object = x3tv, newx = y1), x3tv)

  ret <- data.frame(cbind(y1, delta1, y2, delta2, y_sm,
                          deltaD = if(anyD) delta1D else 0,
                          D_eps, D_eps.time,
                          gamma.true, probD.true,
                          if(p3tv > 0) x3tv,
                          if(pDtv > 0) xDtv))

  if(!is.null(beta3tv.true)){
    attr(ret,which = "p3tv") <- p3tv
    attr(ret,which = "h3tv_degree") <- h3tv_degree
    attr(ret,which = "h3tv_knots") <- h3_knots
  }

  if(anyD & !is.null(betaDtv.true)){
    attr(ret,which = "pDtv") <- pDtv
    attr(ret,which = "Dtv_degree") <- Dtv_degree
    attr(ret,which = "Dtv_knots") <- Dtv_knots
  }

  return(ret)
}
