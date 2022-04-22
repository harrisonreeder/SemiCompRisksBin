#' The function that simulates independent/cluster-correlated semi-competing
#'   risks data under semi-Markov Weibull models with logit submodel for immediate terminal event.
#'
#' @param x1,x2,x3,xD Covariate matrices with \code{n} rows.
#' @param beta1.true,beta2.true,beta3.true,beta3tv.true,betaD.true,betaDtv.true Vectors of true regression parameter values.
#'   The length of each vector should equal the number of columns in the corresponding covariate matrix.
#' @param betaDfrail.true scalar for the coefficient of the shared frailty in the logit submodel.
#' @param beta3tv.true,betaDtv.true Vectors of true regression parameter values for effects of T1 on the logit and h3 submodels.
#' @param alpha1.true,alpha2.true,alpha3.true,kappa1.true,kappa2.true,kappa3.true Vectors of true baseline parameter values.
#' @param theta.true True value for \eqn{\theta}.
#' @param anyD Boolean for whether to allow any "immediate" terminal events
#' @param h3tv_degree,Dtv_degree either the string "cs" indicating restricted cubic spline, or an integer for degree of time-varying hazard/odds ratio basis. (0 is piecewise constant)
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
                  beta1.true, beta2.true, beta3.true, beta3tv.true=NULL,
                  betaD.true=NULL, betaDfrail.true=1, betaDtv.true=NULL,
                  alpha1.true, alpha2.true, alpha3.true,
                  kappa1.true, kappa2.true, kappa3.true,
                  theta.true, h3tv_degree=3,
                  anyD=TRUE, Dtv_degree=3,
                  frailty_type="gamma", cens) {
  # browser()
  n <- dim(x1)[1]
  p1 <- dim(x1)[2]
  p2 <- dim(x2)[2]
  p3 <- dim(x3)[2]
  if(!is.null(xD)){
    pD <- dim(xD)[2]
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
  LP3	<- if(p3>0) as.vector(beta3.true %*% t(x3)) else 0
  LPD	<- if(pD>0) as.vector(betaD.true %*% t(xD)) else 0

  Rind <- NULL
  R <- stats::rweibull(n, shape = alpha1.true, scale = exp(-(log(kappa1.true) +
                                                        LP1 + log(gamma.true))/alpha1.true))
  D <- stats::rweibull(n, shape = alpha2.true, scale = exp(-(log(kappa2.true) +
                                                        LP2 + log(gamma.true))/alpha2.true))
  yesR <- R < D
  Cen <- stats::runif(n, cens[1], cens[2])


  #now, incorporate a possibly time-varying component into LP3
  if(!is.null(beta3tv.true)){ #if we set total number of parameters to 0, then we have no time-varying component.
    p3tv <- length(beta3tv.true)
    if(h3tv_degree == "cs"){ #cubic spline model
      #in cubic spline model, boundary knots are set directly at min/max endpoints,
      #so no need to fix at 0
      h3_quantile_seq <- seq(from = 0,to = 1, length.out = p3tv+1)
      h3_knots <- stats::quantile(pmin(R,Cen)[yesR==1 & R<Cen], h3_quantile_seq)
      x3tv <- splines::ns(x = pmin(R,D,Cen), knots = h3_knots[-c(1,length(h3_knots))],
                          Boundary.knots = h3_knots[c(1,length(h3_knots))],
                          intercept = FALSE)
    } else { #if we don't use restricted cubic, then we are using a regular b-spline with specified degree
      stopifnot(p3tv>=h3tv_degree)
      h3_quantile_seq <- seq(from = 0,to = 1, length.out = p3tv+2-h3tv_degree)[-c(1,p3tv+2-h3tv_degree)]
      #fixing piecewise endpoint at maximum is ok, because splines2 prediction will extrapolate beyond it
      h3_knots <- c(0,stats::quantile(pmin(R,Cen)[yesR==1 & R<Cen],
                                      h3_quantile_seq),max(pmin(R,D,Cen)))
      x3tv <- splines2::bSpline(x = pmin(R,D,Cen), intercept = FALSE, degree = h3tv_degree,
                                knots = h3_knots[-c(1,length(h3_knots))],
                                Boundary.knots = h3_knots[c(1,length(h3_knots))])
    }
    colnames(x3tv) <- paste0("h3tv",1:p3tv)
    LP3 <- LP3 + x3tv %*% beta3tv.true
  } else{
    p3tv <- 0
  }


  D[yesR] <- R[yesR] + stats::rweibull(sum(yesR), shape = alpha3.true,
                                scale = exp(-(log(kappa3.true) + LP3[yesR] + log(gamma.true[yesR]))/alpha3.true))
  delta1 <- rep(NA, n)
  delta2 <- rep(NA, n)
  y1 <- R
  y2 <- D

  #cases where terminal occurs before non-terminal and censoring
  ind01 <- which(D < R & D < Cen)
  y1[ind01] <- D[ind01]
  delta1[ind01] <- 0
  delta2[ind01] <- 1

  #cases where nonterminal occurs, then censoring before terminal
  ind10 <- which(R < D & R < Cen & D >= Cen)
  y2[ind10] <- Cen[ind10]
  delta1[ind10] <- 1
  delta2[ind10] <- 0

  #cases where censoring occurs first
  ind00 <- which(R >= Cen & D >= Cen)
  y1[ind00] <- Cen[ind00]
  y2[ind00] <- Cen[ind00]
  delta1[ind00] <- 0
  delta2[ind00] <- 0

  #cases where nonterminal occurs, then terminal, then censoring
  ind11 <- which(R < Cen & D < Cen & R < D)
  delta1[ind11] <- 1
  delta2[ind11] <- 1

  #this is a gut-check that the values I build the basis for h3tv on above
  #match the values of y1 for all observations that matter (e.g., those with delta1==1)
  stopifnot(all((y1==pmin(R,Cen))[delta1==1]))

  if(anyD){
    LPD <- LPD + betaDfrail.true * log(gamma.true)
    if(!is.null(betaDtv.true)){ #if we set total number of parameters to 0, then we have no time-varying component.
      pDtv <- length(betaDtv.true)
      if(Dtv_degree == "cs"){ #cubic spline model
        #in cubic spline model, boundary knots are set directly at min/max endpoints,
        #so no need to fix at 0
        Dtv_quantile_seq <- seq(from = 0,to = 1, length.out = pDtv+1)
        Dtv_knots <- stats::quantile(y1[delta1==1],Dtv_quantile_seq)
        xDtv <- splines::ns(x = y1, knots = Dtv_knots[-c(1,length(Dtv_knots))],
                            Boundary.knots = Dtv_knots[c(1,length(Dtv_knots))],
                            intercept = FALSE)
      } else { #if we don't use restricted cubic, then we are using a regular b-spline with specified degree
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

    Dimmediate <- if(pD > 0 | pDtv > 0) stats::rbinom(n,size = 1,
                                    prob = stats::plogis(q = LPD)) else numeric(n)
    ind1D <- (delta1 & Dimmediate)
    y2[ind1D] <- y1[ind1D]
    delta2[ind1D] <- 0
    deltaD <- as.numeric(ind1D)
  }

  # #there are no differences between the spline used, and that used
  # all.equal(splines2:::predict.bSpline2(object = x3tv, newx = y1)[delta1==1,],
  #           x3tv[delta1==1,])
  # all.equal(splines2:::predict.bSpline2(object = x3tv, newx = y1), x3tv)
  # all.equal(splines:::predict.ns(object = x3tv, newx = y1)[delta1==1,],
  #           x3tv[delta1==1,])
  # all.equal(splines:::predict.ns(object = x3tv, newx = y1), x3tv)

  ret <- data.frame(cbind(y1, delta1, y2, delta2,
                          deltaD = if(anyD) deltaD else 0,
                          gamma.true,
                          if(p3tv > 0) x3tv,
                          if(pDtv > 0) xDtv))

  if(!is.null(beta3tv.true)){
    attr(ret,which = "p3tv") <- p3tv
    attr(ret,which = "h3tv_degree") <- h3tv_degree
    attr(ret,which = "h3tv_knots") <- h3_knots
  }

  if(!is.null(betaDtv.true)){
    attr(ret,which = "pDtv") <- pDtv
    attr(ret,which = "Dtv_degree") <- Dtv_degree
    attr(ret,which = "Dtv_knots") <- Dtv_knots
  }

  return(ret)
}
