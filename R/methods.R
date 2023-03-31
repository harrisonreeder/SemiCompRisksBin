##methods

#' @export
print.Bayes_HReg2 <- function (x, digits = 3, alpha = 0.05, ...) {
  conf.level = alpha
  nChain = x$setup$nChain

  stopifnot(x$class[2] == "IDlogit" & x$class[3] == "Ind")
  #print header component
  cat("\nAnalysis of independent semi-competing risks data with immediate terminal event \n")
  cat("\nNumber of chains:    ", nChain,"\n")
  cat("Number of scans:     ", x$setup$numReps,"\n")
  cat("Thinning:            ", x$setup$thin,"\n")
  cat("Percentage of burnin: ", x$setup$burninPerc*100, "%\n", sep = "")

  #do the rest later
}

#' @export
summary.Bayes_HReg2 <- function(object, digits=3, alpha=0.05, ...) {
  # browser()
  conf.level <- alpha

  probs <- c(0.5, alpha/2, 1-alpha/2)
  prob_names <- names(stats::quantile(0, probs = probs)) #matches how monitor function gets names
  #generate raw table with estimates
  if(object$setup$nChain>1){
    raw_mat <- t(apply(X = object$samples, MARGIN=3, FUN=stats::quantile,probs=probs))
  } else {
    raw_mat <- t(apply(X = object$samples, MARGIN=2, FUN=stats::quantile,probs=probs))
  }
  # raw_mat <- as.matrix(rstan::monitor(object$samples, warmup = 0, print = FALSE,
  #                    probs = probs, digits_summary = digits))[,c(prob_names,"mean","sd","Rhat","Bulk_ESS"),drop=FALSE]

  if (object$class[2] == "Surv") {
    #baseline estimates
    bh <- log(raw_mat[1:2,])
    dimnames(bh) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"),
                         c("h-PM", "LL", "UL"))
    output.coef <- NULL
    if(object$setup$nCov>0){
      output.coef <- raw_mat[-(1:2),]
      colnames(output.coef) <- c("beta", "LL", "UL")
      rownames(output.coef) <- object$covnames
    }
  } else if (object$class[2] %in% c("ID", "IDlogit")) {
    #theta estimates
    if(object$setup$frailty){
      tbl_theta <- raw_mat[7,prob_names,drop=FALSE]
      dimnames(tbl_theta) <- list("", c( "theta", "LL", "UL"))
    }

    #beta estimates
    beta.names <- unique(c(object$covnames$covNames1, object$covnames$covNames2,
                   object$covnames$covNames3,
                   if(object$class[2]=="IDlogit") object$covnames$covNamesD))
    nP <- length(beta.names)
    output.coef <- NULL
    if(nP > 0){
      ncol_temp <- if(object$class[2]=="IDlogit") 12 else 9
      colnames_temp <- c(c("beta1", "LL", "UL",
                           "beta2", "LL", "UL",
                           "beta3", "LL", "UL"),
           if(object$class[2]=="IDlogit") c("betaD", "LL", "UL"))
      output <- matrix(NA, nrow=nP, ncol=ncol_temp,
                       dimnames = list(beta.names, colnames_temp))
      p1 <- object$setup$nCov[1]; p2 <- object$setup$nCov[2]
      p3 <- object$setup$nCov[3]
      pDnoint <- if(object$class[2]=="IDlogit") object$setup$nCov[4] - 1 else 0
      for(i in 1:nP) {
        if(p1 > 0){
          for(k in 1:p1) if(object$covnames$covNames1[k] == beta.names[i]) output[i,1:3] <- raw_mat[7 + k,prob_names]
        }
        if(p2 > 0){
          for(k in 1:p2) if(object$covnames$covNames2[k] == beta.names[i]) output[i,4:6] <- raw_mat[7 + p1 + k,prob_names]
        }
        if(p3 > 0){
          for(k in 1:p3) if(object$covnames$covNames3[k] == beta.names[i]) output[i,7:9] <- raw_mat[7 + p1 + p2 + k,prob_names]
        }
        if(pDnoint > 0){
          for(k in 1:pDnoint) if(object$covnames$covNamesD[k] == beta.names[i]) output[i,10:12] <- raw_mat[8 + p1 + p2 + p3 + k,prob_names]
        }
      }
      output.coef <- output #SemiCompRisks exponentiates this but I don't!!
    }

    #baseline estimates
    bh <- log(cbind(raw_mat[1:2,],raw_mat[3:4,],raw_mat[5:6,]))
    dimnames(bh) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"),
                         c("h1-PM", "LL", "UL", "h2-PM", "LL", "UL", "h3-PM", "LL", "UL"))
    if(object$class[2]=="IDlogit"){
      tbl_betaD0 <- raw_mat[8+p1+p2+p3,,drop=FALSE]
      dimnames(tbl_betaD0) <- list("",c("betaD0", "LL", "UL"))
    }
  }

  value <- list(classFit=object$class, #psrf=psrf,
                coef=output.coef, h0=bh,
                setup=object$setup,conf.level=conf.level,
                if(object$class[2] %in% c("ID", "IDlogit") & object$setup$frailty) theta=tbl_theta,
                if(object$class[2] == "IDlogit") betaD0=tbl_betaD0)
  class(value) <- "summ.Bayes_HReg2"

  return(value)

  # # convergence diagnostics
  #
  # # theta estimates
  # theta.p <- x$chain1$theta.p
  # if(nChain > 1){
  #   for(i in 2:nChain){
  #     nam <- paste("chain", i, sep="")
  #     theta.p <- rbind(theta.p, x[[nam]]$theta.p)
  #   }
  # }
  # theta.pMed <- apply(theta.p, 2, median)
  # theta.pUb <- apply(theta.p, 2, quantile, prob = (1-conf.level/2))
  # theta.pLb <- apply(theta.p, 2, quantile, prob = conf.level/2)
  # tbl_theta <- matrix(NA, 1, 3)
  # dimnames(tbl_theta) <- list("", c( "theta", "LL", "UL"))
  # tbl_theta[,1]	<- theta.pMed
  # tbl_theta[,2]	<- theta.pLb
  # tbl_theta[,3]	<- theta.pUb
  #
  # # beta estimates
  # beta.names <- unique(c(x$chain1$covNames1, x$chain1$covNames2, x$chain1$covNames3))
  # nP         <- length(beta.names)
  # output <- matrix(NA, nrow=nP, ncol=9)
  # dimnames(output) <- list(beta.names, c("exp(beta1)", "LL", "UL", "exp(beta2)", "LL", "UL", "exp(beta3)", "LL", "UL"))
  # if(length(x$chain1$beta1.p) != 0){
  #   #beta1
  #   p1	= dim(x$chain1$beta1.p)[2]
  #   beta.p <- x$chain1$beta1.p
  #   if(nChain > 1){
  #     for(i in 2:nChain){
  #       nam <- paste("chain", i, sep="")
  #       beta.p <- rbind(beta.p, x[[nam]]$beta1.p)
  #     }
  #   }
  #   beta.pMed <- apply(exp(beta.p), 2, median)
  #   beta.pSd <- apply(exp(beta.p), 2, sd)
  #   beta.pUb <- apply(exp(beta.p), 2, quantile, prob = (1-conf.level/2))
  #   beta.pLb <- apply(exp(beta.p), 2, quantile, prob = conf.level/2)
  #   tbl1 <- matrix(NA, p1, 3)
  #   rownames(tbl1) <- x$chain1$covNames1
  #   tbl1[,1]	<- beta.pMed
  #   tbl1[,2]	<- beta.pLb
  #   tbl1[,3]	<- beta.pUb
  #   for(i in 1:nP) {
  #     for(k in 1:p1) if(x$chain1$covNames1[k] == beta.names[i]) output[i,1:3] <- tbl1[k,]
  #   }
  # }
  # if(length(x$chain1$beta2.p) != 0){
  #   #beta2
  #   p2	= dim(x$chain1$beta2.p)[2]
  #   beta.p <- x$chain1$beta2.p
  #   if(nChain > 1){
  #     for(i in 2:nChain){
  #       nam <- paste("chain", i, sep="")
  #       beta.p <- rbind(beta.p, x[[nam]]$beta2.p)
  #     }
  #   }
  #   beta.pMed <- apply(exp(beta.p), 2, median)
  #   beta.pSd <- apply(exp(beta.p), 2, sd)
  #   beta.pUb <- apply(exp(beta.p), 2, quantile, prob = (1-conf.level/2))
  #   beta.pLb <- apply(exp(beta.p), 2, quantile, prob = conf.level/2)
  #   tbl2 <- matrix(NA, p2, 3)
  #   rownames(tbl2) <- x$chain1$covNames2
  #   tbl2[,1]	<- beta.pMed
  #   tbl2[,2]	<- beta.pLb
  #   tbl2[,3]	<- beta.pUb
  #   for(i in 1:nP) {
  #     for(k in 1:p2) if(x$chain1$covNames2[k] == beta.names[i]) output[i,4:6] <- tbl2[k,]
  #   }
  # }
  # if(length(x$chain1$beta3.p) != 0){
  #   #beta3
  #   p3	= dim(x$chain1$beta3.p)[2]
  #   beta.p <- x$chain1$beta3.p
  #   if(nChain > 1){
  #     for(i in 2:nChain){
  #       nam <- paste("chain", i, sep="")
  #       beta.p <- rbind(beta.p, x[[nam]]$beta3.p)
  #     }
  #   }
  #   beta.pMed <- apply(exp(beta.p), 2, median)
  #   beta.pSd <- apply(exp(beta.p), 2, sd)
  #   beta.pUb <- apply(exp(beta.p), 2, quantile, prob = (1-conf.level/2))
  #   beta.pLb <- apply(exp(beta.p), 2, quantile, prob = conf.level/2)
  #   tbl3 <- matrix(NA, p3, 3)
  #   rownames(tbl3) <- x$chain1$covNames3
  #   tbl3[,1]	<- beta.pMed
  #   tbl3[,2]	<- beta.pLb
  #   tbl3[,3]	<- beta.pUb
  #   for(i in 1:nP) {
  #     for(k in 1:p3) if(x$chain1$covNames3[k] == beta.names[i]) output[i,7:9] <- tbl3[k,]
  #   }
  # }
  #
  # output.coef <- NULL
  # if(nP > 0) {
  #   output.coef <- output
  # }
  #
  # # baseline estimates
  # alpha.p <- x$chain1$alpha1.p
  # if(nChain > 1){
  #   for(i in 2:nChain){
  #     nam <- paste("chain", i, sep="")
  #     alpha.p <- rbind(alpha.p, x[[nam]]$alpha1.p)
  #   }
  # }
  # alpha.pMed <- apply(log(alpha.p), 2, median)
  # alpha.pUb <- apply(log(alpha.p), 2, quantile, prob = (1-conf.level/2))
  # alpha.pLb <- apply(log(alpha.p), 2, quantile, prob = conf.level/2)
  # tbl_a1 <- c(alpha.pMed,alpha.pLb, alpha.pUb)
  # ##
  # alpha.p <- x$chain1$alpha2.p
  # if(nChain > 1){
  #   for(i in 2:nChain){
  #     nam <- paste("chain", i, sep="")
  #     alpha.p <- rbind(alpha.p, x[[nam]]$alpha2.p)
  #   }
  # }
  # alpha.pMed <- apply(log(alpha.p), 2, median)
  # alpha.pUb <- apply(log(alpha.p), 2, quantile, prob = (1-conf.level/2))
  # alpha.pLb <- apply(log(alpha.p), 2, quantile, prob = conf.level/2)
  # tbl_a2 <- c(alpha.pMed,alpha.pLb, alpha.pUb)
  # ##
  # alpha.p <- x$chain1$alpha3.p
  # if(nChain > 1){
  #   for(i in 2:nChain){
  #     nam <- paste("chain", i, sep="")
  #     alpha.p <- rbind(alpha.p, x[[nam]]$alpha3.p)
  #   }
  # }
  # alpha.pMed <- apply(log(alpha.p), 2, median)
  # alpha.pUb <- apply(log(alpha.p), 2, quantile, prob = (1-conf.level/2))
  # alpha.pLb <- apply(log(alpha.p), 2, quantile, prob = conf.level/2)
  # tbl_a3 <- c(alpha.pMed,alpha.pLb, alpha.pUb)
  # ##
  # kappa.p <- x$chain1$kappa1.p
  # if(nChain > 1){
  #   for(i in 2:nChain){
  #     nam <- paste("chain", i, sep="")
  #     kappa.p <- rbind(kappa.p, x[[nam]]$kappa1.p)
  #   }
  # }
  # kappa.pMed <- apply(log(kappa.p), 2, median)
  # kappa.pUb <- apply(log(kappa.p), 2, quantile, prob = (1-conf.level/2))
  # kappa.pLb <- apply(log(kappa.p), 2, quantile, prob = conf.level/2)
  # tbl_k1 <- c(kappa.pMed, kappa.pLb, kappa.pUb)
  # ##
  # kappa.p <- x$chain1$kappa2.p
  # if(nChain > 1){
  #   for(i in 2:nChain){
  #     nam <- paste("chain", i, sep="")
  #     kappa.p <- rbind(kappa.p, x[[nam]]$kappa2.p)
  #   }
  # }
  # kappa.pMed <- apply(log(kappa.p), 2, median)
  # kappa.pUb <- apply(log(kappa.p), 2, quantile, prob = (1-conf.level/2))
  # kappa.pLb <- apply(log(kappa.p), 2, quantile, prob = conf.level/2)
  # tbl_k2 <- c(kappa.pMed, kappa.pLb, kappa.pUb)
  # ##
  # kappa.p <- x$chain1$kappa3.p
  # if(nChain > 1){
  #   for(i in 2:nChain){
  #     nam <- paste("chain", i, sep="")
  #     kappa.p <- rbind(kappa.p, x[[nam]]$kappa3.p)
  #   }
  # }
  # kappa.pMed <- apply(log(kappa.p), 2, median)
  # kappa.pUb <- apply(log(kappa.p), 2, quantile, prob = (1-conf.level/2))
  # kappa.pLb <- apply(log(kappa.p), 2, quantile, prob = conf.level/2)
  # tbl_k3 <- c(kappa.pMed, kappa.pLb, kappa.pUb)
  # bh  <- matrix(c(tbl_k1, tbl_k2, tbl_k3, tbl_a1, tbl_a2, tbl_a3), 2, 9, byrow = T)
  # dimnames(bh) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"), c("h1-PM", "LL", "UL", "h2-PM", "LL", "UL", "h3-PM", "LL", "UL"))
  #
  # value <- list(classFit=x$class, psrf=psrf, theta=tbl_theta, coef=output.coef, h0=bh)
  # #final pieces
  # value$setup <- x$setup
  # value$conf.level <- conf.level
  #
  # class(value) <- "summ.Bayes_HReg"
  #
  # return(value)

}

#' @export
predict.Bayes_HReg2 <- function(object, xnew=NULL,
                                x1new=NULL, x2new=NULL, x3new=NULL,
                                tseq = seq(0,object$setup$ymax,length.out = 100),
                                alpha = 0.05, ...) {
  # browser()
  conf.level = alpha
  yLim <- NULL
  nChain <- object$setup$nChain
  nP = object$setup$nCov
  nP0 = object$setup$nCov0
  value <- list()

  if(object$class[2] == "Surv"){
    expLP <- if(is.null(xnew)) 1 else exp(as.vector(apply(X = object$samples[,,-(1:2)], MARGIN = 2, function(x) x %*% as.matrix(xnew))))
    kappa_vec <- as.vector(object$samples[,,1])
    alpha_vec <- as.vector(object$samples[,,2])
    probs <- c(0.5, conf.level/2, 1-conf.level/2)

    haz_mat <- apply(as.matrix(tseq), MARGIN=1, FUN = function(x) alpha_vec * kappa_vec * x^(alpha_vec - 1) * expLP )
    Surv_mat <- exp(-apply(as.matrix(tseq), MARGIN=1, FUN = function(x) kappa_vec * x^alpha_vec * expLP ))
    haz_out <- cbind(tseq,t(apply(X = haz_mat,MARGIN = 2,FUN = stats::quantile, probs)))
    Surv_out <- cbind(tseq,t(apply(X = Surv_mat,MARGIN = 2,FUN = stats::quantile, probs)))
    colnames(haz_out) <- c("time","h","LL","UL")
    colnames(Surv_out) <- c("time","S","LL","UL")
    if(tseq[1] == 0){
      value[["h"]] <- haz_out[-1,,drop=FALSE]
      value[["S"]] <- Surv_out[-1,,drop=FALSE]
    } else{
      value[["h"]] <- haz_out
      value[["S"]] <- Surv_out
    }
  }

  if(object$class[2] %in% c("ID","IDlogit")){
    #mark off the correct parameter vector start and end indices (for the hazard models)
    nP0_tot <- sum(nP0) + 1 #always leave space for frailty even for non-frailty models
    # nP0_tot <- if (object$setup$frailty == TRUE) sum(nP0) + 1 else sum(nP0)
    nP0_start <- 1 + c(0,nP0[1],nP0[1]+nP0[2])
    nP0_end <- c(nP0[1],nP0[1]+nP0[2],nP0[1]+nP0[2]+nP0[3])
    nP_start <- 1 + c(nP0_tot,nP0_tot+nP[1],nP0_tot+nP[1]+nP[2])
    nP_end <- c(nP0_tot+nP[1],nP0_tot+nP[1]+nP[2],nP0_tot+nP[1]+nP[2]+nP[3])
    x_list <- list(x1new,x2new,x3new)
    for(i in 1:3){

      expLP <- if(is.null(x_list[[i]])) 1 else exp(as.vector(apply(X = object$samples[,,nP_start[i]:nP_end[i]], MARGIN = 2, function(x) x %*% as.matrix(x_list[[i]]))))
      kappa_vec <- as.vector(object$samples[,,nP0_start[i]])
      alpha_vec <- as.vector(object$samples[,,nP0_start[i]+1])
      probs <- c(0.5, conf.level/2, 1-conf.level/2)

      haz_mat <- apply(as.matrix(tseq), MARGIN=1, FUN = function(x) alpha_vec * kappa_vec * x^(alpha_vec - 1) * expLP )
      Surv_mat <- exp(-apply(as.matrix(tseq), MARGIN=1, FUN = function(x) kappa_vec * x^alpha_vec * expLP ))
      haz_out <- cbind(tseq,t(apply(X = haz_mat,MARGIN = 2,FUN = stats::quantile, probs)))
      Surv_out <- cbind(tseq,t(apply(X = Surv_mat,MARGIN = 2,FUN = stats::quantile, probs)))
      colnames(haz_out) <- c("time",paste0("h.",i),paste0("LL.",i),paste0("UL.",i))
      colnames(Surv_out) <- c("time",paste0("S.",i),paste0("LL.",i),paste0("UL.",i))
      if(tseq[1] == 0){
        value[[paste0("h.",i)]] <- haz_out[-1,,drop=FALSE]
        value[[paste0("S.",i)]] <- Surv_out[-1,,drop=FALSE]
      } else{
        value[[paste0("h.",i)]] <- haz_out
        value[[paste0("S.",i)]] <- Surv_out
      }
    }
  }

  value$xnew <- xnew
  value$x1new <- x1new; value$x2new <- x2new; value$x3new <- x3new
  value$tseq <- if(tseq[1] == 0) tseq[-1] else tseq
  value$class <- object$class
  value$setup$model <- object$setup$model
  class(value) <- "pred.Bayes_HReg2"
  return(value)

}

#' function to predict probabilities from bayesian logistic regression
#'
#' function predicts probabilitites
#'
#' @param object object fit by Bayes_Logit
#' @param xnew new matrix of values at which to predict
#' @param alpha significance level for credible interval
#' @param ... extras
#'
#' @return a matrix
#' @export
predict_logit <- function(object, xnew=NULL, alpha = 0.05, ...) {
  # browser()
  conf.level = alpha
  yLim <- NULL
  nChain <- object$setup$nChain
  nP = object$setup$nCov
  nP0 = object$setup$nCov0
  value <- list()

  #always include "index" corresponding to theta, even for non-frailty models
  nP0_tot <- sum(nP0) + 1
  # nP0_tot <- if (object$setup$frailty == TRUE) sum(nP0) + 1 else sum(nP0)
  nP_int <- 1 + nP0_tot+nP[1]+nP[2]+nP[3]
  nP_start <- 1 + nP0_tot+nP[1]+nP[2]+nP[3] + 1
  nP_end <- nP0_tot+nP[1]+nP[2]+nP[3] + nP[4]

  LP <- if(is.null(xnew)) 1 else as.vector(apply(X = object$samples[,,nP_start:nP_end], MARGIN = 2, function(x) x %*% as.vector(xnew)))

  probs <- c(0.5, alpha/2, 1-alpha/2)

  out_vec <- stats::quantile(stats::plogis(q=LP + as.vector(object$samples[,,nP_int])),probs=probs)
  names(out_vec) <- c("p","LL","UL")
  out_vec
}


#' dumb little function to predict outcomes given non-terminal event timing.
#'
#' function predicts probabilitites
#'
#' @param object object fit by Bayes_Logit
#' @param x3new,xDnew new matrix of values at which to predict
#' @param alpha significance level for credible interval
#' @param tseq sequence of values at which to make predictions
#' @param logit_ind augmented model or standard model
#' @param type conditional or marginal
#' @param out format for output
#' @param n_quad number of quadrature points for prediction
#' @param ... extras
#'
#' @return a matrix
#' @export
predict_term <- function(object, x3new=NULL, xDnew=NULL,
                         logit_ind=TRUE, D_eps=0,
                         tseq, type = "conditional", out="long",
                         alpha = 0.05, n_quad=15, ...){
  # browser()
  conf.level = alpha
  probs <- c(0.5, conf.level/2, 1-conf.level/2)
  yLim <- NULL
  nChain <- object$setup$nChain
  nP = object$setup$nCov
  nP0 = object$setup$nCov0
  value <- list()

  if(tseq[1] == 0){ tseq <- tseq[-1]}
  t_len <- length(tseq)
  tseq_D <- pmax(0, tseq-D_eps)

  #basically, we're gonna compute two things: a vector of immediate event probs,
  #and a matrix with rows for each sample, and columns for each future time of
  #cumulative incidence conditional on not experiencing immediate event
  #we can use these two things to compute what we want.

  # logistic model

  #always include "index" corresponding to theta, even for non-frailty models
  nP0_tot <- sum(nP0) + 1
  # nP0_tot <- if (object$setup$frailty == TRUE) sum(nP0) + 1 else sum(nP0)
  # third hazard
  nP3_start <- 1 + nP0_tot+nP[1]+nP[2]
  nP3_end <- nP0_tot+nP[1]+nP[2]+nP[3]
  expLP3 <- if(is.null(x3new)) 1 else exp(as.vector(apply(
    X = object$samples[,,nP3_start:nP3_end], MARGIN = 2, function(x) x %*% as.vector(x3new))))

  if(logit_ind){
    nPD_int <- 1 + nP0_tot+nP[1]+nP[2]+nP[3]
    nPD_start <- 2 + nP0_tot+nP[1]+nP[2]+nP[3]
    nPD_end <- nP0_tot+nP[1]+nP[2]+nP[3] + nP[4]
    LPD <- if(is.null(xDnew)) 0 else as.vector(apply(
      X = object$samples[,,nPD_start:nPD_end], MARGIN = 2, function(x) x %*% as.vector(xDnew)))
    betaD0_vec <- as.vector(object$samples[,,nPD_int])
  }

  kappa_vec <- as.vector(object$samples[,,5])
  alpha_vec <- as.vector(object$samples[,,6])

  #start here and fix these
  if(type=="conditional"){
    S_cond_mat <- exp(-apply(as.matrix(tseq_D), MARGIN=1, FUN = function(x) kappa_vec * x^alpha_vec * expLP3 ))
    p_inst_vec <- if(logit_ind) stats::plogis(q=LPD + betaD0_vec) else 0
    S_mat <- (1-p_inst_vec) * S_cond_mat
    F_mat <- (1-p_inst_vec) * (1-S_cond_mat)
    CIF_mat <- p_inst_vec + (1-p_inst_vec) * (1-S_cond_mat)
  } else{
    S_cond_func <- function(t,kappa_vec,alpha_vec,expLP3,theta_vec,gh_node){
      exp(-kappa_vec * t^alpha_vec * expLP3 * exp(gh_node * sqrt(2 * theta_vec)))
    }
    p_cond_func <- function(betaD0_vec,LPD,theta_vec,gh_node){
      if(logit_ind) stats::plogis(q=LPD + betaD0_vec + gh_node * sqrt(2*theta_vec)) else 0
    }
    S_func <- function(t,kappa_vec,alpha_vec,expLP3,betaD0_vec,LPD,theta_vec,gh_node){
      (1-p_cond_func(betaD0_vec,LPD,theta_vec,gh_node)) *
        S_cond_func(t,kappa_vec,alpha_vec,expLP3,theta_vec,gh_node)
    }
    F_func <- function(t,kappa_vec,alpha_vec,expLP3,betaD0_vec,LPD,theta_vec,gh_node){
      (1-p_cond_func(betaD0_vec,LPD,theta_vec,gh_node)) *
        (1-S_cond_func(t,kappa_vec,alpha_vec,expLP3,theta_vec,gh_node))
    }
    CIF_func <- function(t,kappa_vec,alpha_vec,expLP3,betaD0_vec,LPD,theta_vec,gh_node){
      p_cond_func(betaD0_vec,LPD,theta_vec,gh_node) +
        (1-p_cond_func(betaD0_vec,LPD,theta_vec,gh_node)) *
        (1-S_cond_func(t,kappa_vec,alpha_vec,expLP3,theta_vec,gh_node))
    }

    theta_vec <- as.vector(object$samples[,,7])
    gh_nodes <- get_ghquad_pointsweights(n_quad=n_quad)$points
    gh_weights <- get_ghquad_pointsweights(n_quad=n_quad)$weights
    S_cond_mat <- S_mat <- F_mat <- CIF_mat <-
      matrix(data = 0, nrow=length(theta_vec), ncol=length(tseq_D))
    p_inst_vec <- numeric(length(theta_vec))
    for(x in 1:n_quad){
      S_cond_mat <- S_cond_mat + gh_weights[x] / sqrt(pi) *
        apply(as.matrix(tseq_D), MARGIN=1,
          FUN = S_cond_func, kappa_vec=kappa_vec,alpha_vec=alpha_vec,
                  expLP3=expLP3,theta_vec=theta_vec,gh_node=gh_nodes[x])
      p_inst_vec <- p_inst_vec + gh_weights[x] / sqrt(pi) *
        p_cond_func(betaD0_vec = betaD0_vec,LPD = LPD,theta_vec = theta_vec,gh_node = gh_nodes[x])
      S_mat <- S_mat + gh_weights[x] / sqrt(pi) *
        apply(as.matrix(tseq_D), MARGIN=1,
              FUN = S_func, kappa_vec=kappa_vec,alpha_vec=alpha_vec,
              expLP3=expLP3,betaD0_vec = betaD0_vec,LPD = LPD,
              theta_vec=theta_vec,gh_node=gh_nodes[x])
      F_mat <- F_mat + gh_weights[x] / sqrt(pi) *
        apply(as.matrix(tseq_D), MARGIN=1,
              FUN = F_func, kappa_vec=kappa_vec,alpha_vec=alpha_vec,
              expLP3=expLP3,betaD0_vec = betaD0_vec,LPD = LPD,
              theta_vec=theta_vec,gh_node=gh_nodes[x])
      CIF_mat <- CIF_mat + gh_weights[x] / sqrt(pi) *
        apply(as.matrix(tseq_D), MARGIN=1,
              FUN = CIF_func, kappa_vec=kappa_vec,alpha_vec=alpha_vec,
              expLP3=expLP3,betaD0_vec = betaD0_vec,LPD = LPD,
              theta_vec=theta_vec,gh_node=gh_nodes[x])
    }
  }

  #equivalent to what is reported in the baseline plots
  p_inst_out <- mean(p_inst_vec)
  S_cond_out <- apply(X = S_cond_mat,MARGIN = 2,FUN = mean)
  S_out <- apply(X = S_mat,MARGIN = 2,FUN = mean)
  F_out <- apply(X = F_mat,MARGIN = 2,FUN = mean)
  CIF_out <- apply(X = CIF_mat,MARGIN = 2,FUN = mean)

  #turned out that medians didn't add to 1, which was kind of a bummer!
  # p_inst_out <- stats::quantile(p_inst_vec,mean)
  # S_cond_out <- t(apply(X = S_cond_mat,MARGIN = 2,FUN = mean, probs))
  # S_out <- t(apply(X = S_mat,MARGIN = 2,FUN = stats::quantile, probs))
  # F_out <- t(apply(X = F_mat,MARGIN = 2,FUN = stats::quantile, probs))
  # CIF_out <- t(apply(X = CIF_mat,MARGIN = 2,FUN = stats::quantile, probs))
  if(out=="simple"){
    out_mat <- cbind(p_inst_out,F_out,S_out)
    # out_mat <- cbind(p_inst_out[1],F_out[,1],S_out[,1])
    colnames(out_mat) <- c("p_both_inst","p_both_noinst","p_ntonly")
    out_mat
  } else{
    # out_mat <- cbind(p_inst=p_inst_out,F=F_out,S=S_out,S_cond=S_cond_out,CIF=CIF_out)

    # names(p_inst_out) <- c("p_inst","p_inst_LL","p_inst_UL")
    # colnames(S_cond_out) <- c("S_cond","S_cond_LL","S_cond_UL")
    # colnames(S_out) <- c("S","S_LL","S_UL")
    # colnames(F_out) <- c("F","F_LL","F_UL")
    # colnames(CIF_out) <- c("CIF","CIF_LL","CIF_UL")
    list(tseq=tseq,tseq_D=tseq_D,p_inst=p_inst_out,
         S_cond=S_cond_out,S=S_out,F=F_out,CIF=CIF_out)
  }

}

