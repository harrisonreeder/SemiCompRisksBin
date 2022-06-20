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
  prob_names <- names(quantile(0, probs = probs)) #matches how monitor function gets names
  #generate raw table with estimates
  if(object$setup$nChain>1){
    raw_mat <- t(apply(X = object$samples, MARGIN=3, FUN=quantile,probs=probs))
  } else {
    raw_mat <- t(apply(X = object$samples, MARGIN=2, FUN=quantile,probs=probs))
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
  } else if (object$class[2] == "IDlogit") {
    #theta estimates
    if(object$setup$frailty){
      tbl_theta <- raw_mat[7,prob_names,drop=FALSE]
      dimnames(tbl_theta) <- list("", c( "theta", "LL", "UL"))
    }

    #beta estimates
    beta.names <- unique(c(object$covnames$covNames1, object$covnames$covNames2,
                           object$covnames$covNames3, object$covnames$covNamesD))
    nP <- length(beta.names)
    output.coef <- NULL
    if(nP > 0){
      output <- matrix(NA, nrow=nP, ncol=12,
                       dimnames = list(beta.names, c("beta1", "LL", "UL",
                                                     "beta2", "LL", "UL",
                                                     "beta3", "LL", "UL",
                                                     "betaD", "LL", "UL")))
      p1 <- object$setup$nCov[1]; p2 <- object$setup$nCov[2]
      p3 <- object$setup$nCov[3]; pDnoint <- object$setup$nCov[4] - 1
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
    tbl_betaD0 <- raw_mat[8+p1+p2+p3,,drop=FALSE]
    dimnames(tbl_betaD0) <- list("",c("betaD0", "LL", "UL"))
  }

  value <- list(classFit=object$class, #psrf=psrf,
                coef=output.coef, h0=bh,
                setup=object$setup,conf.level=conf.level,
                if(object$class[2] == "IDlogit") theta=tbl_theta,
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


pred_helper_bayes <- function(tseq, xnew, pred_type,
                        conf.level, phi_array, beta_array, hazard_lab){

  #assume to start that everything is an array.

  # this little experiment shows how to get single vector of LPs from array
  # in order of chain1,1 - chain1,n ; chain2,1 - chain2,n ; chain3,1 - chain3,n
  # temp_array <- array(0,dim = c(2,3,4)) + 1:24
  # temp_vec <- c(1,1,1,1)
  # as.vector(apply(X = temp_array, MARGIN = 2, function(x) x %*% temp_vec))

  #for each "slice", do matrix multiply and then arrange into a long vector
  LP <- if(is.null(x1new)) 0 else as.vector(apply(X = beta_array, MARGIN = 2, function(x) x %*% xnew))

  kappa_vec <- as.vector(phi_array[,,1])
  alpha_vec <- as.vector(phi_array[,,2])

  haz_mat <- sapply(tseq, FUN = function(x) alpha_vec * kappa_vec * x^(alpha_vec - 1) * LP )
  Surv_mat <- sapply(tseq, FUN = function(x) alpha_vec * kappa_vec * x^(alpha_vec - 1) * LP )

  if(is.null(x$chain1$beta1.p)) {
    if(!is.null(x1new)) {
      stop("x1new does not have the same length as the x1 specified in the Formula.")
    }else {
      LP1 <- 0
    }
  } else if(length(x$chain1$beta1.p) != 0) {
    if(!is.null(x1new)) {
      #beta1
      p1 = dim(x$chain1$beta1.p)[2]
      beta1.p <- x$chain1$beta1.p
      if(nChain > 1) {
        for(i in 2:nChain) {
          nam <- paste("chain", i, sep="")
          beta1.p <- rbind(beta1.p, x[[nam]]$beta1.p)
        }
      }
      LP1 <- rowSums(beta1.p * matrix(x1new, nrow=dim(beta1.p)[1], ncol = dim(beta1.p)[2], byrow=T))
    } else {
      LP1 <- 0
    }
  }


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
    haz_out <- cbind(tseq,t(apply(X = haz_mat,MARGIN = 2,FUN = quantile, probs)))
    Surv_out <- cbind(tseq,t(apply(X = Surv_mat,MARGIN = 2,FUN = quantile, probs)))
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

  if(object$class[2] == "IDlogit"){
    #mark off the correct parameter vector start and end indices (for the hazard models)
    nP0_tot <- if (object$setup$frailty == TRUE) sum(nP0) + 1 else sum(nP0)
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
      haz_out <- cbind(tseq,t(apply(X = haz_mat,MARGIN = 2,FUN = quantile, probs)))
      Surv_out <- cbind(tseq,t(apply(X = Surv_mat,MARGIN = 2,FUN = quantile, probs)))
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

  # p1 <- object$setup$nCov[1]; p2 <- object$setup$nCov[2]
  # p3 <- object$setup$nCov[3]; pDnoint <- object$setup$nCov[4] - 1
  #
  # if (!is.null(xnew)) {
  #   stop("'xnew' is for univariate models so it must be specified as NULL for semi-competing risks models")
  # }
  # #
  # nP0_tot <- if (obj$frailty == TRUE) sum(nP0) + 1 else sum(nP0)
  # nP0_start <- 1 + c(0,nP0[1],nP0[1]+nP0[2])
  # nP0_end <- c(nP0[1],nP0[1]+nP0[2],nP0[1]+nP0[2]+nP0[3])
  # nP_start <- 1 + c(nP0_tot,nP0_tot+nP[1],nP0_tot+nP[1]+nP[2])
  # nP_end <- c(nP0_tot+nP[1],nP0_tot+nP[1]+nP[2],nP0_tot+nP[1]+nP[2]+nP[3])
  # x_list <- list(x1new,x2new,x3new)
  # value <- list()
  # for(type in c("h","S")){
  #   for(i in 1:3){
  #     value[[paste0(type,".",i)]] <-
  #       pred_helper(tseq = tseq,xnew = x_list[[i]],
  #                   pred_type = type,conf.level = conf.level,
  #                   phi = obj$estimate[nP0_start[i]:nP0_end[i]],
  #                   beta = obj$estimate[nP_start[i]:nP_end[i]],
  #                   hess_sub = obj$Finv[c(nP0_start[i]:nP0_end[i],
  #                                         if(!is.null(x_list[[i]])) nP_start[i]:nP_end[i]),
  #                                       c(nP0_start[i]:nP0_end[i],
  #                                         if(!is.null(x_list[[i]])) nP_start[i]:nP_end[i])],
  #                   hazard_lab = obj$class[4], knots_vec=obj$knots_list[[i]],
  #                   n_quad=n_quad,quad_method=quad_method)
  #     colnames(value[[paste0(type,".",i)]]) <- c("time",paste0(type,".",i),
  #                                                paste0("LL.",i),
  #                                                paste0("UL.",i))
  #   }
  # }
  #
  #
  #
  #
  # LP1
  #
  # #first, get LP1, LP2 and LP3
  # if(is.null(x$chain1$beta1.p)) {
  #   if(!is.null(x1new)) {
  #     stop("x1new does not have the same length as the x1 specified in the Formula.")
  #   }else {
  #     LP1 <- 0
  #   }
  # } else if(length(x$chain1$beta1.p) != 0) {
  #   if(!is.null(x1new)) {
  #     #beta1
  #     p1 = dim(x$chain1$beta1.p)[2]
  #     beta1.p <- x$chain1$beta1.p
  #     if(nChain > 1) {
  #       for(i in 2:nChain) {
  #         nam <- paste("chain", i, sep="")
  #         beta1.p <- rbind(beta1.p, x[[nam]]$beta1.p)
  #       }
  #     }
  #     LP1 <- rowSums(beta1.p * matrix(x1new, nrow=dim(beta1.p)[1], ncol = dim(beta1.p)[2], byrow=T))
  #   } else {
  #     LP1 <- 0
  #   }
  # }
  # if(is.null(x$chain1$beta2.p)) {
  #   if(!is.null(x2new)) {
  #     stop("x2new does not have the same length as the x2 specified in the Formula.")
  #   } else {
  #     LP2 <- 0
  #   }
  # } else if(length(x$chain1$beta2.p) != 0){
  #   if(!is.null(x2new)) {
  #     #beta2
  #     p2    = dim(x$chain1$beta2.p)[2]
  #     beta2.p <- x$chain1$beta2.p
  #     if(nChain > 1){
  #       for(i in 2:nChain){
  #         nam <- paste("chain", i, sep="")
  #         beta2.p <- rbind(beta2.p, x[[nam]]$beta2.p)
  #       }
  #     }
  #     LP2 <- rowSums(beta2.p * matrix(x2new, nrow=dim(beta2.p)[1], ncol = dim(beta2.p)[2], byrow=T))
  #   } else {
  #     LP2 <- 0
  #   }
  # }
  # if(is.null(x$chain1$beta3.p)) {
  #   if(!is.null(x3new)) {
  #     stop("x3new does not have the same length as the x3 specified in the Formula.")
  #   } else {
  #     LP3 <- 0
  #   }
  # } else if(length(x$chain1$beta3.p) != 0){
  #   if(!is.null(x3new)) {
  #     #beta3
  #     p3    = dim(x$chain1$beta3.p)[2]
  #     beta3.p <- x$chain1$beta3.p
  #     if(nChain > 1){
  #       for(i in 2:nChain){
  #         nam <- paste("chain", i, sep="")
  #         beta3.p <- rbind(beta3.p, x[[nam]]$beta3.p)
  #       }
  #     }
  #     LP3 <- rowSums(beta3.p * matrix(x3new, nrow=dim(beta3.p)[1], ncol = dim(beta3.p)[2], byrow=T))
  #   } else {
  #     LP3 <- 0
  #   }
  # }
  #
  # if(x$class[4] == "WB") {
  #   time1 <- time2 <- time3 <- seq(from=min(tseq), to=max(tseq), length=100)
  #   nStore <- length(x$chain1$alpha1.p)
  #   numSpl <- nStore * nChain
  #
  #   basehaz1 <- matrix(NA, numSpl, length(time1))
  #   basehaz2 <- matrix(NA, numSpl, length(time2))
  #   basehaz3 <- matrix(NA, numSpl, length(time3))
  #   basesurv1 <- matrix(NA, numSpl, length(time1))
  #   basesurv2 <- matrix(NA, numSpl, length(time2))
  #   basesurv3 <- matrix(NA, numSpl, length(time3))
  #
  #   alpha1.p    <- x$chain1$alpha1.p
  #   alpha2.p    <- x$chain1$alpha2.p
  #   alpha3.p    <- x$chain1$alpha3.p
  #   kappa1.p    <- x$chain1$kappa1.p
  #   kappa2.p    <- x$chain1$kappa2.p
  #   kappa3.p    <- x$chain1$kappa3.p
  #
  #   if(nChain > 1){
  #     for(i in 2:nChain){
  #       nam <- paste("chain", i, sep="")
  #       alpha1.p <- c(alpha1.p, x[[nam]]$alpha1.p)
  #       alpha2.p <- c(alpha2.p, x[[nam]]$alpha2.p)
  #       alpha3.p <- c(alpha3.p, x[[nam]]$alpha3.p)
  #
  #       kappa1.p <- c(kappa1.p, x[[nam]]$kappa1.p)
  #       kappa2.p <- c(kappa2.p, x[[nam]]$kappa2.p)
  #       kappa3.p <- c(kappa3.p, x[[nam]]$kappa3.p)
  #     }
  #   }
  #
  #   for(i in 1:numSpl){
  #     basehaz1[i, ] <- alpha1.p[i] * kappa1.p[i] * time1^(alpha1.p[i] - 1)
  #     basehaz2[i, ] <- alpha2.p[i] * kappa2.p[i] * time2^(alpha2.p[i] - 1)
  #     basehaz3[i, ] <- alpha3.p[i] * kappa3.p[i] * time3^(alpha3.p[i] - 1)
  #
  #     basesurv1[i, ] <- kappa1.p[i] * time1^(alpha1.p[i])
  #     basesurv2[i, ] <- kappa2.p[i] * time2^(alpha2.p[i])
  #     basesurv3[i, ] <- kappa3.p[i] * time3^(alpha3.p[i])
  #   }
  #
  #   basehaz1 <- basehaz1 * exp(matrix(LP1, nrow=dim(basehaz1)[1], ncol=dim(basehaz1)[2], byrow=F))
  #   basehaz2 <- basehaz2 * exp(matrix(LP2, nrow=dim(basehaz2)[1], ncol=dim(basehaz2)[2], byrow=F))
  #   basehaz3 <- basehaz3 * exp(matrix(LP3, nrow=dim(basehaz3)[1], ncol=dim(basehaz3)[2], byrow=F))
  #
  #   basesurv1 <- exp(-basesurv1 * exp(matrix(LP1, nrow=dim(basehaz1)[1], ncol=dim(basehaz1)[2], byrow=F)))
  #   basesurv2 <- exp(-basesurv2 * exp(matrix(LP2, nrow=dim(basehaz2)[1], ncol=dim(basehaz2)[2], byrow=F)))
  #   basesurv3 <- exp(-basesurv3 * exp(matrix(LP3, nrow=dim(basehaz3)[1], ncol=dim(basehaz3)[2], byrow=F)))
  #
  #   time1hz <- time1
  #   time2hz <- time2
  #   time3hz <- time3
  #
  #   if(tseq[1] == 0){
  #     time1hz <- time1[-1]
  #     time2hz <- time2[-1]
  #     time3hz <- time3[-1]
  #
  #     basehaz1 <- basehaz1[,-1]
  #     basehaz2 <- basehaz2[,-1]
  #     basehaz3 <- basehaz3[,-1]
  #   }
  #
  #   BH1Med <- apply(basehaz1, 2, median)
  #   BH1Ub <- apply(basehaz1, 2, quantile, prob = (1-conf.level/2))
  #   BH1Lb <- apply(basehaz1, 2, quantile, prob = conf.level/2)
  #   BH2Med <- apply(basehaz2, 2, median)
  #   BH2Ub <- apply(basehaz2, 2, quantile, prob = (1-conf.level/2))
  #   BH2Lb <- apply(basehaz2, 2, quantile, prob = conf.level/2)
  #   BH3Med <- apply(basehaz3, 2, median)
  #   BH3Ub <- apply(basehaz3, 2, quantile, prob = (1-conf.level/2))
  #   BH3Lb <- apply(basehaz3, 2, quantile, prob = conf.level/2)
  #
  #   BS1Med <- apply(basesurv1, 2, median)
  #   BS1Ub <- apply(basesurv1, 2, quantile, prob = (1-conf.level/2))
  #   BS1Lb <- apply(basesurv1, 2, quantile, prob = conf.level/2)
  #   BS2Med <- apply(basesurv2, 2, median)
  #   BS2Ub <- apply(basesurv2, 2, quantile, prob = (1-conf.level/2))
  #   BS2Lb <- apply(basesurv2, 2, quantile, prob = conf.level/2)
  #   BS3Med <- apply(basesurv3, 2, median)
  #   BS3Ub <- apply(basesurv3, 2, quantile, prob = (1-conf.level/2))
  #   BS3Lb <- apply(basesurv3, 2, quantile, prob = conf.level/2)
  # }
  #
  # BH1_tbl <- data.frame(time1=time1hz, h.1=BH1Med, LL.1=BH1Lb, UL.1=BH1Ub)
  # BH2_tbl <- data.frame(time2=time2hz, h.2=BH2Med, LL.2=BH2Lb, UL.2=BH2Ub)
  # BH3_tbl <- data.frame(time3=time3hz, h.3=BH3Med, LL.3=BH3Lb, UL.3=BH3Ub)
  #
  # BS1_tbl <- data.frame(time1=time1, S.1=BS1Med, LL.1=BS1Lb, UL.1=BS1Ub)
  # BS2_tbl <- data.frame(time2=time2, S.2=BS2Med, LL.2=BS2Lb, UL.2=BS2Ub)
  # BS3_tbl <- data.frame(time3=time3, S.3=BS3Med, LL.3=BS3Lb, UL.3=BS3Ub)
  #
  # value <- list(h.1=BH1_tbl, h.2=BH2_tbl, h.3=BH3_tbl, S.1=BS1_tbl, S.2=BS2_tbl, S.3=BS3_tbl)
  # value$xnew <- xnew
  # value$x1new <- x1new
  # value$x2new <- x2new
  # value$x3new <- x3new
  # value$tseq <- tseq
  # value$setup$model <- x$setup$model
  # value$class <- x$class
  #
  # class(value) <- "pred.Bayes_HReg"
  # return(value)
}

#' @export
predict_logit <- function(object, xnew=NULL, alpha = 0.05, ...) {
  # browser()
  conf.level = alpha
  yLim <- NULL
  nChain <- object$setup$nChain
  nP = object$setup$nCov
  nP0 = object$setup$nCov0
  value <- list()

  nP0_tot <- if (object$setup$frailty == TRUE) sum(nP0) + 1 else sum(nP0)
  nP_int <- 1 + nP0_tot+nP[1]+nP[2]+nP[3]
  nP_start <- 1 + nP0_tot+nP[1]+nP[2]+nP[3] + 1
  nP_end <- nP0_tot+nP[1]+nP[2]+nP[3] + nP[4]

  LP <- if(is.null(xnew)) 1 else as.vector(apply(X = object$samples[,,nP_start:nP_end], MARGIN = 2, function(x) x %*% as.vector(xnew)))

  probs <- c(0.5, alpha/2, 1-alpha/2)

  out_vec <- quantile(plogis(q=LP + as.vector(object$samples[,,nP_int])),probs=probs)
  names(out_vec) <- c("p","LL","UL")
  out_vec
}
