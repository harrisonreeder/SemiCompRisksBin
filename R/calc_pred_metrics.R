

#' Get matrix of observed outcome categories
#'
#' This function returns a matrix giving the observed outcome categories of each observation at various
#'   time cutoffs.
#'
#' @param y1,y2 Numeric vectors of length \eqn{n} with (possibly censored) non-terminal and terminal event times
#' @param delta1,delta2 Numeric vectors of length \eqn{n}  with indicators of 1 if the event was observed and 0 otherwise
#' @inheritParams calc_risk
#'
#' @return a matrix or array.
#' @export
get_outcome_mat <- function(y1, y2, delta1, delta2, t_cutoff){

  n <- length(y1)
  t_length <- length(t_cutoff)
  if(n > 1){
    if(t_length > 1){
      out_mat <- array(dim=c(t_length,4,n),dimnames = list(paste0("t",t_cutoff),c("ntonly","both","tonly","neither"),paste0("i",1:n)))
    } else{
      out_mat <- matrix(nrow=n,ncol=4,dimnames = list(paste0("i",1:n),c("ntonly","both","tonly","neither")))
    }
  } else{
    out_mat <- matrix(nrow=t_length,ncol=4,dimnames = list(paste0("t",t_cutoff),c("ntonly","both","tonly","neither")))
  }

  for(t_ind in 1:t_length){
    #For cases where y=t_cutoff, I consider events that happened exactly at t_cutoff in categorization.
    neither <- t_cutoff[t_ind] < y1 | #neither
      y2 <= t_cutoff[t_ind] & delta1 == 0 & delta2 == 0 #neither
    ntonly <- y1 <= t_cutoff[t_ind] & t_cutoff[t_ind] < y2 | #ntonly
      y2 <= t_cutoff[t_ind] & delta1 == 1 & delta2 == 0 #ntonly
    tonly <- y2 <= t_cutoff[t_ind] & delta1 == 0 & delta2 == 1 #tonly
    both <- y2 <= t_cutoff[t_ind] & delta1 == 1 & delta2 == 1 #both

    out_temp <- cbind(ntonly=as.numeric(ntonly),
                      both=as.numeric(both),
                      tonly=as.numeric(tonly),
                      neither=as.numeric(neither))

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




#' Get matrix of observed outcome categories
#'
#' This function returns a matrix giving the observed outcome categories of each observation at various
#'   time cutoffs.
#'
#' @inheritParams get_outcome_mat
#' @inheritParams calc_risk
#' @param deltaD indicator that is 1 if subject has delta1==1 and observed immediate event, 0 else
#'
#' @return a matrix or array.
#' @export
get_outcome_mat_logit <- function(y1, y2, delta1, delta2, deltaD, t_cutoff){

  # browser()
  n <- length(y1)
  t_length <- length(t_cutoff)
  if(n > 1){
    if(t_length > 1){
      out_mat <- array(dim=c(t_length,5,n),dimnames = list(paste0("t",t_cutoff),c("ntonly","both_noinst","both_inst","tonly","neither"),paste0("i",1:n)))
    } else{
      out_mat <- matrix(nrow=n,ncol=5,dimnames = list(paste0("i",1:n),c("ntonly","both_noinst","both_inst","tonly","neither")))
    }
  } else{
    out_mat <- matrix(nrow=t_length,ncol=5,dimnames = list(paste0("t",t_cutoff),c("ntonly","both_noinst","both_inst","tonly","neither")))
  }

  for(t_ind in 1:t_length){
    #For cases where y=t_cutoff, I consider events that happened exactly at t_cutoff in categorization.
    neither <- t_cutoff[t_ind] < y1 | #neither
      y2 <= t_cutoff[t_ind] & delta1 == 0 & delta2 == 0 #neither
    ntonly <- y1 <= t_cutoff[t_ind] & t_cutoff[t_ind] < y2 | #ntonly
      y2 <= t_cutoff[t_ind] & delta1 == 1 & deltaD == 0 & delta2 == 0 #ntonly
    tonly <- y2 <= t_cutoff[t_ind] & delta1 == 0 & delta2 == 1 #tonly
    both_noinst <- y2 <= t_cutoff[t_ind] & delta1 == 1 &  deltaD == 0 & delta2 == 1 #both
    both_inst <- y2 <= t_cutoff[t_ind] & delta1 == 1 & deltaD == 1 #both

    out_temp <- cbind(ntonly=as.numeric(ntonly),
                      both_noinst=as.numeric(both_noinst),
                      both_inst=as.numeric(both_inst),
                      tonly=as.numeric(tonly),
                      neither=as.numeric(neither))

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


#' Get inverse probability of censoring weights
#'
#' This function returns a vector of inverse probability of censoring weights from an unadjusted Cox model
#'   for censoring times.
#'
#' @inheritParams get_outcome_mat
#' @inheritParams calc_risk
#'
#' @return a vector.
#' @export
get_ipcw_mat <- function(y2,delta2,t_cutoff){

  # browser()
  n <- length(y2)
  t_length <- length(t_cutoff)

  #this is Ghat, a non-parametric model of the 'survival distribution' of censoring var C.
  sfcens <- survival::survfit(survival::Surv(y2, delta2==0) ~ 1)

  #* want to get Ghat(z) where (following Graf 1999 'three categories')
  #* Category 1:
  #* z=y2- if y2<=s and delta2=1,
  #* Category 2:
  #* z=s if s<y2
  #* Category 3:
  #* z=Inf if y2<s and delta2=0, (aka, 1/Ghat(z)=0, I know they're subtly different)
  #* z=s if s=y2 and delta2=0, #this one situation is what I'm unsure about
  #* because graf and spitoni would say that if s=y2 and delta2=0, then z=Inf and result should be tossed.
  #* below I defer to them, but I'm just noting that to me there is logic in the other direction.
  #*
  #* so, we define
  #* z = min(y2,s)
  #* then change z=y2- if y2<=s and delta2=1
  #* then change z=Inf if y2<=s and delta2=0 (again, technically change 1/Ghat(z)=0 for these obs but still)
  #*
  #* then change z=Inf if y2< s and delta2=0 (again, technically change 1/Ghat(z)=0 for these obs but still)
  #* and that should do it! (I hope)

  ipcw_mat <- matrix(nrow=n,ncol=t_length,dimnames = list(paste0("i",1:n),paste0("t",t_cutoff)))
  for(t_ind in 1:t_length){
    #vector of min(ttilde,t)
    y_last <- pmin(y2,t_cutoff[t_ind])
    if(sum(y2 <= t_cutoff[t_ind] & delta2==1)>0){
      y_last[y2 <= t_cutoff[t_ind] & delta2==1] <- y_last[y2 <= t_cutoff[t_ind] & delta2==1] - 1e-8
    }
    y_last_cens <- rep(NA,n)
    y_last_cens[order(y_last)] <- summary(sfcens, times = y_last, extend=TRUE)$surv

    #now, to eliminate the the 'censored' observations
    if(sum(y_last <= t_cutoff[t_ind] & delta2==0) > 0){
      y_last_cens[y2 <= t_cutoff[t_ind] & delta2==0] <- Inf
    }
    #below is my former definition, but I will defer to graf and spitoni above
    # if(sum(y_last < t_cutoff[t_ind] & delta2==0) > 0){
    #   y_last_cens[y2 < t_cutoff[t_ind] & delta2==0] <- Inf
    # }
    ipcw_mat[,t_ind] <- 1/y_last_cens
  }

  ipcw_mat

}


#' Get inverse probability of censoring weights
#'
#' This function returns a vector of inverse probability of censoring weights from an unadjusted Cox model
#'   for censoring times.
#'
#' @inheritParams get_outcome_mat_logit
#' @inheritParams calc_risk
#'
#' @return a vector.
#' @export
get_ipcw_mat_logit <- function(y2,delta1,delta2,deltaD,t_cutoff){

  # browser()
  n <- length(y2)
  t_length <- length(t_cutoff)

  #updated event indicator for whether t2 is observed or not (aka, 1-censoring_indicator)
  delta2_fixed <- delta2 == 1 | (delta1==1 & deltaD==1)

  #this is Ghat, a non-parametric model of the 'survival distribution' of censoring var C.
  sfcens <- survival::survfit(survival::Surv(y2, delta2_fixed==0) ~ 1)

  #* want to get Ghat(z) where (following Graf 1999 'three categories')
  #* Category 1:
  #* z=y2- if y2<=s and delta2=1,
  #* Category 2:
  #* z=s if s<y2
  #* Category 3:
  #* z=Inf if y2<s and delta2=0, (aka, 1/Ghat(z)=0, I know they're subtly different)
  #* z=s if s=y2 and delta2=0, #this one situation is what I'm unsure about
  #* because graf and spitoni would say that if s=y2 and delta2=0, then z=Inf and result should be tossed.
  #* below I defer to them, but I'm just noting that to me there is logic in the other direction.
  #*
  #* so, we define
  #* z = min(y2,s)
  #* then change z=y2- if y2<=s and delta2=1
  #* then change z=Inf if y2<=s and delta2=0 (again, technically change 1/Ghat(z)=0 for these obs but still)
  #*
  #* then change z=Inf if y2< s and delta2=0 (again, technically change 1/Ghat(z)=0 for these obs but still)
  #* and that should do it! (I hope)

  ipcw_mat <- matrix(nrow=n,ncol=t_length,dimnames = list(paste0("i",1:n),paste0("t",t_cutoff)))
  for(t_ind in 1:t_length){
    #vector of min(ttilde,t)
    y_last <- pmin(y2,t_cutoff[t_ind])
    if(sum(y2 <= t_cutoff[t_ind] & delta2_fixed==1)>0){
      y_last[y2 <= t_cutoff[t_ind] & delta2_fixed==1] <- y_last[y2 <= t_cutoff[t_ind] & delta2_fixed==1] - 1e-8
    }
    y_last_cens <- rep(NA,n)
    y_last_cens[order(y_last)] <- summary(sfcens, times = y_last, extend=TRUE)$surv

    #now, to eliminate the the 'censored' observations
    if(sum(y_last <= t_cutoff[t_ind] & delta2_fixed==0) > 0){
      y_last_cens[y2 <= t_cutoff[t_ind] & delta2_fixed==0] <- Inf
    }
    #below is my former definition, but I will defer to graf and spitoni above
    # if(sum(y_last < t_cutoff[t_ind] & delta2_fixed==0) > 0){
    #   y_last_cens[y2 < t_cutoff[t_ind] & delta2_fixed==0] <- Inf
    # }
    ipcw_mat[,t_ind] <- 1/y_last_cens
  }

  ipcw_mat

}


#' Compute prediction performance score
#'
#' This function takes in all of the ingredients needed for prediction validation,
#'   and returns the corresponding scores.
#'
#' @param outcome_mat Output from get_outcome_mat function
#' @param pred_mat Output from calc_risks function
#' @param ipcw_mat Output from get_ipcw_mat function
#' @param score String indicating whether 'brier' score, or 'entropy' should be computed.
#'
#' @return a vector.
#' @export
compute_score <- function(outcome_mat, pred_mat, ipcw_mat, score="brier"){
  #this function is for brier and kl scores (aka cross entropy)
  #one weird thing is that in spitoni, the authors divide by n instead of by the sum of weights, is that right?
  # browser()
  if(length(dim(outcome_mat))==3){
    if(tolower(score) %in% c("brier")){
      out <- apply( t(apply((outcome_mat - pred_mat)^2,
                            MARGIN = c(1,3),
                            FUN = sum)) *
                      ipcw_mat, MARGIN = 2, FUN = mean)
    } else{
      out <- apply( t(apply(-outcome_mat*log(pred_mat),
                            MARGIN = c(1,3),
                            FUN = sum)) *
                      ipcw_mat, MARGIN = 2, FUN = mean)
    }
  } else{ #this must mean that there is only a single time cutoff, so the input mats are n by 4 and weights is just a vector
    if(tolower(score) %in% c("brier")){
      out <- colMeans( apply((outcome_mat - pred_mat)^2,
                            MARGIN = c(1),
                            FUN = sum) * ipcw_mat)
    } else{
      out <- colMeans( apply(-outcome_mat*log(pred_mat),
                         MARGIN = c(1),
                         FUN = sum) * ipcw_mat)
    }
  }
  return(out)

}



#' Compute IPCW-adjusted Hypervolume Under the Manifold
#'
#' This function computes the IPCW-adjusted Hypervolume under the manifold, as in Lee (unpublished).
#'
#' @inheritParams compute_score
#'
#' @return a scalar HUM
#' @export
compute_hum <- function(outcome_mat, pred_mat, ipcw_mat){
  #based on code from https://github.com/gaoming96/mcca/blob/master/R/hum.R
  #which I find a little hard to read because they put things in terms of the kronecker operation
  #but ultimately, things work out it seems
  #I added in the inverse probability weighting
  #maybe eventually I'll rewrite this so that it's neater

  if(length(dim(outcome_mat))==3){
    stop("for now, hum can only be computed at a single t_cutoff point")
  }

  pp <- pred_mat
  ipcw_vec <- as.vector(ipcw_mat)

  #first, let's just consider the case of t_cutoff being a single point, so everything here is matrices
  n <- nrow(outcome_mat)
  #n is the sample size
  a <- matrix(0,n,4)
  one1=a
  one1[,1]=1
  one2=a
  one2[,2]=1
  one3=a
  one3[,3]=1
  one4=a
  one4[,4]=1

  x1=which(outcome_mat[,1]==1)
  x2=which(outcome_mat[,2]==1)
  x3=which(outcome_mat[,3]==1)
  x4=which(outcome_mat[,4]==1)
  n1=length(x1)
  n2=length(x2)
  n3=length(x3)
  n4=length(x4)

  dd1=pp-one1
  dd2=pp-one2
  dd3=pp-one3
  dd4=pp-one4

  #here, the 'brier score' used to assess correct categorization is computed on the exp scale, and each term is square rooted
  jd1=exp(sqrt(dd1[,1]^2+dd1[,2]^2+dd1[,3]^2+dd1[,4]^2))
  jd2=exp(sqrt(dd2[,1]^2+dd2[,2]^2+dd2[,3]^2+dd2[,4]^2))
  jd3=exp(sqrt(dd3[,1]^2+dd3[,2]^2+dd3[,3]^2+dd3[,4]^2))
  jd4=exp(sqrt(dd4[,1]^2+dd4[,2]^2+dd4[,3]^2+dd4[,4]^2))

  # resulting matrix pattern is as follows (consider simple case of 2 individuals in each of 4 categories to establish pattern):
  #   1111  1211
  #   1112  1212
  #   1121  1221
  #   1122  1222
  #   2111  2211
  #   2112  2212
  #   2121  2221
  #   2122  2222
  mt1=kronecker(kronecker(jd1[x1]%*%t(jd2[x2]),jd3[x3]),jd4[x4]);
  mt7=kronecker(kronecker(jd1[x1]%*%t(jd2[x2]),jd4[x3]),jd3[x4]);
  mt2=kronecker(kronecker(jd1[x1]%*%t(jd3[x2]),jd2[x3]),jd4[x4]);
  mt8=kronecker(kronecker(jd1[x1]%*%t(jd3[x2]),jd4[x3]),jd2[x4]);
  mt9 =kronecker(kronecker(jd1[x1]%*%t(jd4[x2]),jd2[x3]),jd3[x4]);
  mt10=kronecker(kronecker(jd1[x1]%*%t(jd4[x2]),jd3[x3]),jd2[x4]);
  mt3 =kronecker(kronecker(jd2[x1]%*%t(jd1[x2]),jd3[x3]),jd4[x4]);
  mt11=kronecker(kronecker(jd2[x1]%*%t(jd1[x2]),jd4[x3]),jd3[x4]);
  mt4 =kronecker(kronecker(jd2[x1]%*%t(jd3[x2]),jd1[x3]),jd4[x4]);
  mt12=kronecker(kronecker(jd2[x1]%*%t(jd3[x2]),jd4[x3]),jd1[x4]);
  mt13=kronecker(kronecker(jd2[x1]%*%t(jd4[x2]),jd3[x3]),jd1[x4]);
  mt14=kronecker(kronecker(jd2[x1]%*%t(jd4[x2]),jd1[x3]),jd3[x4]);
  mt5 =kronecker(kronecker(jd3[x1]%*%t(jd2[x2]),jd1[x3]),jd4[x4]);
  mt15=kronecker(kronecker(jd3[x1]%*%t(jd2[x2]),jd4[x3]),jd1[x4]);
  mt6 =kronecker(kronecker(jd3[x1]%*%t(jd1[x2]),jd2[x3]),jd4[x4]);
  mt16=kronecker(kronecker(jd3[x1]%*%t(jd1[x2]),jd4[x3]),jd2[x4]);
  mt23=kronecker(kronecker(jd3[x1]%*%t(jd4[x2]),jd1[x3]),jd2[x4]);
  mt24=kronecker(kronecker(jd3[x1]%*%t(jd4[x2]),jd2[x3]),jd1[x4]);
  mt17=kronecker(kronecker(jd4[x1]%*%t(jd1[x2]),jd2[x3]),jd3[x4]);
  mt18=kronecker(kronecker(jd4[x1]%*%t(jd1[x2]),jd3[x3]),jd2[x4]);
  mt19=kronecker(kronecker(jd4[x1]%*%t(jd2[x2]),jd1[x3]),jd3[x4]);
  mt20=kronecker(kronecker(jd4[x1]%*%t(jd2[x2]),jd3[x3]),jd1[x4]);
  mt21=kronecker(kronecker(jd4[x1]%*%t(jd3[x2]),jd1[x3]),jd2[x4]);
  mt22=kronecker(kronecker(jd4[x1]%*%t(jd3[x2]),jd2[x3]),jd1[x4]);

  #now, to compute the product of the weights for each of the four individuals at each entry
  weight_mat=kronecker(kronecker(ipcw_vec[x1]%*%t(ipcw_vec[x2]),ipcw_vec[x3]),ipcw_vec[x4]);

  #finally, to compute numerator, multiply by 1 to turn it into numeric matrix
  # num <- as.numeric(1e-10 > abs(mt1-pmin(mt7,pmin(mt8,pmin(mt9,pmin(mt10,pmin(mt11,pmin(mt12,pmin(mt13,pmin(mt14,pmin(mt15,pmin(mt16,pmin(mt17,pmin(mt18,pmin(mt19,pmin(mt20,pmin(mt21,pmin(mt22,pmin(mt23,pmin(mt24,pmin(pmin(pmin(pmin(pmin(mt1, mt2), mt3), mt4), mt5), mt6)))))))))))))))))))))
  num <- 1*(mt1==pmin(mt7,pmin(mt8,pmin(mt9,pmin(mt10,pmin(mt11,pmin(mt12,pmin(mt13,pmin(mt14,pmin(mt15,pmin(mt16,pmin(mt17,pmin(mt18,pmin(mt19,pmin(mt20,pmin(mt21,pmin(mt22,pmin(mt23,pmin(mt24,pmin(pmin(pmin(pmin(pmin(mt1, mt2), mt3), mt4), mt5), mt6))))))))))))))))))))
  hum <- sum(num*weight_mat)/sum(weight_mat)

  # cr=sum(mt1==pmin(mt7,pmin(mt8,pmin(mt9,pmin(mt10,pmin(mt11,pmin(mt12,pmin(mt13,pmin(mt14,pmin(mt15,pmin(mt16,pmin(mt17,pmin(mt18,pmin(mt19,pmin(mt20,pmin(mt21,pmin(mt22,pmin(mt23,pmin(mt24,pmin(pmin(pmin(pmin(pmin(mt1, mt2), mt3), mt4), mt5), mt6))))))))))))))))))));
  #hypervolume under ROC manifold
  # hum=sum(num)/(n1*n2*n3*n4);

  return(hum)
}


#' Compute IPCW-adjusted Polytomous Discrimination Index
#'
#' This function computes the IPCW-adjusted polytomous discrimination index. Proposed in Van Caster and presented in Li
#'
#' @inheritParams compute_score
#'
#' @return a scalar HUM
#' @export
compute_pdi <- function(outcome_mat, pred_mat, ipcw_mat){
  #based on code from https://github.com/gaoming96/mcca/blob/master/R/pdi.R
  #which I find a little hard to read because they put things in terms of the kronecker operation
  #but ultimately, things work out it seems
  #I added in the inverse probability weighting
  #maybe eventually I'll rewrite this so that it's neater
  if(length(dim(outcome_mat))==3){
    stop("for now, pdi can only be computed at a single t_cutoff point")
  }
  #flag which subjects ended in which outcome category
  x1=which(outcome_mat[,1]==1)
  x2=which(outcome_mat[,2]==1)
  x3=which(outcome_mat[,3]==1)
  x4=which(outcome_mat[,4]==1)
  n1=length(x1)
  n2=length(x2)
  n3=length(x3)
  n4=length(x4)

  #vector of ipcw weights, and then individual vectors corresponding to each outcome category
  ipcw_vec <- as.vector(ipcw_mat)
  ipcw1 <- ipcw_vec[x1]
  ipcw2 <- ipcw_vec[x2]
  ipcw3 <- ipcw_vec[x3]
  ipcw4 <- ipcw_vec[x4]

  #individual specific vectors with 0 weighted elements removed (to speed up computation)
  ipcw1_nozero <- ipcw1[ipcw1 != 0]
  ipcw2_nozero <- ipcw2[ipcw2 != 0]
  ipcw3_nozero <- ipcw3[ipcw3 != 0]
  ipcw4_nozero <- ipcw4[ipcw4 != 0]

  #now, goal is to sum the product of every combination of subjects from the four categories
  pdi_denom <- sum(kronecker(kronecker(ipcw1_nozero %*% t(ipcw2_nozero),ipcw3_nozero),ipcw4_nozero))
  # ipcw234_mat <- as.matrix(expand.grid(ipcw2_nozero,ipcw3_nozero,ipcw4_nozero))
  # ipcw134_mat <- as.matrix(expand.grid(ipcw1_nozero,ipcw3_nozero,ipcw4_nozero))
  # ipcw124_mat <- as.matrix(expand.grid(ipcw1_nozero,ipcw2_nozero,ipcw4_nozero))
  # ipcw123_mat <- as.matrix(expand.grid(ipcw1_nozero,ipcw2_nozero,ipcw3_nozero))
  # sum(kronecker(X = ipcw1_nozero,Y = apply(X = ipcw234_mat,MARGIN = 1,FUN = prod)))

  #separate matrices for subjects in each of the four observed categories
  #use "drop = FALSE" to suppress conversion of single row to vector
  #https://stackoverflow.com/questions/12601692/r-and-matrix-with-1-row
  pv=pred_mat
  pv1=pv[x1,,drop = FALSE]
  pv2=pv[x2,,drop = FALSE]
  pv3=pv[x3,,drop = FALSE]
  pv4=pv[x4,,drop = FALSE]

  #approach here is to compute PDI_i as in Van Caster (2012) SIM paper
  # pdi1 <- pdi2 <- pdi3 <- pdi4 <- 0
  pdi4_num <- pdi3_num <- pdi2_num <- pdi1_num <- 0

  # there are a total of n_1 * n_2 * n_3 * n_4 combinations of members of the four groups
  # out of all of those comparisons, in how many is p1 for the person in category 1 the highest?
  # in order for the p1 from cat 1 to be the highest, it must be higher than each other one in turn,
  # so, one computational simplification is to multiply the number of people in cat 2 who are lower, cat 3, and cat 4
  # because multiplying corresponds with the number of total combinations, and thus, the total count.
  # for(i in 1:n1){
  #   pdi1=pdi1+sum(pv1[i,1]>pv2[,1])*sum(pv1[i,1]>pv3[,1])*sum(pv1[i,1]>pv4[,1])
  # }
  # for(i in 1:n2){
  #   pdi2=pdi2+sum(pv2[i,2]>pv1[,2])*sum(pv2[i,2]>pv3[,2])*sum(pv2[i,2]>pv4[,2])
  # }
  # for(i in 1:n3){
  #   pdi3=pdi3+sum(pv3[i,3]>pv1[,3])*sum(pv3[i,3]>pv2[,3])*sum(pv3[i,3]>pv4[,3])
  # }
  # for(i in 1:n4){
  #   pdi4=pdi4+sum(pv4[i,4]>pv1[,4])*sum(pv4[i,4]>pv2[,4])*sum(pv4[i,4]>pv3[,4])
  # }

  #for IPCW, we instead compute the weighted numerators separately
  for(i in 1:n1){
    #create four vectors, with the ipcw weights for each of
    #the relevant observations from each category that contribute to the PDI computation
    outvec1 <- ipcw1[i]
    if(outvec1 == 0){next}
    outvec2 <- ipcw2[pv1[i,1]>pv2[,1]]
    outvec2 <- outvec2[outvec2 != 0]
    outvec3 <- ipcw3[pv1[i,1]>pv3[,1]]
    outvec3 <- outvec3[outvec3 != 0]
    outvec4 <- ipcw4[pv1[i,1]>pv4[,1]]
    outvec4 <- outvec4[outvec4 != 0]

    #next, rather than simply multiplying the number of contributors from each category as before,
    #we need to actually multiply the weights for each combination of contributors, and sum those quantities up
    pdi1_num <- pdi1_num + outvec1*sum(apply(X = as.matrix(expand.grid(outvec2,outvec3,outvec4)),
                                     MARGIN = 1,
                                     FUN = prod))
  }
  for(i in 1:n2){
    outvec2 <- ipcw2[i]
    if(outvec2 == 0){next}
    outvec1 <- ipcw1[pv2[i,2]>pv1[,2]]
    outvec1 <- outvec1[outvec1 != 0]
    outvec3 <- ipcw3[pv2[i,2]>pv3[,2]]
    outvec3 <- outvec3[outvec3 != 0]
    outvec4 <- ipcw4[pv2[i,2]>pv4[,2]]
    outvec4 <- outvec4[outvec4 != 0]
    pdi2_num <- pdi2_num + outvec2*sum(apply(X = as.matrix(expand.grid(outvec1,outvec3,outvec4)),
                                             MARGIN = 1,
                                             FUN = prod))
  }
  for(i in 1:n3){
    outvec3 <- ipcw3[i]
    if(outvec3 == 0){next}
    outvec1 <- ipcw1[pv3[i,3]>pv1[,3]]
    outvec1 <- outvec1[outvec1 != 0]
    outvec2 <- ipcw2[pv3[i,3]>pv2[,3]]
    outvec2 <- outvec2[outvec2 != 0]
    outvec4 <- ipcw4[pv3[i,3]>pv4[,3]]
    outvec4 <- outvec4[outvec4 != 0]
    pdi3_num <- pdi3_num + outvec3*sum(apply(X = as.matrix(expand.grid(outvec1,outvec2,outvec4)),
                                             MARGIN = 1,
                                             FUN = prod))
  }
  for(i in 1:n4){
    outvec4 <- ipcw4[i]
    if(outvec4 == 0){next}
    outvec1 <- ipcw1[pv4[i,4]>pv1[,4]]
    outvec1 <- outvec1[outvec1 != 0]
    outvec2 <- ipcw2[pv4[i,4]>pv2[,4]]
    outvec2 <- outvec2[outvec2 != 0]
    outvec3 <- ipcw3[pv4[i,4]>pv3[,4]]
    outvec3 <- outvec3[outvec3 != 0]
    pdi4_num <- pdi4_num + outvec4*sum(apply(X = as.matrix(expand.grid(outvec1,outvec2,outvec3)),
                                             MARGIN = 1,
                                             FUN = prod))
  }

  # pdi<-(pdi1+pdi2+pdi3+pdi4)/(4*n1*n2*n3*n4)
  pdi <- (pdi1_num + pdi2_num + pdi3_num + pdi4_num) / (4*pdi_denom)

  return(pdi)
}


#' Compute IPCW-adjusted Correct Classification Probability
#'
#' This function computes the IPCW-adjusted Correct Classification Probability. Presented in Li (2019)
#'
#' @inheritParams compute_score
#'
#' @return a scalar CCP
#' @export
compute_ccp <- function(outcome_mat, pred_mat, ipcw_mat){

  if(length(dim(outcome_mat))==3){
    stop("for now, ccp can only be computed at a single t_cutoff point")
  }

  #compute the ccp directly as the average of correct classification overall
    #this implicitly weights the categories by their prevalence, which is how the CCP was
  #We could instead first compute ccp_m for each of the categories,
  #and then take a different weighted average of those, but we don't.

  #left of == is n-length vector of predicted prob for observed outcome
  #right of == is n-length max predicted prob for each obs
  #result is an n-length logical vector
  correct_ind <- 1*(t(pred_mat)[t(outcome_mat)==1] == apply(X=pred_mat,MARGIN = 1,FUN = max))
  ipcw_vec <- as.vector(ipcw_mat)
  stopifnot(length(correct_ind)==length(ipcw_vec))
  #proportion of correct classifications, weighted by ipcw!
  return(mean(ipcw_vec*correct_ind))

  #y is the quandr-nomial response, i.e., a single vector taking three distinct values, can be nominal or numerical
  #d is the continuous marker, turn out to be the probability matrix when option="prob"
  #flag which subjects ended in which outcome category
  # y <- apply(X = outcome_mat,MARGIN = 1,FUN = function(x){which(x==1)}) #generate vector of categories

  #flag for each subject which category has the maximum probability
  # pv <- max.col(pred_mat)
  # #figure out how to handle ties in the predicted probabilities....
  # l=max.col(pred_mat)
  # for(i in 1:nrow(pred_mat)){
  #   if (length(which(max(pred_mat[i,])==pred_mat[i,]))>1){
  #     cat("WARNING: there exists two same max probability in one sample.\n")
  #     break
  #   }
  # }
  # ccp=(sum(pv==1 & y==1)+ sum(pv==2 & y==2)+sum(pv==3 & y==3)+sum(pv==4 & y==4))/nn
  # return(ccp)

}

#' Compute IPCW-adjusted Net Reclassification Improvement
#'
#' This function computes the IPCW-adjusted Net Reclassification Probability. Presented as S in Li (2013) and Wang (2020)
#'
#' @inheritParams compute_score
#' @param pred_mat1 matrix of predicted probabilities with model 1 ("old" model)
#' @param pred_mat2 matrix of predicted probabilities with model 2 ("new" model)
#'
#' @return a scalar NRI
#' @export
compute_nri <- function(outcome_mat, pred_mat1, pred_mat2, ipcw_mat){

  #for a prevalence-weighted result, theres the standard difference of ccp (this is presented in Li 2019 and mcca package)
  # ccp_model1 <- compute_ccp(outcome_mat=outcome_mat, pred_mat=pred_mat1,ipcw_mat=ipcw_mat)
  # ccp_model2 <- compute_ccp(outcome_mat=outcome_mat, pred_mat=pred_mat2,ipcw_mat=ipcw_mat)
  # nri <- ccp_model2-ccp_model1

  #for an equal-weighted version, need to compute ccp separately across categories
  #left of == is n-length vector of predicted prob for observed outcome
  #right of == is n-length max predicted prob for each obs
  #result is an n-length logical vector
  correct_ind1 <- 1*(t(pred_mat1)[t(outcome_mat)==1] == apply(X=pred_mat1,MARGIN = 1,FUN = max))
  correct_ind2 <- 1*(t(pred_mat2)[t(outcome_mat)==1] == apply(X=pred_mat2,MARGIN = 1,FUN = max))
  ipcw_vec <- as.vector(ipcw_mat)
  stopifnot(length(correct_ind1)==length(ipcw_vec))
  y <- apply(X = outcome_mat,MARGIN = 1,FUN = function(x){which(x==1)}) #generate vector of categories
  ccp_vec_model1 <- tapply(X = correct_ind1*ipcw_vec,INDEX = y,FUN = mean)
  ccp_vec_model2 <- tapply(X = correct_ind2*ipcw_vec,INDEX = y,FUN = mean)
  nri <- mean(ccp_vec_model2-ccp_vec_model1)

  return(nri)

}


#' Compute AUC for Non-terminal and Terminal outcomes
#'
#' Function to compute univariate prediction statistics
#'
#'
#' @inheritParams compute_score
#' @inheritParams calc_risk
#' @param dat Data
#'
#' @return a vector with the non-terminal and terminal AUC.
#' @export
compute_auc <- function(dat,t_cutoff, pred_mat){
  #For now, this one is only implemented for when these are just matrices, aka for a single choice of t.
  #need pred_mat to be from same value as t_cutoff, however!


  # browser()
  nonterm_comp_risk_time <- ifelse(dat$y1 < dat$y2, dat$y1, dat$y2)
  comp_risk_event <- ifelse( (dat$y1 == dat$y2) & dat$delta2==1,2,ifelse(dat$y1 == dat$y2 & dat$delta2==0,0,1))
  outcomes <- cbind(dat[,c("y1","delta1","y2","delta2")],
                    nonterm_comp_risk_time=nonterm_comp_risk_time,
                    comp_risk_event=comp_risk_event)

  # #This former approach used dplyr but we don't need it!
  # #Now, I'm just testing another commit
  # outcomes <- dat %>% dplyr::select(y1,delta1,y2,delta2) %>%
  #   dplyr::mutate(nonterm_comp_risk_time = ifelse(y1 < y2, y1, y2),
  #          comp_risk_event = ifelse( (y1 == y2) & delta2==1,2,ifelse(y1 == y2 & delta2==0,0,1))
  #   )

  #treats terminal outcome as a competing risk
  ROC_nonterm <- timeROC::timeROC(T=outcomes$nonterm_comp_risk_time,
                           delta=outcomes$comp_risk_event,
                           marker=pred_mat[,"p_ntonly"] + pred_mat[,"p_both"],
                           cause=1,weighting="marginal",
                           times=t_cutoff,
                           iid=TRUE)

  ROC_term <- timeROC::timeROC(T=outcomes$y2,
                        delta=outcomes$delta2,
                        marker=pred_mat[,"p_tonly"] + pred_mat[,"p_both"],
                        cause=1,weighting="marginal",
                        times=t_cutoff,
                        iid=TRUE)

  return(c(AUC_nonterm = ROC_nonterm$AUC_1[2],
           AUC_term = ROC_term$AUC[2]))
}


#' Compute AUC for Non-terminal and Terminal outcomes
#'
#' Function to compute univariate prediction statistics
#'
#'
#' @inheritParams compute_score
#' @inheritParams calc_risk
#' @param dat Data
#'
#' @return a vector with the non-terminal and terminal AUC.
#' @export
compute_auc_logit <- function(dat,t_cutoff, pred_mat){
  #For now, this one is only implemented for when these are just matrices, aka for a single choice of t.
  #need pred_mat to be from same value as t_cutoff, however!

  # browser()

  #non-terminal event time treating terminal as competing risk
  nonterm_comp_risk_time <- pmin(dat$y1, dat$y2)
  comp_risk_event <- ifelse( dat$delta1==1, 1,
                             ifelse(dat$delta2==1, 2,
                             0))

  delta2_fixed <- as.numeric(dat$delta2 == 1 | (dat$delta1==1 & dat$deltaD==1))

  # #This former approach used dplyr but we don't need it!
  # #Now, I'm just testing another commit
  # outcomes <- dat %>% dplyr::select(y1,delta1,y2,delta2) %>%
  #   dplyr::mutate(nonterm_comp_risk_time = ifelse(y1 < y2, y1, y2),
  #          comp_risk_event = ifelse( (y1 == y2) & delta2==1,2,ifelse(y1 == y2 & delta2==0,0,1))
  #   )

  #treats terminal outcome as a competing risk
  ROC_nonterm <- timeROC::timeROC(T=nonterm_comp_risk_time,
                                  delta=comp_risk_event,
                                  marker=pred_mat[,"p_ntonly"] + pred_mat[,"p_both_noinst"]  + pred_mat[,"p_both_inst"],
                                  cause=1,weighting="marginal",
                                  times=t_cutoff,
                                  iid=TRUE)


  ROC_term <- timeROC::timeROC(T=dat$y2,
                               delta=delta2_fixed,
                               marker=pred_mat[,"p_tonly"] + pred_mat[,"p_both_noinst"]  + pred_mat[,"p_both_inst"],
                               cause=1,weighting="marginal",
                               times=t_cutoff,
                               iid=TRUE)

  return(c(AUC_nonterm = ROC_nonterm$AUC_1[2],
           AUC_term = ROC_term$AUC[2]))
}

