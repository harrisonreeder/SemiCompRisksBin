#' Bayesian Illness-Death Model with Weibull Baselines and Semi-Markov structure
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
Bayes_Uni <- function(Formula, data, na.action="na.fail", subset=NULL,
                      hyperParams, mcmcParams, n_chains, start_mat=NULL){
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
  time <- Formula::model.part(form2, data=data, lhs=1)
  y <- time[[1]]
  delta <- time[[2]]
  #Create covariate matrices for each of three transition hazards
  Xmat <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=1),
                                        data=data))
  p <- ncol(Xmat)
  n <- length(y)

  ####SET HYPERPARAMETERS####
  hyper_vec <- as.vector(c(hyperParams$WB$WB.ab, hyperParams$WB$WB.cd))
  tuning_vec <- as.vector(c(mhProp_alpha_var=mcmcParams$tuning$mhProp_alpha_var))

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
    start_vec <- c(kappa=0.1,alpha=1, if(p>0) stats::runif(p, -0.1, 0.1))
    start_mat <- matrix(data=start_vec,ncol=n_chains,nrow=length(start_vec))
    rownames(start_mat) <- c("kappa","alpha",if(p>0) paste0("beta_",1:p))
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
    mcmcRet <- WeibUnimcmc(y, delta, Xmat, hyper_vec, tuning_vec, start_mat[,i],
                           n_burnin, n_sample, thin)

    out_list[["samples"]][,i,1] <- mcmcRet[["samples"]][["kappa"]]
    out_list[["samples"]][,i,2] <- mcmcRet[["samples"]][["alpha"]]
    if(p>0) out_list[["samples"]][,i,3:(2+p)] <- mcmcRet[["samples"]][["beta"]]
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
