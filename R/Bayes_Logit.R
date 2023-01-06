#' Bayesian Logistic Regression using Polya-Gamma data augmentation
#'
#' @param y an n-length binary vector of outcomes
#' @param Xmat an n by p matrix of covariates
#' @param m0 prior mean vector
#' @param P0 prior variance matrix
#' @param numReps numeric total number of iterations
#' @param thin numeric number of iterations to take before each sample
#' @param burninPerc numeric proportion of numReps to toss as burn-in
#' @param n_chains numeric number of chains to initialize, if start_list is null
#' @param start_mat a matrix of start values, with each column corresponding to a chain. If null, n_chains sets number of chains
#'
#' @return a list with outputs
#' @export
Bayes_Logit <- function(y, Xmat,
                            m0=rep(0, ncol(Xmat)),
                            P0=matrix(0, nrow=ncol(Xmat), ncol=ncol(Xmat)),
                            numReps, thin, burninPerc,n_chains,start_mat) {
  # browser()

  #set index values
  n	<- length(y)
  p	<- ncol(Xmat) #for now, there must be at least one for the intercept

  #set mcmc parameters
  n_burnin <- numReps * burninPerc
  n_sample <- numReps - n_burnin
  n_store      <- numReps/thin * (1 - burninPerc)
  stopifnot(n_burnin %% 1 == 0)
  stopifnot(n_sample > 0)
  if(n_store %% 1 != 0){ stop("numReps * burninPerc  must be divisible by thin")}

  ####ASSIGN START VALUES####
  if(is.null(start_mat)){
    if(is.null(n_chains)){stop("either supply non-null start_list, or specify n_chains.")}
    #create an empty list with as many empty lists as there are chains
    start_mat <- matrix(stats::runif(n=n_chains*p,-0.1,0.1),nrow = p,ncol=n_chains)
  }

  out_list <- list(samples = array(dim = c(n_store, n_chains, p),
                                   dimnames = list(as.character(1:n_store),
                                                   paste0("chain:",1:n_chains),
                                                   rownames(start_mat))),
                   accept = vector(mode = "list",length=n_chains),
                   setup = list(m0=m0,P0=P0,
                                start_mat = start_mat,
                                mcmc_para = c(n_burnin=n_burnin,
                                              n_sample=n_sample,thin=thin)))
  names(out_list$accept) <- paste0("chain",1:n_chains)

  # #TODO: parallelize this loop (I know it can be done!)
  for(i in 1:n_chains){
    print(paste0("Chain: ", i))
    mcmcRet <- Logitmcmc_PG(y = y,
                         Xmat = Xmat,
                         m0 = m0, P0 = P0,
                         start_vec = start_mat[,i],
                         n_burnin = n_burnin,
                         n_sample	= n_sample,
                         thin = thin)

    out_list[["samples"]][,i,] <- mcmcRet[["samples"]][["beta"]]
    out_list[["accept"]][[paste0("chain",i)]] <- mcmcRet$accept
  }

  #for now, my plan is going to be to leverage the bayesplot package to
  #make visuals

  return(out_list)

}
