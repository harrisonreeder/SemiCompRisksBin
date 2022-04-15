#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include "PG_utilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List Logitmcmc_PG(const arma::vec& y, const arma::mat& Xmat,
                     const arma::vec& m0, const arma::mat& P0,
                     const arma::vec& start_vec,
                     int n_burnin, int n_sample, int thin){
  //timekeeping objects
  std::time_t newt;

  //set constants
  int p = Xmat.n_cols;
  int n = Xmat.n_rows;
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'
  int move; //index for which parameter to update
  int M, i; //counters for MCMC sampler and for inner loops
  int StoreInx; //index for where to store a sample, post-thinning

  //initialize starting values
  arma::vec beta = start_vec;
  arma::vec mean_const = Xmat.t() * (y-0.5) + P0 * m0;
  arma::vec omega = arma::ones(n);
  arma::vec eta = arma::zeros(n);
  arma::vec beta_mean = arma::zeros(p);
  arma::mat beta_var = arma::eye(p,p);

  //create storage objects
  arma::mat sample_beta = arma::mat(p,n_store,arma::fill::zeros); //store betas column by column for "speed" ?
  int accept_beta = 0;

  for(M = 0; M < n_iter; M++){

    eta = Xmat * beta;

    //draw auxiliary omega variables
    for(i = 0; i < n; i++){
      omega(i) = rpg1z_devroye(eta(i));
    }

    //compute mean and variance vectors for beta draw
    beta_var = inv_sympd(Xmat.t() * (Xmat.each_col() % omega) + P0);
    beta_mean = beta_var * mean_const;

    //draw beta (inplace)
    mvrnorm_inplace(beta,beta_mean,beta_var);
    accept_beta++; //this is dumb but it's a holdover from MH algorithm

    /* Storing posterior samples */
    if( ( (M+1) % thin ) == 0 && (M+1) > n_burnin){
      StoreInx = (M+1 - n_burnin)/thin;
      sample_beta.col(StoreInx - 1) = beta;
    }

    //check for user interruption, and report iterations
    if( ( (M+1) % 10000 ) == 0){
      newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      Rcpp::Rcout << "iteration: " << M+1 << ": " << ctime(&newt) << "\n";
      Rcpp::checkUserInterrupt(); //cheks if the user hit the "stop" icon to cancel running sampler.
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("samples") = Rcpp::List::create(
      Rcpp::Named("beta") = sample_beta.t()),
      Rcpp::Named("accept") = Rcpp::List::create(
        Rcpp::Named("beta") = accept_beta)
  );
}
