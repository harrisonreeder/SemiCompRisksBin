#include <RcppArmadillo.h>
#include <chrono>
#include <ctime>
#include "PG_utilities.h"
// [[Rcpp::depends(RcppArmadillo)]]


#define Pi 3.141592653589793238462643383280

//logsumexp trick

double logsumexpnew(double x,double y){
  double c = fmax(x,y);
  return c + log(exp(x-c) + exp(y-c));
}

void Logit_update_beta_frail(arma::vec& beta, arma::vec& eta, const arma::mat& Xmat,
                             const arma::vec& frail, const arma::vec& mean_const,
                             const arma::vec& P0diag, int& accept_beta){

  int n = Xmat.n_rows;
  int p = Xmat.n_cols;
  //initialize starting values
  arma::vec omega = arma::ones(n);
  arma::vec beta_mean = arma::zeros(p);
  arma::mat beta_var = arma::eye(p,p);
  arma::mat beta_prec = arma::eye(p,p);

  //draw auxiliary omega variables
  for(int i = 0; i < n; i++){
    omega(i) = rpg1z_devroye(eta(i) + log(frail(i))); //log-frailty included in ith linear predictor
  }

  //compute mean and variance vectors for beta draw
  beta_prec = Xmat.t() * (Xmat.each_col() % omega);
  beta_prec.diag() = beta_prec.diag() + P0diag;
  //beta_var = inv_sympd(Xmat.t() * (Xmat.each_col() % omega));
  //beta_var = inv_sympd(Xmat.t() * (Xmat.each_col() % omega) + P0); //if we had a prior P0
  beta_var = inv_sympd(beta_prec);
  beta_mean = beta_var * (mean_const - Xmat.t() * (omega % arma::log(frail)));

  //draw beta (inplace)
  mvrnorm_inplace(beta,beta_mean,beta_var);
  //update values
  eta = Xmat * beta;
  accept_beta++;
}








double logLikWB_uni(const arma::vec& y, const arma::uvec& delta,
                     double alpha, double kappa,
                     const arma::vec& eta, const arma::vec& frail){
  //delta * (log lambda) + log(survival)
  double obj_val = arma::accu( delta % (log(alpha) + log(kappa) + (alpha - 1) * arma::log(y)
                                        + eta + arma::log(frail))
                                - ( kappa * arma::pow(y,alpha) ) % arma::exp(eta) % frail );
  return(obj_val);
}

double logLikWB_uni_i(const double& yi, const int& deltai,
                    double alpha, double kappa,
                    const double& etai, const double& fraili){
  //delta * (log lambda) + log(survival)
  double obj_val = - ( kappa * pow(yi,alpha) ) * exp(etai) * fraili ;
  if(deltai==1){
    obj_val += log(alpha) + log(kappa) + (alpha - 1) * log(yi)
                + etai + log(fraili);
  }
  return(obj_val);
}



void Bweib_update_betaj(arma::vec& beta, arma::vec& eta,
                         double &alpha, double &kappa, const arma::vec& frail,
                         const arma::vec& y, const arma::uvec& delta,
                         const arma::mat& Xmat, arma::vec &accept_beta, double &curr_loglik){

  int p = Xmat.n_cols;
  int j = (int) R::runif(0, p);
  //for(j = 0; j < p; j++){
    //temp represents negative log lambda
    arma::vec temp = - kappa * frail % arma::pow(y,alpha) % arma::exp(eta) % Xmat.col(j);
    double D2 = arma::dot(temp, Xmat.col(j));
    double D1 = arma::accu( delta % Xmat.col(j) + temp );
    double beta_prop_me    = beta(j) - D1/D2;
    double beta_prop_var   = - 2.4*2.4/D2;
    double betaj_prop = R::rnorm(beta_prop_me, sqrt(beta_prop_var));

    arma::vec eta_prop = eta + Xmat.col(j) * (betaj_prop - beta(j));
    temp = temp % arma::exp(eta_prop - eta); //now temp becomes proposal negative log Lambda
    double D2_prop = arma::dot(temp, Xmat.col(j));
    double D1_prop = arma::accu( delta % Xmat.col(j) + temp );
    double loglik_prop = logLikWB_uni(y, delta, alpha, kappa, eta_prop, frail);

    double beta_prop_me_prop   = betaj_prop - D1_prop/D2_prop;
    double beta_prop_var_prop  = - 2.4*2.4/D2_prop;
    double logProp_IniToProp = R::dnorm(betaj_prop, beta_prop_me, sqrt(beta_prop_var), 1);
    double logProp_PropToIni = R::dnorm(beta(j), beta_prop_me_prop, sqrt(beta_prop_var_prop), 1);
    double logR = loglik_prop - curr_loglik + logProp_PropToIni - logProp_IniToProp;
    if( log(R::unif_rand()) < logR ){
      beta(j) = betaj_prop;
      eta = eta_prop;
      accept_beta(j) = accept_beta(j) + 1;
      curr_loglik = loglik_prop;
    }
  //}
  return;
}

void Bweib_update_alpha(const arma::vec& eta, double &alpha, double &kappa,
                         const arma::vec& frail, const arma::vec& y, const arma::uvec& delta,
                         const double &mhProp_alpha_var, const double &alpha_a, const double &alpha_b,
                         int &accept_alpha, double &curr_loglik){
  double alpha_prop = R::rgamma( alpha*alpha / mhProp_alpha_var,
                                 mhProp_alpha_var/alpha);
  double loglik_prop = logLikWB_uni(y, delta, alpha_prop, kappa, eta, frail);

  double logPrior = R::dgamma(alpha, alpha_a, 1/alpha_b, 1);
  double logPrior_prop = R::dgamma(alpha_prop, alpha_a, 1/alpha_b, 1);
  double logProp_PropToIni = R::dgamma(alpha, alpha_prop*alpha_prop/mhProp_alpha_var,
                                       mhProp_alpha_var/alpha_prop, 1);
  double logProp_IniToProp = R::dgamma(alpha_prop, alpha*alpha/mhProp_alpha_var,
                                       mhProp_alpha_var/alpha, 1);
  double logR = loglik_prop - curr_loglik + logPrior_prop - logPrior
                + logProp_PropToIni - logProp_IniToProp;
  if( log(R::unif_rand()) < logR ){
    alpha = alpha_prop;
    accept_alpha++;
    curr_loglik = loglik_prop;
  }
  return;
}

void Bweib_update_kappa(const arma::vec& eta, double &alpha, double &kappa,
                        const arma::vec& frail, const arma::vec& y,
                        const arma::uvec& delta, double &delta_sum,
                        const double &kappa_a, const double &kappa_b,
                        int &accept_kappa, double &curr_loglik){
  kappa = R::rgamma(delta_sum + kappa_a,
                     1/(arma::dot(frail, arma::pow(y,alpha) % arma::exp(eta)) + kappa_b) );
  curr_loglik = logLikWB_uni(y, delta, alpha, kappa, eta, frail);
  accept_kappa++;
}




// [[Rcpp::export]]
Rcpp::List WeibUnimcmc(const arma::vec &y, const arma::uvec &delta,
                       const arma::mat &Xmat,
                       const arma::vec &hyper_vec,
                       const arma::vec &tuning_vec,
                       const arma::vec &start_vec,
                       int n_burnin,
                       int n_sample,
                       int thin){

  //I'm trying something new, which is explicitly passing in data as three sets
  //of vectors, corresponding with the three transitions...

  //timekeeping objects
  std::time_t newt;

  //SET CONSTANTS
  int p = Xmat.n_cols;
  int n = y.n_rows;
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'

  //INITIALIZE HYPERPARAMETERS
  double alpha_a = hyper_vec[0];
  double alpha_b = hyper_vec[1];
  double kappa_a = hyper_vec[2];
  double kappa_b = hyper_vec[3];

  //INITIALIZE TUNING PARAMETERS
  double mhProp_alpha_var = tuning_vec[0];

  // Rcpp::Rcout << "finished setting MCMC tuning params" << "\n";


  //INITIALIZE STARTING VALUES
  double kappa=start_vec(0);
  double alpha=start_vec(1);
  arma::vec beta;
  arma::vec eta(n,arma::fill::zeros);
  if(p>0){
    beta = start_vec(arma::span(2,1+p));
    eta = Xmat * beta;
  }

  //CREATE STORAGE VECTORS/MATRICES FOR SAMPLING
  arma::vec sample_kappa = arma::vec(n_store,arma::fill::zeros);
  arma::vec sample_alpha = arma::vec(n_store,arma::fill::zeros);
  arma::mat sample_beta = arma::mat(p,n_store,arma::fill::zeros);

  int accept_kappa = 0;
  int accept_alpha = 0;
  arma::vec accept_beta = arma::vec(p,arma::fill::zeros);

  int StoreInx=0; //index for where to store a sample, post-thinning

  //INITIALIZE WORKING VALUES
  double delta_sum = arma::accu(delta);
  arma::vec frail = arma::ones(n);
  double curr_loglik = logLikWB_uni(y, delta, alpha, kappa, eta, frail);

  int numUpdate = 2;
  if(p > 0) numUpdate += 1;
  double prob_kappa = (double) 1/numUpdate;
  double prob_alpha = (double) 1/numUpdate;
  double prob_beta = 1 - prob_kappa - prob_alpha;
  int move; //index for which parameter to update

  //RUN MCMC
  // Rcpp::Rcout << "begin mcmc" << "\n";
  arma::uvec move_vec = arma::uvec(n_iter,arma::fill::zeros);
  for(int M = 0; M < n_iter; ++M){
    // Rcpp::Rcout << "iteration: " << M << "\n";

    move = (int) R::runif(0, numUpdate);
    move_vec(M) = move;

    //KAPPA updates
    if(move==0){
      kappa = R::rgamma(delta_sum + kappa_a,
                         1/(arma::dot(arma::pow(y,alpha), arma::exp(eta)) + kappa_b) );
      curr_loglik = logLikWB_uni(y, delta, alpha, kappa, eta, frail);
      accept_kappa++;
      // Rcpp::Rcout << "updated kappa: " << kappa << "\n";
    } else if(move==1){
      Bweib_update_alpha(eta, alpha, kappa, frail,    y,   delta,
                         mhProp_alpha_var, alpha_a, alpha_b,
                         accept_alpha, curr_loglik);
      // Rcpp::Rcout << "updated alpha: " << alpha << "\n";
    } else if(move==2){
      //BETA updates
      Bweib_update_betaj(beta, eta, alpha, kappa, frail, y, delta,
                         Xmat, accept_beta, curr_loglik);
      // Rcpp::Rcout << "updated beta: " << beta1.t() << "\n";
    }

    // Storing posterior samples
    if( ( (M+1) % thin ) == 0 && (M+1) > n_burnin)
    {
      StoreInx = (M+1 - n_burnin)/thin;

      sample_alpha(StoreInx - 1) = alpha;
      sample_kappa(StoreInx - 1) = kappa;
      sample_beta.col(StoreInx - 1) = beta;
    }

    if( ( (M+1) % 10000 ) == 0){
      newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      Rcpp::Rcout << "iteration: " << M+1 << ": " << ctime(&newt) << "\n";
      Rcpp::checkUserInterrupt(); //checks if the user hit the "stop" icon to cancel running sampler.
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("samples") = Rcpp::List::create(
      Rcpp::Named("alpha") = sample_alpha,
      Rcpp::Named("kappa") = sample_kappa,
      Rcpp::Named("beta") = sample_beta.t()),
    Rcpp::Named("accept") = Rcpp::List::create(
      Rcpp::Named("move") = move_vec,
      Rcpp::Named("alpha") = accept_alpha,
      Rcpp::Named("kappa") = accept_kappa,
      Rcpp::Named("beta") = accept_beta));
}




void BweibScrSM_update_theta_gamma(double &theta, const arma::vec &frail,
                                   const double &mhProp_theta_var,
                                   const double &theta_a, const double &theta_b,
                                   int &accept_theta){
  int n = frail.n_rows;
  double xi = 1/theta;

  // //alternative approach uses a newton step to center proposal distro
  // double logPost = arma::accu( arma::log(frail) - frail ) * xi
  //   + (n * xi + theta_a - 1) * log(xi)
  //   - xi * theta_b - n * std::lgamma(xi);
  // double D1 = arma::accu(arma::log(frail) - frail) + n * log(xi)
  //   + (n*xi + theta_a - 1)/xi - theta_b - n * R::digamma(xi);
  // double D2 = n/xi - (theta_a - 1)/(xi*xi) - n * R::trigamma(xi);
  // double xi_prop_me    = xi - D1/D2;
  // double xi_prop_var   = - 2.4*2.4/D2;
  // double xi_prop = R::rgamma( xi_prop_me*xi_prop_me/xi_prop_var,
  //                             xi_prop_var/xi_prop_me);
  // double logPost_prop = arma::accu( arma::log(frail) - frail ) * xi_prop
  //   + (n * xi_prop + theta_a - 1) * log(xi_prop)
  //   - xi_prop * theta_b - n * std::lgamma(xi_prop);
  // double D1_prop = arma::accu(arma::log(frail) - frail) + n * log(xi_prop)
  //   + (n*xi_prop + theta_a - 1)/xi_prop - theta_b - n * R::digamma(xi_prop);
  // double D2_prop = n/xi_prop - (theta_a - 1)/(xi_prop*xi_prop) - n * R::trigamma(xi_prop);
  // double xi_prop_me_prop    = xi_prop - D1_prop/D2_prop;
  // double xi_prop_var_prop   = - 2.4*2.4/D2_prop;
  // double logProp_PropToIni = R::dgamma(xi, xi_prop_me_prop*xi_prop_me_prop/xi_prop_var_prop,
  //                                      xi_prop_var_prop/xi_prop_me_prop, 1);
  // double logProp_IniToProp = R::dgamma(xi_prop, xi_prop_me*xi_prop_me/xi_prop_var,
  //                                      xi_prop_var/xi_prop_me, 1);

  double xi_prop = R::rgamma(xi*xi/mhProp_theta_var, mhProp_theta_var/xi);
  double logPost = arma::accu( arma::log(frail) - frail ) * xi
    + (n * xi + theta_a - 1) * log(xi)
    - xi * theta_b - n * std::lgamma(xi);
    double logPost_prop = arma::accu( arma::log(frail) - frail ) * xi_prop
      + (n * xi_prop + theta_a - 1) * log(xi_prop)
      - xi_prop * theta_b - n * std::lgamma(xi_prop);
      double logProp_PropToIni = R::dgamma(xi, xi_prop*xi_prop/mhProp_theta_var,
                                           mhProp_theta_var/xi_prop, 1);
      double logProp_IniToProp = R::dgamma(xi_prop, xi*xi/mhProp_theta_var,
                                           mhProp_theta_var/xi, 1);
      double logR = logPost_prop - logPost + logProp_PropToIni - logProp_IniToProp;
      if( log(R::unif_rand()) < logR ){
        theta = 1/xi_prop;
        accept_theta++;
      }
      return;
}



void BweibScrSM_update_frail_gamma(arma::vec &frail, arma::vec &frail_sm, const double &theta,
                             const arma::vec &eta1, const arma::vec &eta2, const arma::vec &eta3,
                             const double &kappa1, const double &alpha1,
                             const double &kappa2, const double &alpha2,
                             const double &kappa3, const double &alpha3,
                             const arma::vec &y1, const arma::vec &y_sm,
                             const arma::uvec &sm_ind, const arma::uvec &sm_ind_long,
                             const arma::uvec &delta1, const arma::uvec &delta1noD,
                             const arma::uvec &delta_cr,
                             const arma::uvec &delta_sm, int &accept_frail,
                             double &curr_loglik1, double &curr_loglik_cr, double &curr_loglik_sm){
  int n = y1.n_rows;
  double temp, temp2;
  for(int i = 0; i < n; i++){
    // Rcpp::Rcout << "gamma loop i: " << i << "\n";
    //temp is sum of cumulative cause-specific hazards
    temp = kappa1 * pow(y1(i), alpha1) * exp(eta1(i))
    + kappa2 * pow(y1(i), alpha2) * exp(eta2(i)) ;
    // Rcpp::Rcout << "temp done " << "\n";
    //temp2 is sum of how many events person experienced
    temp2 = delta1(i) + delta_cr(i);
    // Rcpp::Rcout << "temp2 done " << "\n";
    if(delta1noD(i)>0){ //if the non-terminal event has occurred and terminal event not immediate
      // Rcpp::Rcout << "temp update happening " << "\n";
      temp += kappa3 * pow(y_sm(sm_ind_long(i)), alpha3) * exp(eta3(sm_ind_long(i)));
      // Rcpp::Rcout << "temp2 update happening " << "\n";
      temp2 += delta_sm(sm_ind_long(i));
      // Rcpp::Rcout << "temp updates done " << "\n";
    }
    // Rcpp::Rcout << "temp: " << temp << "\n";
    // Rcpp::Rcout << "temp2: " << temp2 << "\n";
    frail(i) = R::rgamma( temp2 + 1/theta, 1/(temp + 1/theta));
  }
  frail_sm = frail(sm_ind);
  curr_loglik1 = logLikWB_uni(y1, delta1, alpha1, kappa1, eta1, frail);
  curr_loglik_cr = logLikWB_uni(y1, delta_cr, alpha2, kappa2, eta2, frail);
  curr_loglik_sm = logLikWB_uni(y_sm, delta_sm, alpha3, kappa3, eta3, frail_sm);
  accept_frail++;
}


void BweibScrSM_update_theta_ln(double &theta, const arma::vec &frail,
                                const double &theta_a, const double &theta_b,
                                int &accept_theta){
  int n = frail.n_rows;
  double temp = arma::norm(arma::log(frail));
  theta = 1 / R::rgamma(theta_a + n/2, 1/(theta_b + temp*temp/2));
  accept_theta++;
}

void BweibScrSM_update_frail_ln(arma::vec &frail, arma::vec &frail_sm, const double &theta,
                                   const arma::vec &eta1, const arma::vec &eta2, const arma::vec &eta3,
                                   const double &kappa1, const double &alpha1,
                                   const double &kappa2, const double &alpha2,
                                   const double &kappa3, const double &alpha3,
                                   const arma::vec &y1, const arma::vec &y_sm,
                                   const arma::uvec &sm_ind, const arma::uvec &sm_ind_long,
                                   const arma::uvec &delta1, const arma::uvec &delta1noD,
                                   const arma::uvec &delta_cr,
                                   const arma::uvec &delta_sm, const double& frail_prop_var, arma::vec &accept_frail,
                                   double &curr_loglik1, double &curr_loglik_cr, double &curr_loglik_sm){
  int n = y1.n_rows;
  double logliki,logliki_prop,fraili_prop,logprior,logprior_prop, logR;
  for(int i = 0; i < n; i++){
    //frail is on positives, so normal is drawing the log-frailty
    fraili_prop = exp(R::rnorm(log(frail(i)), sqrt(frail_prop_var)));

    logliki = logLikWB_uni_i(y1(i), delta1(i), alpha1, kappa1, eta1(i), frail(i))
    + logLikWB_uni_i(y1(i), delta_cr(i), alpha2, kappa2, eta2(i), frail(i));
    logliki_prop = logLikWB_uni_i(y1(i), delta1(i), alpha1, kappa1, eta1(i), fraili_prop)
      + logLikWB_uni_i(y1(i), delta_cr(i), alpha2, kappa2, eta2(i), fraili_prop);
    if(delta1noD(i)>0){ //if the non-terminal event has occurred
      logliki += logLikWB_uni_i(y_sm(sm_ind_long(i)), delta_sm(sm_ind_long(i)),
                                alpha3, kappa3, eta3(sm_ind_long(i)), frail(i));
      logliki_prop += logLikWB_uni_i(y_sm(sm_ind_long(i)), delta_sm(sm_ind_long(i)),
                                alpha3, kappa3, eta3(sm_ind_long(i)), fraili_prop);
    }

    logprior = R::dnorm(log(frail(i)), 0, sqrt(theta), 1);
    logprior_prop = R::dnorm(log(fraili_prop), 0, sqrt(theta), 1);

    logR = logliki_prop - logliki + logprior_prop - logprior;
    if( log(R::unif_rand()) < logR ){
      frail(i) = fraili_prop;
      accept_frail(i) = accept_frail(i) + 1;
    }
  }

  frail_sm = frail(sm_ind);
  curr_loglik1 = logLikWB_uni(y1, delta1, alpha1, kappa1, eta1, frail);
  curr_loglik_cr = logLikWB_uni(y1, delta_cr, alpha2, kappa2, eta2, frail);
  curr_loglik_sm = logLikWB_uni(y_sm, delta_sm, alpha3, kappa3, eta3, frail_sm);
}

void BweibScrSMlogit_update_frail_ln(arma::vec &frail, arma::vec &frail_sm, arma::vec &frailD, const double &theta,
                                const arma::vec &eta1, const arma::vec &eta2, const arma::vec &eta3, const arma::vec &etaD,
                                const double &kappa1, const double &alpha1,
                                const double &kappa2, const double &alpha2,
                                const double &kappa3, const double &alpha3,
                                const arma::vec &y1, const arma::vec &y_sm,
                                const arma::uvec &delta1, const arma::uvec &delta1noD,
                                const arma::uvec &delta_cr, const arma::uvec &delta_sm,
                                const arma::uvec &sm_ind, const arma::uvec &sm_ind_long, const arma::uvec &delta1_ind_long,
                                const double& frail_prop_var, arma::vec &accept_frail,
                                double &curr_loglik1, double &curr_loglik_cr, double &curr_loglik_sm){
  int n = y1.n_rows;
  double logliki,logliki_prop,fraili_prop,logprior,logprior_prop, logR;
  for(int i = 0; i < n; i++){
    //frail is on positives, so normal is drawing the log-frailty, then we exponentiate
    fraili_prop = exp(R::rnorm(log(frail(i)), sqrt(frail_prop_var)));

    logliki = logLikWB_uni_i(y1(i), delta1(i), alpha1, kappa1, eta1(i), frail(i))
      + logLikWB_uni_i(y1(i), delta_cr(i), alpha2, kappa2, eta2(i), frail(i));
    logliki_prop = logLikWB_uni_i(y1(i), delta1(i), alpha1, kappa1, eta1(i), fraili_prop)
      + logLikWB_uni_i(y1(i), delta_cr(i), alpha2, kappa2, eta2(i), fraili_prop);
    if(delta1(i)>0){ //nonterminal event has occurred
      //I'm writing binary log-likelihood contribution as in
      //slide 22 of https://www.stat.rutgers.edu/home/pingli/papers/Logit.pdf
      logliki += -log1p( exp(etaD(delta1_ind_long(i)) + log(frail(i))) );
      logliki_prop += -log1p( exp(etaD(delta1_ind_long(i)) + log(fraili_prop)) );
      if(delta1noD(i)>0){ //if the non-terminal event has occurred without immediate death
        logliki +=      logLikWB_uni_i(y_sm(sm_ind_long(i)), delta_sm(sm_ind_long(i)),
                                       alpha3, kappa3, eta3(sm_ind_long(i)), frail(i));
        logliki_prop += logLikWB_uni_i(y_sm(sm_ind_long(i)), delta_sm(sm_ind_long(i)),
                                       alpha3, kappa3, eta3(sm_ind_long(i)), fraili_prop);
      } else{ //then non-terminal event has occurred followed by immediate death
        //etaD cancels out from these, so just add the log-frailties
        logliki += log(frail(i));
        logliki_prop += log(fraili_prop);
      }
    }

    logprior = R::dnorm(log(frail(i)), 0, sqrt(theta), 1);
    logprior_prop = R::dnorm(log(fraili_prop), 0, sqrt(theta), 1);
    logR = logliki_prop - logliki + logprior_prop - logprior;
    if( log(R::unif_rand()) < logR ){
      frail(i) = fraili_prop;
      accept_frail(i) = accept_frail(i) + 1;
    }
  }

  //Rcpp::Rcout << "actually did all the frailty updates!" << "\n";

  frailD = frail( arma::find(delta1) );
  frail_sm = frail( sm_ind );

  //Rcpp::Rcout << "subsetted the frailties!" << "\n";

  curr_loglik1 = logLikWB_uni(y1, delta1, alpha1, kappa1, eta1, frail);
  curr_loglik_cr = logLikWB_uni(y1, delta_cr, alpha2, kappa2, eta2, frail);
  curr_loglik_sm = logLikWB_uni(y_sm, delta_sm, alpha3, kappa3, eta3, frail_sm);
}

//function to just compute the likelihood contributions of every subject in place

void BweibScrSMlogit_logLH_vec(arma::vec &logLH_vec, const arma::vec &frail,
                               const arma::vec &eta1, const arma::vec &eta2,
                               const arma::vec &eta3, const arma::vec &etaD,
                               const double &kappa1, const double &alpha1,
                               const double &kappa2, const double &alpha2,
                               const double &kappa3, const double &alpha3,
                               const arma::vec &y1, const arma::vec &y_sm,
                               const arma::uvec &delta1, const arma::uvec &delta1noD,
                               const arma::uvec &delta_cr, const arma::uvec &delta_sm,
                               const arma::uvec &sm_ind_long, const arma::uvec &delta1_ind_long){
  int n = y1.n_rows;
  double logliki;
  for(int i = 0; i < n; i++){
    logliki = logLikWB_uni_i(y1(i), delta1(i), alpha1, kappa1, eta1(i), frail(i))
      + logLikWB_uni_i(y1(i), delta_cr(i), alpha2, kappa2, eta2(i), frail(i));
    if(delta1(i)>0){ //nonterminal event has occurred
      //I'm writing binary log-likelihood contribution as in
      //slide 22 of https://www.stat.rutgers.edu/home/pingli/papers/Logit.pdf
      logliki += -log1p( exp(etaD(delta1_ind_long(i)) + log(frail(i))) );
      if(delta1noD(i)>0){ //if the non-terminal event has occurred without immediate death
        logliki +=      logLikWB_uni_i(y_sm(sm_ind_long(i)), delta_sm(sm_ind_long(i)),
                                       alpha3, kappa3, eta3(sm_ind_long(i)), frail(i));
      } else{ //then non-terminal event has occurred followed by immediate death
        //etaD cancels out from these, so just add the log-frailties
        logliki += log(frail(i));
      }
    }
    //update ith contribution
    logLH_vec(i) = logliki;
  }
}

//version that takes a single logfrailty value, that's a helper for below

void BweibScrSMlogit_logLH_marg_vec(arma::vec &logLH_marg_vec,
                               const arma::vec &eta1, const arma::vec &eta2,
                               const arma::vec &eta3, const arma::vec &etaD,
                               const double &kappa1, const double &alpha1,
                               const double &kappa2, const double &alpha2,
                               const double &kappa3, const double &alpha3, const double &theta,
                               const arma::vec &y1, const arma::vec &y_sm,
                               const arma::uvec &delta1, const arma::uvec &delta1noD,
                               const arma::uvec &delta_cr, const arma::uvec &delta_sm,
                               const arma::uvec &sm_ind_long, const arma::uvec &delta1_ind_long,
                               const arma::vec &gh_nodes, const arma::vec &gh_weights){
  int n = y1.n_rows;
  double logliki, logliki_marg, logfrail, pisqrt;

  //nodes are scaled by sqrt(2) * sigma * x + mu, but here sigma = sqrt(theta) and mu=0
  arma::vec gh_nodes_adj = sqrt(theta * 2) * gh_nodes;
  //note in gauss-hermite derivation there's an extra pi^{-.5} to fold into weights
  arma::vec log_gh_weights_adj = arma::log(gh_weights) - 0.5 * log(Pi);

  //loop through subjects
  for(int i = 0; i < n; i++){
    logliki_marg = 0;

    //loop through nodes
    for(int j = 0; j < gh_nodes.n_rows; j++){
      logfrail = gh_nodes_adj(j);

      logliki = logLikWB_uni_i(y1(i), delta1(i), alpha1, kappa1, eta1(i), exp(logfrail))
      + logLikWB_uni_i(y1(i), delta_cr(i), alpha2, kappa2, eta2(i), exp(logfrail));
      if(delta1(i)>0){ //nonterminal event has occurred
        //I'm writing binary log-likelihood contribution as in
        //slide 22 of https://www.stat.rutgers.edu/home/pingli/papers/Logit.pdf
        logliki += -log1p( exp(etaD(delta1_ind_long(i)) + logfrail) );
        if(delta1noD(i)>0){ //if the non-terminal event has occurred without immediate death
          logliki +=      logLikWB_uni_i(y_sm(sm_ind_long(i)), delta_sm(sm_ind_long(i)),
                                         alpha3, kappa3, eta3(sm_ind_long(i)), exp(logfrail));
        } else{ //then non-terminal event has occurred followed by immediate death
          //etaD cancels out from these, so just add the log-frailties
          logliki += logfrail;
        }
      }

      //now to accumulate the result using logsumexp trick
      //log Lmarg = log( sum_j( exp( log w_j + log Lcond_j ) ) )
      if(j == 0){
        logliki_marg = logliki + log_gh_weights_adj(j);
      } else{
        logliki_marg = logsumexpnew(logliki_marg, logliki + log_gh_weights_adj(j));
      }
    }
    //update ith contribution
    logLH_marg_vec(i) = logliki_marg;
  }

}






// [[Rcpp::export]]
Rcpp::List WeibSCRmcmc(const arma::vec &y1, const arma::vec &y_sm,
                       const arma::uvec &delta1,
                       const arma::uvec &delta1noD,
                       const arma::uvec &delta_cr,
                       const arma::uvec &delta_sm,
                       const arma::mat &Xmat1,
                       const arma::mat &Xmat2,
                       const arma::mat &Xmat3,
                       const arma::vec &hyper_vec,
                       const arma::vec &tuning_vec,
                       const arma::vec &start_vec,
                       int n_burnin,
                       int n_sample,
                       int thin,
                       const std::string frail_path = ""){

  //I'm trying something new, which is explicitly passing in data as three sets
  //of vectors, corresponding with the three transitions...

  //timekeeping objects
  std::time_t newt;

  //SET CONSTANTS
  int p1 = Xmat1.n_cols;
  int p2 = Xmat2.n_cols;
  int p3 = Xmat3.n_cols;
  int n = y1.n_rows;
  int n_sm = y_sm.n_rows;
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'

  //INITIALIZE HYPERPARAMETERS
  double alpha1_a = hyper_vec[0];
  double alpha1_b = hyper_vec[1];
  double alpha2_a = hyper_vec[2];
  double alpha2_b = hyper_vec[3];
  double alpha3_a = hyper_vec[4];
  double alpha3_b = hyper_vec[5];
  double kappa1_a = hyper_vec[6];
  double kappa1_b = hyper_vec[7];
  double kappa2_a = hyper_vec[8];
  double kappa2_b = hyper_vec[9];
  double kappa3_a = hyper_vec[10];
  double kappa3_b = hyper_vec[11];
  double theta_a  = hyper_vec[12];
  double theta_b  = hyper_vec[13];

  //INITIALIZE TUNING PARAMETERS
  double mhProp_alpha1_var = tuning_vec[0];
  double mhProp_alpha2_var = tuning_vec[1];
  double mhProp_alpha3_var = tuning_vec[2];
  double mhProp_theta_var  = tuning_vec[3];

  double frail_prop_var = 0.3;

  // Rcpp::Rcout << "finished setting MCMC tuning params" << "\n";

  //INITIALIZE STARTING VALUES
  double kappa1=start_vec(0);
  double alpha1=start_vec(1);
  double kappa2=start_vec(2);
  double alpha2=start_vec(3);
  double kappa3=start_vec(4);
  double alpha3=start_vec(5);
  double theta=start_vec(6);
  arma::vec beta1, beta2, beta3;
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n_sm,arma::fill::zeros);
  if(p1>0){
    beta1 = start_vec(arma::span(7,6+p1));
    eta1 = Xmat1 * beta1;
  }
  if(p2>0){
    beta2 = start_vec(arma::span(7+p1,6+p1+p2));
    eta2 = Xmat2 * beta2;
  }
  if(p3>0){
    beta3 = start_vec(arma::span(7+p1+p2,6+p1+p2+p3));
    eta3 = Xmat3 * beta3;
  }
  arma::vec frail = start_vec(arma::span(7+p1+p2+p3,6+p1+p2+p3+n));

  //grab subset of frails corresponding with h3
  arma::uvec sm_ind = arma::find( delta1noD );
  arma::vec frail_sm = frail(sm_ind);
  //now make a lookup vector connecting the index out of n with the index out of n_sm
  arma::uvec sm_ind_long = arma::uvec(n,arma::fill::zeros);
  // Rcpp::Rcout << "sm_ind_long(sm_ind): " << sm_ind_long(sm_ind) << "\n";
  // Rcpp::Rcout << "arma::linspace<arma::uvec>(0, n_sm-1,n_sm): " << arma::linspace<arma::uvec>(0, n_sm-1,n_sm) << "\n";
  sm_ind_long(sm_ind) = arma::linspace<arma::uvec>(0, n_sm-1,n_sm);

  // for(int j=0; j<n;j++){
  //   Rcpp::Rcout << "iteration j: " << j << "\n";
  //   Rcpp::Rcout << "delta1 status: " << delta1(j) << "\n";
  //   if(delta1(j)==1){
  //     Rcpp::Rcout << "corresponding y_sm value: " << y_sm(sm_ind_long(j)) << "\n";
  //   }
  // }

  // Rcpp::Rcout << "cbind(delta1, sm_ind): " << arma::join_horiz(delta1,sm_ind_long) << "\n";
  // Rcpp::Rcout << "cbind(delta1, frail): " << arma::join_horiz(arma::conv_to<arma::vec>::from(delta1(arma::span(0,20))),frail(arma::span(0,20))) << "\n";
  // Rcpp::Rcout << "frail_sm: " << frail_sm(arma::span(0,20)) << "\n";



  // Rcpp::Rcout << "cbind(delta1, sm_ind): " << arma::join_horiz(delta1,sm_ind_long) << "\n";
  //Rcpp::Rcout << "sm_ind: " << sm_ind(arma::span(0,10)) << "\n";
  //Rcpp::Rcout << "frail_sm: " << frail_sm(sm_ind(arma::span(0,10))) << "\n";
  // Rcpp::Rcout << "finished setting starting values" << "\n";


  //CREATE STORAGE VECTORS/MATRICES FOR SAMPLING
  arma::vec sample_kappa1 = arma::vec(n_store,arma::fill::zeros);
  arma::vec sample_alpha1 = arma::vec(n_store,arma::fill::zeros);
  arma::mat sample_beta1 = arma::mat(p1,n_store,arma::fill::zeros);
  arma::vec sample_kappa2 = arma::vec(n_store,arma::fill::zeros);
  arma::vec sample_alpha2 = arma::vec(n_store,arma::fill::zeros);
  arma::mat sample_beta2 = arma::mat(p2,n_store,arma::fill::zeros);
  arma::vec sample_kappa3 = arma::vec(n_store,arma::fill::zeros);
  arma::vec sample_alpha3 = arma::vec(n_store,arma::fill::zeros);
  arma::mat sample_beta3 = arma::mat(p3,n_store,arma::fill::zeros);
  arma::vec sample_theta = arma::vec(n_store,arma::fill::zeros);
  arma::mat sample_frail;
  if(frail_path.size() > 0){
    sample_frail = arma::mat(n,n_store,arma::fill::zeros);
  }

  int accept_kappa1 = 0;
  int accept_alpha1 = 0;
  int accept_kappa2 = 0;
  int accept_alpha2 = 0;
  int accept_kappa3 = 0;
  int accept_alpha3 = 0;
  int accept_theta = 0;
  arma::vec accept_beta1 = arma::vec(p1,arma::fill::zeros);
  arma::vec accept_beta2 = arma::vec(p2,arma::fill::zeros);
  arma::vec accept_beta3 = arma::vec(p3,arma::fill::zeros);

  arma::vec accept_frail = arma::vec(n,arma::fill::zeros);
  //int accept_frail = 0;

  int StoreInx=0; //index for where to store a sample, post-thinning

  //INITIALIZE WORKING VALUES
  double delta1_sum = arma::accu(delta1);
  double delta_cr_sum = arma::accu(delta_cr);
  double delta_sm_sum = arma::accu( delta_sm );
  //keep running log-likelihoods for each submodel
  double curr_loglik1 = logLikWB_uni(y1, delta1, alpha1, kappa1, eta1, frail);
  double curr_loglik_cr = logLikWB_uni(y1, delta_cr, alpha2, kappa2, eta2, frail);
  double curr_loglik_sm = logLikWB_uni(y_sm, delta_sm, alpha3, kappa3, eta3, frail_sm);
  // double curr_loglik = curr_loglik1 + curr_loglik_cr + curr_loglik_sm;
  // double temp = 0; //used for updating frails
  // double temp2 = 0; //used for updating frails

  int numUpdate = 8;
  if(p1 > 0) numUpdate += 1;
  if(p2 > 0) numUpdate += 1;
  if(p3 > 0) numUpdate += 1;
  double prob_beta1 = (p1 > 0) ? (double) 1/numUpdate : 0;
  double prob_beta2 = (p2 > 0) ? (double) 1/numUpdate : 0;
  double prob_beta3 = (p3 > 0) ? (double) 1/numUpdate : 0;
  double prob_kappa1 = (double) 1/numUpdate;
  double prob_kappa2 = (double) 1/numUpdate;
  double prob_kappa3 = (double) 1/numUpdate;
  double prob_alpha1 = (double) 1/numUpdate;
  double prob_alpha2 = (double) 1/numUpdate;
  double prob_alpha3 = (double) 1/numUpdate;
  double prob_frail  = (double) 1/numUpdate;
  double prob_theta = 1 - prob_beta1 - prob_beta2 - prob_beta3
                        - prob_kappa1 - prob_kappa2 - prob_kappa3
                        - prob_alpha1 - prob_alpha2 - prob_alpha3
                        - prob_frail;

  int move; //index for which parameter to update

  // Rcpp::Rcout << "Path: " << frail_path << "\n" ;
  // Rcpp::Rcout << "Length of path: " << frail_path.size() << "\n" ;


  //RUN MCMC
  newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  Rcpp::Rcout << "Begin MCMC: " << ctime(&newt) << "\n";

  arma::uvec move_vec = arma::uvec(n_iter,arma::fill::zeros);
  for(int M = 0; M < n_iter; ++M){
    // Rcpp::Rcout << "iteration: " << M << "\n";

    move = (int) R::runif(0, numUpdate);
    move_vec(M) = move;

    //KAPPA updates
    if(move==0){ //kappa1
      Bweib_update_kappa(eta1, alpha1, kappa1, frail, y1,
                         delta1, delta1_sum, kappa1_a, kappa1_b,
                         accept_kappa1, curr_loglik1);

    // Rcpp::Rcout << "updated kappa1: " << kappa1 << "\n";
    } else if(move==1){ //kappa2
      Bweib_update_kappa(eta2, alpha2, kappa2, frail, y1,
                         delta_cr, delta_cr_sum, kappa2_a, kappa2_b,
                         accept_kappa2, curr_loglik_cr);

    // Rcpp::Rcout << "updated kappa2: " << kappa2 << "\n";
    } else if(move==2){ //kappa3
      Bweib_update_kappa(eta3, alpha3, kappa3, frail_sm, y_sm,
                         delta_sm, delta_sm_sum, kappa3_a, kappa3_b,
                         accept_kappa3, curr_loglik_sm);
    // Rcpp::Rcout << "updated kappa3: " << kappa3 << "\n";
    } else if(move==3){ //alpha1
      //ALPHA updates (in-place)
      Bweib_update_alpha(eta1, alpha1, kappa1, frail,    y1,   delta1,
                              mhProp_alpha1_var, alpha1_a, alpha1_b,
                              accept_alpha1, curr_loglik1);
      // Rcpp::Rcout << "updated alpha1: " << alpha1 << "\n";
    } else if(move==4){ //alpha2
      Bweib_update_alpha(eta2, alpha2, kappa2, frail,    y1,   delta_cr,
                              mhProp_alpha2_var, alpha2_a, alpha2_b,
                              accept_alpha2, curr_loglik_cr);
      // Rcpp::Rcout << "updated alpha2: " << alpha2 << "\n";
    } else if(move==5){ //alpha3
      Bweib_update_alpha(eta3, alpha3, kappa3, frail_sm, y_sm, delta_sm,
                              mhProp_alpha3_var, alpha3_a, alpha3_b,
                              accept_alpha3, curr_loglik_sm);
      // Rcpp::Rcout << "updated alpha3: " << alpha3 << "\n";
    } else if(move==6){ //frailties
      //This is the update for gamma-distributed frailties
      // BweibScrSM_update_frail_gamma(frail, frail_sm, theta, eta1, eta2, eta3,
      //                               kappa1, alpha1, kappa2, alpha2, kappa3, alpha3,
      //                               y1, y_sm, sm_ind, sm_ind_long,
      //                               delta1, delta1noD, delta_cr, delta_sm, accept_frail,
      //                               curr_loglik1, curr_loglik_cr, curr_loglik_sm);

      //This is the update for lognormal-distributed frailties
      BweibScrSM_update_frail_ln(frail, frail_sm, theta, eta1, eta2, eta3,
                                    kappa1, alpha1, kappa2, alpha2, kappa3, alpha3,
                                    y1, y_sm, sm_ind, sm_ind_long,
                                    delta1, delta1noD, delta_cr, delta_sm, frail_prop_var, accept_frail,
                                    curr_loglik1, curr_loglik_cr, curr_loglik_sm);

    // Rcpp::Rcout << "updated frail: " << frail(arma::span(0,10)) << "\n";
    } else if(move==7){ //theta (frailty variance)
      //This is update for gamma-distributed frailty variance
      //BweibScrSM_update_theta_gamma(theta, frail, mhProp_theta_var, theta_a, theta_b, accept_theta);

      //This is update for lognormal distributed frailty variance
      BweibScrSM_update_theta_ln(theta, frail, theta_a, theta_b, accept_theta);

      // Rcpp::Rcout << "updated theta: " << theta << "\n";
    } else if(move==8){//beta1
      Bweib_update_betaj(beta1, eta1, alpha1, kappa1, frail, y1, delta1,
                              Xmat1, accept_beta1, curr_loglik1);
      // Rcpp::Rcout << "updated beta1: " << beta1.t() << "\n";
    } else if(move==9){//beta2
      Bweib_update_betaj(beta2, eta2, alpha2, kappa2, frail, y1, delta_cr,
                              Xmat2, accept_beta2, curr_loglik_cr);
      // Rcpp::Rcout << "updated beta2: " << beta2.t() << "\n";
    } else if(move==10){//beta3
      Bweib_update_betaj(beta3, eta3, alpha3, kappa3, frail_sm, y_sm, delta_sm,
                              Xmat3, accept_beta3, curr_loglik_sm);
      // Rcpp::Rcout << "updated beta3: " << beta3.t() << "\n";
    }

    // Storing posterior samples
    if( ( (M+1) % thin ) == 0 && (M+1) > n_burnin) {
      StoreInx = (M+1 - n_burnin)/thin;

      sample_alpha1(StoreInx - 1) = alpha1;
      sample_alpha2(StoreInx - 1) = alpha2;
      sample_alpha3(StoreInx - 1) = alpha3;
      sample_kappa1(StoreInx - 1) = kappa1;
      sample_kappa2(StoreInx - 1) = kappa2;
      sample_kappa3(StoreInx - 1) = kappa3;
      sample_theta(StoreInx - 1) = theta;
      sample_beta1.col(StoreInx - 1) = beta1;
      sample_beta2.col(StoreInx - 1) = beta2;
      sample_beta3.col(StoreInx - 1) = beta3;
      if(frail_path.size() > 0){
        sample_frail.col(StoreInx - 1) = frail;
      }
    }

    if( ( (M+1) % 10000 ) == 0){
      newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      Rcpp::Rcout << "iteration: " << M+1 << ": " << ctime(&newt) << "\n";
      Rcpp::checkUserInterrupt(); //checks if the user hit the "stop" icon to cancel running sampler.
    }
  }

  if(frail_path.size() > 0){
    // Rcpp::Rcout << "printing frailties to the address: " << frail_path << "\n" ;
    sample_frail.save(frail_path, arma::csv_ascii);
  }

  return Rcpp::List::create(
    Rcpp::Named("samples") = Rcpp::List::create(
      Rcpp::Named("alpha1") = sample_alpha1,
      Rcpp::Named("alpha2") = sample_alpha2,
      Rcpp::Named("alpha3") = sample_alpha3,
      Rcpp::Named("kappa1") = sample_kappa1,
      Rcpp::Named("kappa2") = sample_kappa2,
      Rcpp::Named("kappa3") = sample_kappa3,
      Rcpp::Named("theta") = sample_theta,
      Rcpp::Named("beta1") = sample_beta1.t(),
      Rcpp::Named("beta2") = sample_beta2.t(),
      Rcpp::Named("beta3") = sample_beta3.t()),
    Rcpp::Named("accept") = Rcpp::List::create(
      Rcpp::Named("move") = move_vec,
      Rcpp::Named("alpha1") = accept_alpha1,
      Rcpp::Named("alpha2") = accept_alpha2,
      Rcpp::Named("alpha3") = accept_alpha3,
      Rcpp::Named("kappa1") = accept_kappa1,
      Rcpp::Named("kappa2") = accept_kappa2,
      Rcpp::Named("kappa3") = accept_kappa3,
      Rcpp::Named("theta") = accept_theta,
      Rcpp::Named("beta1") = accept_beta1,
      Rcpp::Named("beta2") = accept_beta2,
      Rcpp::Named("beta3") = accept_beta3,
      Rcpp::Named("gamma") = accept_frail));

}


// [[Rcpp::export]]
void WeibSCRlogitmcmc(const arma::vec &y1, const arma::vec &y_sm,
                       const arma::uvec &delta1, const arma::uvec &delta1noD, const arma::uvec &delta_cr,
                       const arma::uvec &delta_sm, const arma::uvec &delta1D_sub,
                       const arma::mat &Xmat1, const arma::mat &Xmat2,
                       const arma::mat &Xmat3, const arma::mat &XmatD,
                       const arma::vec &hyper_vec,
                       const arma::vec &tuning_vec,
                       const arma::vec &start_vec,
                       arma::vec &sample_alpha1,
                       arma::vec &sample_alpha2,
                       arma::vec &sample_alpha3,
                       arma::vec &sample_kappa1,
                       arma::vec &sample_kappa2,
                       arma::vec &sample_kappa3,
                       arma::mat &sample_beta1,
                       arma::mat &sample_beta2,
                       arma::mat &sample_beta3,
                       arma::mat &sample_betaD,
                       arma::mat &sample_frail,
                       arma::vec &sample_theta,
                       arma::vec &accept_base,
                       arma::vec &accept_frail,
                       arma::vec &accept_beta1,
                       arma::vec &accept_beta2,
                       arma::vec &accept_beta3,
                       arma::vec &LH_marg_mean_vec,
                       arma::vec &invLH_marg_mean_vec,
                       arma::vec &sample_logLH_marg,
                       arma::mat &sample_logLHi_marg,
                       arma::vec &LH_mean_vec,
                       arma::vec &invLH_mean_vec,
                       arma::vec &sample_logLH,
                       arma::mat &sample_logLHi,
                       arma::vec &move_vec,
                       int n_burnin, int n_sample, int thin, int frail_ind,
                       int nGam_save, int nlogLHi_save,
                       const arma::vec &gh_nodes, const arma::vec &gh_weights){

  //I'm trying something new, which is explicitly passing in data as three sets
  //of vectors, corresponding with the three transitions...

  //timekeeping objects
  std::time_t newt;

  //SET CONSTANTS
  int p1 = Xmat1.n_cols;
  int p2 = Xmat2.n_cols;
  int p3 = Xmat3.n_cols;
  int pD = XmatD.n_cols;
  int n = y1.n_rows; //everyone
  int nD = XmatD.n_rows; //everyone who experiences non-terminal event
  int n_sm = y_sm.n_rows; //everyone who experiences non-terminal without immediate terminal
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'

  //INITIALIZE HYPERPARAMETERS
  double alpha1_a = hyper_vec[0];
  double alpha1_b = hyper_vec[1];
  double alpha2_a = hyper_vec[2];
  double alpha2_b = hyper_vec[3];
  double alpha3_a = hyper_vec[4];
  double alpha3_b = hyper_vec[5];
  double kappa1_a = hyper_vec[6];
  double kappa1_b = hyper_vec[7];
  double kappa2_a = hyper_vec[8];
  double kappa2_b = hyper_vec[9];
  double kappa3_a = hyper_vec[10];
  double kappa3_b = hyper_vec[11];

  //for now, I'm going to assume that the hyperparameters have inputs here
  //even if a frailty is not ultimately included. Then the indexing doesn't get jammed up.
  double theta_a  = hyper_vec[12];
  double theta_b  = hyper_vec[13];

  arma::vec m0 = hyper_vec(arma::span(14,13+pD));
  arma::vec P0diag = hyper_vec(arma::span(14+pD,13+pD+pD));

  //INITIALIZE TUNING PARAMETERS
  double mhProp_alpha1_var = tuning_vec[0];
  double mhProp_alpha2_var = tuning_vec[1];
  double mhProp_alpha3_var = tuning_vec[2];

  //again, assume these exist
  double mhProp_theta_var  = tuning_vec[3];
  double frail_prop_var = 0.3;

  // Rcpp::Rcout << "finished setting MCMC tuning params" << "\n";

  //INITIALIZE STARTING VALUES
  double kappa1=start_vec(0);
  double alpha1=start_vec(1);
  double kappa2=start_vec(2);
  double alpha2=start_vec(3);
  double kappa3=start_vec(4);
  double alpha3=start_vec(5);

  //again, assume that there's a start value for theta
  double theta=start_vec(6);

  arma::vec beta1, beta2, beta3, betaD;
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n_sm,arma::fill::zeros);
  arma::vec etaD(nD,arma::fill::zeros);
  if(p1>0){
    beta1 = start_vec(arma::span(7,6+p1));
    eta1 = Xmat1 * beta1;
  }
  if(p2>0){
    beta2 = start_vec(arma::span(7+p1,6+p1+p2));
    eta2 = Xmat2 * beta2;
  }
  if(p3>0){
    beta3 = start_vec(arma::span(7+p1+p2,6+p1+p2+p3));
    eta3 = Xmat3 * beta3;
  }
  if(pD>0){
    betaD = start_vec(arma::span(7+p1+p2+p3,6+p1+p2+p3+pD));
    etaD = XmatD * betaD;
  }

  //again, assume that these exist.
  arma::vec frail = start_vec(arma::span(7+p1+p2+p3+pD,6+p1+p2+p3+pD+n));

  //grab subset of frails corresponding with risk of immediate death
  arma::uvec delta1_ind = arma::find( delta1 );
  arma::vec frailD = frail(delta1_ind);
  //now make a lookup vector connecting the index out of n with the index out of nD
  arma::uvec delta1_ind_long = arma::uvec(n,arma::fill::zeros);
  delta1_ind_long(delta1_ind) = arma::linspace<arma::uvec>(0, nD-1,nD);

  //grab subset of frails corresponding with risk in h3 (i.e., non-term but no immediate death)
  //"sm" meaning "semi-markov"
  arma::uvec sm_ind = arma::find( delta1noD );
  arma::vec frail_sm = frail( sm_ind );
  //now make a lookup vector connecting the index out of n with the index out of n_sm
  arma::uvec sm_ind_long = arma::uvec(n,arma::fill::zeros);
  sm_ind_long(sm_ind) = arma::linspace<arma::uvec>(0, n_sm-1,n_sm);

  // Rcpp::Rcout << "cbind(delta1noD, frail): " << arma::join_horiz(arma::conv_to<arma::vec>::from(delta1noD(arma::span(0,20))),frail(arma::span(0,20))) << "\n";
  // Rcpp::Rcout << "frail_sm: " << frail_sm(arma::span(0,20)) << "\n";
  // Rcpp::Rcout << "cbind(delta1, frail): " << arma::join_horiz(arma::conv_to<arma::vec>::from(delta1(arma::span(0,20))),frail(arma::span(0,20))) << "\n";
  // Rcpp::Rcout << "frailD: " << frailD(arma::span(0,20)) << "\n";
  // Rcpp::Rcout << "nrow delta1_ind " << delta1_ind.n_rows << "\n";
  // Rcpp::Rcout << "nrow sm_ind " << sm_ind.n_rows << "\n";

  //structures to store running likelihood contributions
  arma::vec logLH_vec = arma::vec(n,arma::fill::zeros);
  arma::vec logLH_marg_vec = arma::vec(n,arma::fill::zeros);
  arma::vec logLH_temp_vec = arma::vec(n,arma::fill::zeros); //placeholder used for computing marginal

  //temporary storage of MH acceptance counts
  int accept_kappa1 = 0;
  int accept_alpha1 = 0;
  int accept_kappa2 = 0;
  int accept_alpha2 = 0;
  int accept_kappa3 = 0;
  int accept_alpha3 = 0;
  int accept_theta = 0;
  int accept_betaD = 0; //because all logit betas are updated at once

  int StoreInx=0; //index for where to store a sample, post-thinning

  //INITIALIZE WORKING VALUES
  double delta1_sum = arma::accu(delta1);
  double delta_cr_sum = arma::accu(delta_cr);
  double delta_sm_sum = arma::accu( delta_sm );
  arma::vec mean_const = XmatD.t() * (arma::conv_to<arma::vec>::from(delta1D_sub)-0.5) + P0diag % m0;
  //keep running log-likelihoods for each submodel
  double curr_loglik1 = logLikWB_uni(y1, delta1, alpha1, kappa1, eta1, frail);
  double curr_loglik_cr = logLikWB_uni(y1, delta_cr, alpha2, kappa2, eta2, frail);
  double curr_loglik_sm = logLikWB_uni(y_sm, delta_sm, alpha3, kappa3, eta3, frail_sm);

  int numUpdate = 6; //3 kappas, 3 alphas
  if(frail_ind > 0) numUpdate += 2; //update theta and frailties
  if(p1 > 0) numUpdate += 1;
  if(p2 > 0) numUpdate += 1;
  if(p3 > 0) numUpdate += 1;
  if(pD > 0) numUpdate += 1;

  // double prob_beta1 = (p1 > 0) ? (double) 1/numUpdate : 0;
  // double prob_beta2 = (p2 > 0) ? (double) 1/numUpdate : 0;
  // double prob_beta3 = (p3 > 0) ? (double) 1/numUpdate : 0;
  // double prob_betaD = (pD > 0) ? (double) 1/numUpdate : 0;
  // double prob_kappa1 = (double) 1/numUpdate;
  // double prob_kappa2 = (double) 1/numUpdate;
  // double prob_kappa3 = (double) 1/numUpdate;
  // double prob_alpha1 = (double) 1/numUpdate;
  // double prob_alpha2 = (double) 1/numUpdate;
  // double prob_alpha3 = (double) 1/numUpdate;
  // double prob_frail  = (double) 1/numUpdate;
  // double prob_theta = 1 - prob_beta1 - prob_beta2 - prob_beta3 - prob_betaD
  //   - prob_kappa1 - prob_kappa2 - prob_kappa3
  //   - prob_alpha1 - prob_alpha2 - prob_alpha3
  //   - prob_frail;

  //how about this:
  //about 0.2 is devoted to betas
  //at least 0.45 is devoted to theta + frailties
  //at least 0.35 is devoted to baseline params

  double total_betas = p1 + p2 + p3;
  if(pD > 0){
    total_betas = total_betas + 1;
  }
  //I'm going to reweight the sampling scheme to oversample theta
  double prob_beta1 = (p1 > 0) ? (double) 0.2 * p1 / total_betas : 0;
  double prob_beta2 = (p2 > 0) ? (double) 0.2 * p2 / total_betas : 0;
  double prob_beta3 = (p3 > 0) ? (double) 0.2 * p3 / total_betas : 0;
  double prob_betaD = (pD > 0) ? (double) 0.2 * 1 / total_betas : 0; //smallest bc updates in a block
  double prob_remaining = 1 - prob_beta1 - prob_beta2 - prob_beta3 - prob_betaD;

  double prob_frail = (frail_ind > 0) ? (double) 0.15 * prob_remaining : 0;
  double prob_theta = (frail_ind > 0) ? (double) 0.35 * prob_remaining : 0;
  prob_remaining = prob_remaining - prob_theta - prob_frail;

  //now, just divide up "proportions" of the remaining probability.
  double prob_kappa1 = (double) prob_remaining / 6;
  double prob_kappa2 = (double) prob_remaining / 6;
  double prob_kappa3 = (double) prob_remaining / 6;
  double prob_alpha1 = (double) prob_remaining / 6;
  double prob_alpha2 = (double) prob_remaining / 6;
  double prob_alpha3 = (double) prob_remaining / 6;
  arma::vec probs = {prob_kappa1,prob_kappa2,prob_kappa3,
                     prob_alpha1,prob_alpha2,prob_alpha3,
                     prob_frail,prob_theta,
                     prob_beta1,prob_beta2,prob_beta3,prob_betaD};
  arma::vec cumprobs = arma::cumsum(probs);

  // Rcpp::Rcout << "Sampling Probabilities:" << "\n"
  //             << "prob_beta1: \t" << prob_beta1 << "\n"
  //             << "prob_beta2: \t" << prob_beta2 << "\n"
  //             << "prob_beta3: \t" << prob_beta3 << "\n"
  //             << "prob_betaD: \t" << prob_betaD << "\n"
  //             << "prob_alpha1: \t" << prob_alpha1 << "\n"
  //             << "prob_kappa1: \t" << prob_kappa1 << "\n"
  //             << "prob_alpha2: \t" << prob_alpha2 << "\n"
  //             << "prob_kappa2: \t" << prob_kappa2 << "\n"
  //             << "prob_alpha3: \t" << prob_alpha3 << "\n"
  //             << "prob_kappa3: \t" << prob_kappa3 << "\n"
  //             << "prob_frail: \t" << prob_frail << "\n"
  //             << "prob_theta: \t" << prob_theta << "\n";


  double move_unit; //random number on the interval (0,1) deciding what to sample

  //RUN MCMC
  newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  Rcpp::Rcout << "Begin MCMC: " << ctime(&newt) << "\n";
  for(int M = 0; M < n_iter; ++M){
    // Rcpp::Rcout << "iteration: " << M << "\n";

    move_unit = R::unif_rand();

    //KAPPA updates
    if(move_unit < cumprobs[0]){ // kappa1
    move_vec(M) = 0;
      Bweib_update_kappa(eta1, alpha1, kappa1, frail, y1,
                         delta1, delta1_sum, kappa1_a, kappa1_b,
                         accept_kappa1, curr_loglik1);

      // Rcpp::Rcout << "updated kappa1: " << kappa1 << "\n";
    } else if(move_unit < cumprobs[1]){ //kappa2
      move_vec(M) = 1;
      Bweib_update_kappa(eta2, alpha2, kappa2, frail, y1,
                         delta_cr, delta_cr_sum, kappa2_a, kappa2_b,
                         accept_kappa2, curr_loglik_cr);

      // Rcpp::Rcout << "updated kappa2: " << kappa2 << "\n";
    } else if(move_unit < cumprobs[2]){ //kappa3
      move_vec(M) = 2;
      Bweib_update_kappa(eta3, alpha3, kappa3, frail_sm, y_sm,
                         delta_sm, delta_sm_sum, kappa3_a, kappa3_b,
                         accept_kappa3, curr_loglik_sm);
      // Rcpp::Rcout << "updated kappa3: " << kappa3 << "\n";
    } else if(move_unit < cumprobs[3]){ //alpha1
      move_vec(M) = 3;
      Bweib_update_alpha(eta1, alpha1, kappa1, frail,    y1,   delta1,
                         mhProp_alpha1_var, alpha1_a, alpha1_b,
                         accept_alpha1, curr_loglik1);
      // Rcpp::Rcout << "updated alpha1: " << alpha1 << "\n";
    } else if(move_unit < cumprobs[4]){ //alpha2
      move_vec(M) = 4;
      Bweib_update_alpha(eta2, alpha2, kappa2, frail,    y1,   delta_cr,
                         mhProp_alpha2_var, alpha2_a, alpha2_b,
                         accept_alpha2, curr_loglik_cr);
      // Rcpp::Rcout << "updated alpha2: " << alpha2 << "\n";
    } else if(move_unit < cumprobs[5]){ //alpha3
      move_vec(M) = 5;
      Bweib_update_alpha(eta3, alpha3, kappa3, frail_sm, y_sm, delta_sm,
                         mhProp_alpha3_var, alpha3_a, alpha3_b,
                         accept_alpha3, curr_loglik_sm);
      // Rcpp::Rcout << "updated alpha3: " << alpha3 << "\n";
    } else if(move_unit < cumprobs[6]){ //frailties
      move_vec(M) = 6;
      // Rcpp::Rcout << "updating frailties... " << "\n";

      //This is the update for gamma-distributed frailties
      // BweibScrSM_update_frail_gamma(frail, frail_sm, theta, eta1, eta2, eta3,
      //                               kappa1, alpha1, kappa2, alpha2, kappa3, alpha3,
      //                               y1, y_sm, sm_ind, sm_ind_long,
      //                               delta1, delta_cr, delta_sm, accept_frail,
      //                               curr_loglik1, curr_loglik_cr, curr_loglik_sm);

      //This is the update for lognormal-distributed frailties
      BweibScrSMlogit_update_frail_ln(frail, frail_sm, frailD, theta, eta1, eta2, eta3, etaD,
                                      kappa1, alpha1, kappa2, alpha2, kappa3, alpha3,
                                      y1, y_sm, delta1, delta1noD, delta_cr, delta_sm,
                                      sm_ind, sm_ind_long, delta1_ind_long,
                                      frail_prop_var, accept_frail,
                                      curr_loglik1, curr_loglik_cr, curr_loglik_sm);

      // Rcpp::Rcout << "updated frail: " << frail(arma::span(0,10)) << "\n";
    } else if(move_unit < cumprobs[7]){ //theta (frailty variance)
      move_vec(M) = 7;

      //This is update for gamma-distributed frailty variance
      //BweibScrSM_update_theta_gamma(theta, frail, mhProp_theta_var, theta_a, theta_b, accept_theta);

      //This is update for lognormal distributed frailty variance
      BweibScrSM_update_theta_ln(theta, frail, theta_a, theta_b, accept_theta);

      // Rcpp::Rcout << "updated theta: " << theta << "\n";
    } else if(move_unit < cumprobs[8]){ //beta1
      move_vec(M) = 8;

      Bweib_update_betaj(beta1, eta1, alpha1, kappa1, frail, y1, delta1,
                         Xmat1, accept_beta1, curr_loglik1);
      // Rcpp::Rcout << "updated beta1: " << beta1.t() << "\n";
    } else if(move_unit < cumprobs[9]){ //beta2
      move_vec(M) = 9;
      Bweib_update_betaj(beta2, eta2, alpha2, kappa2, frail, y1, delta_cr,
                         Xmat2, accept_beta2, curr_loglik_cr);
      // Rcpp::Rcout << "updated beta2: " << beta2.t() << "\n";
    } else if(move_unit < cumprobs[10]){ //beta3
      move_vec(M) = 10;
      Bweib_update_betaj(beta3, eta3, alpha3, kappa3, frail_sm, y_sm, delta_sm,
                         Xmat3, accept_beta3, curr_loglik_sm);
      // Rcpp::Rcout << "updated beta3: " << beta3.t() << "\n";
    } else { //betaD
      move_vec(M) = 11;

      Logit_update_beta_frail(betaD, etaD, XmatD, frailD, mean_const, P0diag, accept_betaD);

      // Rcpp::Rcout << "updated betaD: " << betaD.t() << "\n";
    }

    // Storing posterior samples
    if( ( (M+1) % thin ) == 0 && (M+1) > n_burnin) {
      StoreInx = (M+1 - n_burnin)/thin;

      sample_alpha1(StoreInx - 1) = alpha1;
      sample_alpha2(StoreInx - 1) = alpha2;
      sample_alpha3(StoreInx - 1) = alpha3;
      sample_kappa1(StoreInx - 1) = kappa1;
      sample_kappa2(StoreInx - 1) = kappa2;
      sample_kappa3(StoreInx - 1) = kappa3;
      sample_beta1.col(StoreInx - 1) = beta1;
      sample_beta2.col(StoreInx - 1) = beta2;
      sample_beta3.col(StoreInx - 1) = beta3;
      sample_betaD.col(StoreInx - 1) = betaD;
      if(frail_ind>0){
        sample_theta(StoreInx - 1) = theta;

        if(nGam_save>0){
          sample_frail.col(StoreInx - 1) = frail.head(nGam_save);
        }

        // Rcpp::Rcout << "iter: " << M << "\n";
        BweibScrSMlogit_logLH_marg_vec(logLH_marg_vec, eta1, eta2, eta3, etaD,
                                       kappa1, alpha1, kappa2, alpha2, kappa3, alpha3, theta,
                                       y1, y_sm, delta1, delta1noD, delta_cr, delta_sm,
                                       sm_ind_long, delta1_ind_long, gh_nodes, gh_weights);
        sample_logLH_marg(StoreInx - 1) = arma::accu(logLH_marg_vec);
        LH_marg_mean_vec = ((StoreInx - 1) * LH_marg_mean_vec + arma::exp(logLH_marg_vec)) / StoreInx;
        invLH_marg_mean_vec = ((StoreInx - 1) * invLH_marg_mean_vec + arma::exp(-logLH_marg_vec)) / StoreInx;
        if(nlogLHi_save>0){
          sample_logLHi_marg.col(StoreInx - 1) = logLH_marg_vec.head(nlogLHi_save);
        }
      }

      //deviance information (might be better to compute "marginal" version here, millar 2009)
      BweibScrSMlogit_logLH_vec(logLH_vec, frail, eta1, eta2, eta3, etaD,
                                kappa1, alpha1, kappa2, alpha2, kappa3, alpha3,
                                y1, y_sm, delta1, delta1noD, delta_cr, delta_sm,
                                sm_ind_long, delta1_ind_long);
      if(nlogLHi_save>0){
        sample_logLHi.col(StoreInx - 1) = logLH_vec.head(nlogLHi_save);
      }
      //store overall log likelihood sample
      sample_logLH(StoreInx - 1) = arma::accu(logLH_vec);
      //update running mean likelihood contributions
      LH_mean_vec = ((StoreInx - 1) * LH_mean_vec + arma::exp(logLH_vec)) / StoreInx;
      invLH_mean_vec = ((StoreInx - 1) * invLH_mean_vec + arma::exp(-logLH_vec)) / StoreInx;
    }

    if( ( (M+1) % 10000 ) == 0){
      newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      Rcpp::Rcout << "iteration: " << M+1 << ": " << ctime(&newt) << "\n";
      Rcpp::checkUserInterrupt(); //checks if the user hit the "stop" icon to cancel running sampler.
    }
  }

  accept_base(0) = accept_alpha1;
  accept_base(1) = accept_alpha2;
  accept_base(2) = accept_alpha3;
  accept_base(3) = accept_kappa1;
  accept_base(4) = accept_kappa2;
  accept_base(5) = accept_kappa3;
  accept_base(6) = accept_theta;
  accept_base(7) = accept_betaD;

  return;

}







// ADDING A COEFFICIENT IN LOGIT MODEL FOR THE FRAILTY EFFECT

void BweibScrSMlogit_update_frail_ln_coef(arma::vec &frail, arma::vec &frail_sm, const double &theta,
                                     const arma::vec &eta1, const arma::vec &eta2, const arma::vec &eta3,
                                     arma::vec &etaD, const arma::vec &betaD,
                                     const double &kappa1, const double &alpha1,
                                     const double &kappa2, const double &alpha2,
                                     const double &kappa3, const double &alpha3,
                                     const arma::vec &y1, const arma::vec &y_sm,
                                     const arma::uvec &delta1, const arma::uvec &delta1noD,
                                     const arma::uvec &delta_cr, const arma::uvec &delta_sm,
                                     const arma::uvec &sm_ind, const arma::uvec &sm_ind_long,
                                     const arma::uvec &delta1_ind_long, arma::mat &XmatD,
                                     const double& frail_prop_var, arma::vec &accept_frail,
                                     double &curr_loglik1, double &curr_loglik_cr, double &curr_loglik_sm){

  //in this formulation, etaD already contains betaD(0) * log(frail)

  int n = y1.n_rows;
  double logliki,logliki_prop,fraili_prop,logprior,logprior_prop, logR;
  for(int i = 0; i < n; i++){
    //frail is on positives, so normal is drawing the log-frailty, then we exponentiate
    fraili_prop = exp(R::rnorm(log(frail(i)), sqrt(frail_prop_var)));

    logliki = logLikWB_uni_i(y1(i), delta1(i), alpha1, kappa1, eta1(i), frail(i))
      + logLikWB_uni_i(y1(i), delta_cr(i), alpha2, kappa2, eta2(i), frail(i));
    logliki_prop = logLikWB_uni_i(y1(i), delta1(i), alpha1, kappa1, eta1(i), fraili_prop)
      + logLikWB_uni_i(y1(i), delta_cr(i), alpha2, kappa2, eta2(i), fraili_prop);
    if(delta1(i)>0){ //nonterminal event has occurred
      //I'm writing binary log-likelihood contribution as in
      //slide 22 of https://www.stat.rutgers.edu/home/pingli/papers/Logit.pdf
      logliki += -log1p( exp(etaD(delta1_ind_long(i))) ); //etaD already contains current log-frailty
      logliki_prop += -log1p( exp(etaD(delta1_ind_long(i))
                                  + betaD(0) * (log(fraili_prop) - log(frail(i))) ) );
      if(delta1noD(i)>0){ //if the non-terminal event has occurred without immediate death
        logliki +=      logLikWB_uni_i(y_sm(sm_ind_long(i)), delta_sm(sm_ind_long(i)),
                                       alpha3, kappa3, eta3(sm_ind_long(i)), frail(i));
        logliki_prop += logLikWB_uni_i(y_sm(sm_ind_long(i)), delta_sm(sm_ind_long(i)),
                                       alpha3, kappa3, eta3(sm_ind_long(i)), fraili_prop);
      } else{ //then non-terminal event has occurred followed by immediate death
        //rest of etaD cancels out from these, so just add the log-frailties
        logliki += betaD(0) * log(frail(i));
        logliki_prop += betaD(0) * log(fraili_prop);
      }
    }

    logprior = R::dnorm(log(frail(i)), 0, sqrt(theta), 1);
    logprior_prop = R::dnorm(log(fraili_prop), 0, sqrt(theta), 1);
    logR = logliki_prop - logliki + logprior_prop - logprior;
    if( log(R::unif_rand()) < logR ){
      frail(i) = fraili_prop;
      accept_frail(i) = accept_frail(i) + 1;
      if(delta1(i)>0){
        etaD(delta1_ind_long(i)) = etaD(delta1_ind_long(i)) + betaD(0) * (log(fraili_prop) - log(frail(i)));
      }
    }
  }

  // Rcpp::Rcout << "actually did all the frailty updates!" << "\n";

  XmatD.col(0) = arma::log(frail( arma::find(delta1) ));
  frail_sm = frail( sm_ind );

  //Rcpp::Rcout << "subsetted the frailties!" << "\n";

  curr_loglik1 = logLikWB_uni(y1, delta1, alpha1, kappa1, eta1, frail);
  curr_loglik_cr = logLikWB_uni(y1, delta_cr, alpha2, kappa2, eta2, frail);
  curr_loglik_sm = logLikWB_uni(y_sm, delta_sm, alpha3, kappa3, eta3, frail_sm);
}

void Logit_update_beta_frail_coef(arma::vec& beta, arma::vec& eta, const arma::mat& Xmat,
                                  const arma::uvec& delta1D_sub, int& accept_beta){

  //frailties are already contained in first column of Xmat
  int n = Xmat.n_rows;
  int p = Xmat.n_cols;
  //initialize starting values
  arma::vec omega = arma::ones(n);
  arma::vec beta_mean = arma::zeros(p);
  arma::mat beta_var = arma::eye(p,p);

  //draw auxiliary omega variables
  for(int i = 0; i < n; i++){
    omega(i) = rpg1z_devroye(eta(i)); //log-frailty with coefficient already included in ith linear predictor
  }

  // Rcpp::Rcout << "Current XmatD head: " << Xmat.rows(0,10) << "\n";

  // Rcpp::Rcout << "Beta precision: " << Xmat.t() * (Xmat.each_col() % omega) << "\n";


  //compute mean and variance vectors for beta draw
  beta_var = inv_sympd(Xmat.t() * (Xmat.each_col() % omega));
  // beta_var = inv_sympd(Xmat.t() * (Xmat.each_col() % omega) + P0); //if we had a prior P0
  beta_mean = beta_var * (Xmat.t() * (arma::conv_to<arma::vec>::from(delta1D_sub)-0.5));

  //draw beta (inplace)
  mvrnorm_inplace(beta,beta_mean,beta_var);
  //update values
  eta = Xmat * beta;
  accept_beta++;
}


// [[Rcpp::export]]
Rcpp::List WeibSCRlogitmcmc2(const arma::vec &y1, const arma::vec &y_sm,
                            const arma::uvec &delta1, const arma::uvec &delta1noD, const arma::uvec &delta_cr,
                            const arma::uvec &delta_sm, const arma::uvec &delta1D_sub,
                            const arma::mat &Xmat1, const arma::mat &Xmat2,
                            const arma::mat &Xmat3, arma::mat &XmatD,
                            const arma::vec &hyper_vec,
                            const arma::vec &tuning_vec,
                            const arma::vec &start_vec,
                            int n_burnin,
                            int n_sample,
                            int thin){

  //this version includes as the first element of betaD a coefficient for the frailties

  //timekeeping objects
  std::time_t newt;

  //SET CONSTANTS
  int p1 = Xmat1.n_cols;
  int p2 = Xmat2.n_cols;
  int p3 = Xmat3.n_cols;
  int pD = XmatD.n_cols;
  int n = y1.n_rows; //everyone
  int nD = XmatD.n_rows; //everyone who experiences non-terminal event
  int n_sm = y_sm.n_rows; //everyone who experiences non-terminal without immediate terminal
  int n_store = n_sample / thin; //tested in wrapper function that these numbers 'work'
  int n_iter = n_burnin + n_sample; //tested elsewhere that this 'fits'

  //INITIALIZE HYPERPARAMETERS
  double alpha1_a = hyper_vec[0];
  double alpha1_b = hyper_vec[1];
  double alpha2_a = hyper_vec[2];
  double alpha2_b = hyper_vec[3];
  double alpha3_a = hyper_vec[4];
  double alpha3_b = hyper_vec[5];
  double kappa1_a = hyper_vec[6];
  double kappa1_b = hyper_vec[7];
  double kappa2_a = hyper_vec[8];
  double kappa2_b = hyper_vec[9];
  double kappa3_a = hyper_vec[10];
  double kappa3_b = hyper_vec[11];
  double theta_a  = hyper_vec[12];
  double theta_b  = hyper_vec[13];

  //INITIALIZE TUNING PARAMETERS
  double mhProp_alpha1_var = tuning_vec[0];
  double mhProp_alpha2_var = tuning_vec[1];
  double mhProp_alpha3_var = tuning_vec[2];
  double mhProp_theta_var  = tuning_vec[3];

  double frail_prop_var = 0.3;

  // Rcpp::Rcout << "finished setting MCMC tuning params" << "\n";


  //INITIALIZE STARTING VALUES
  double kappa1=start_vec(0);
  double alpha1=start_vec(1);
  double kappa2=start_vec(2);
  double alpha2=start_vec(3);
  double kappa3=start_vec(4);
  double alpha3=start_vec(5);
  double theta=start_vec(6);
  arma::vec beta1, beta2, beta3, betaD;
  arma::vec eta1(n,arma::fill::zeros);
  arma::vec eta2(n,arma::fill::zeros);
  arma::vec eta3(n_sm,arma::fill::zeros);
  arma::vec etaD(nD,arma::fill::zeros);
  if(p1>0){
    beta1 = start_vec(arma::span(7,6+p1));
    eta1 = Xmat1 * beta1;
  }
  if(p2>0){
    beta2 = start_vec(arma::span(7+p1,6+p1+p2));
    eta2 = Xmat2 * beta2;
  }
  if(p3>0){
    beta3 = start_vec(arma::span(7+p1+p2,6+p1+p2+p3));
    eta3 = Xmat3 * beta3;
  }

  arma::vec frail = start_vec(arma::span(7+p1+p2+p3+pD,6+p1+p2+p3+pD+n));
  //test model with fixed frailties, verify that that is working ok!
  // arma::vec frail = arma::ones(n);
  // theta=0.5;

  //grab subset of frails corresponding with risk of immediate death
  arma::uvec delta1_ind = arma::find( delta1 );
  XmatD.col(0) = arma::log(frail(delta1_ind)); //first column of XmatD holds log-frailties
  //now make a lookup vector connecting the index out of n with the index out of nD
  arma::uvec delta1_ind_long = arma::uvec(n,arma::fill::zeros);
  delta1_ind_long(delta1_ind) = arma::linspace<arma::uvec>(0, nD-1,nD);

  //grab subset of frails corresponding with risk in h3 (i.e., non-term but no immediate death)
  //"sm" meaning "semi-markov"
  arma::uvec sm_ind = arma::find( delta1noD );
  arma::vec frail_sm = frail( sm_ind );
  //now make a lookup vector connecting the index out of n with the index out of n_sm
  arma::uvec sm_ind_long = arma::uvec(n,arma::fill::zeros);
  sm_ind_long(sm_ind) = arma::linspace<arma::uvec>(0, n_sm-1,n_sm);

  // Rcpp::Rcout << "cbind(delta1noD, frail): " << arma::join_horiz(arma::conv_to<arma::vec>::from(delta1noD(arma::span(0,20))),frail(arma::span(0,20))) << "\n";
  // Rcpp::Rcout << "frail_sm: " << frail_sm(arma::span(0,20)) << "\n";
  // Rcpp::Rcout << "cbind(delta1, frail): " << arma::join_horiz(arma::conv_to<arma::vec>::from(delta1(arma::span(0,20))),frail(arma::span(0,20))) << "\n";
  // Rcpp::Rcout << "frailD: " << frailD(arma::span(0,20)) << "\n";
  // Rcpp::Rcout << "nrow delta1_ind " << delta1_ind.n_rows << "\n";
  // Rcpp::Rcout << "nrow sm_ind " << sm_ind.n_rows << "\n";

  //now that we've initialized frailties, we can compute the logit linear predictor
  if(pD>0){
    betaD = start_vec(arma::span(7+p1+p2+p3,6+p1+p2+p3+pD));
    etaD = XmatD * betaD;
  }

  //CREATE STORAGE VECTORS/MATRICES FOR SAMPLING
  arma::vec sample_kappa1 = arma::vec(n_store,arma::fill::zeros);
  arma::vec sample_alpha1 = arma::vec(n_store,arma::fill::zeros);
  arma::mat sample_beta1 = arma::mat(p1,n_store,arma::fill::zeros);
  arma::vec sample_kappa2 = arma::vec(n_store,arma::fill::zeros);
  arma::vec sample_alpha2 = arma::vec(n_store,arma::fill::zeros);
  arma::mat sample_beta2 = arma::mat(p2,n_store,arma::fill::zeros);
  arma::vec sample_kappa3 = arma::vec(n_store,arma::fill::zeros);
  arma::vec sample_alpha3 = arma::vec(n_store,arma::fill::zeros);
  arma::mat sample_beta3 = arma::mat(p3,n_store,arma::fill::zeros);
  arma::mat sample_frail = arma::mat(n,n_store,arma::fill::zeros);
  arma::vec sample_theta = arma::vec(n_store,arma::fill::zeros);
  arma::mat sample_betaD = arma::mat(pD,n_store,arma::fill::zeros);

  int accept_kappa1 = 0;
  int accept_alpha1 = 0;
  int accept_kappa2 = 0;
  int accept_alpha2 = 0;
  int accept_kappa3 = 0;
  int accept_alpha3 = 0;
  int accept_theta = 0;
  int accept_betaD = 0; //because all logit betas are updated at once
  arma::vec accept_beta1 = arma::vec(p1,arma::fill::zeros);
  arma::vec accept_beta2 = arma::vec(p2,arma::fill::zeros);
  arma::vec accept_beta3 = arma::vec(p3,arma::fill::zeros);

  arma::vec accept_frail = arma::vec(n,arma::fill::zeros);
  //int accept_frail = 0;

  int StoreInx=0; //index for where to store a sample, post-thinning

  //INITIALIZE WORKING VALUES
  double delta1_sum = arma::accu(delta1);
  double delta_cr_sum = arma::accu(delta_cr);
  double delta_sm_sum = arma::accu( delta_sm );

  //keep running log-likelihoods for each submodel
  double curr_loglik1 = logLikWB_uni(y1, delta1, alpha1, kappa1, eta1, frail);
  double curr_loglik_cr = logLikWB_uni(y1, delta_cr, alpha2, kappa2, eta2, frail);
  double curr_loglik_sm = logLikWB_uni(y_sm, delta_sm, alpha3, kappa3, eta3, frail_sm);

  int numUpdate = 8; //3 kappas, 3 alphas, theta and frailties
  if(p1 > 0) numUpdate += 1;
  if(p2 > 0) numUpdate += 1;
  if(p3 > 0) numUpdate += 1;
  if(pD > 0) numUpdate += 1;
  double prob_beta1 = (p1 > 0) ? (double) 1/numUpdate : 0;
  double prob_beta2 = (p2 > 0) ? (double) 1/numUpdate : 0;
  double prob_beta3 = (p3 > 0) ? (double) 1/numUpdate : 0;
  double prob_betaD = (pD > 0) ? (double) 1/numUpdate : 0;
  double prob_kappa1 = (double) 1/numUpdate;
  double prob_kappa2 = (double) 1/numUpdate;
  double prob_kappa3 = (double) 1/numUpdate;
  double prob_alpha1 = (double) 1/numUpdate;
  double prob_alpha2 = (double) 1/numUpdate;
  double prob_alpha3 = (double) 1/numUpdate;
  double prob_frail  = (double) 1/numUpdate;
  double prob_theta = 1 - prob_beta1 - prob_beta2 - prob_beta3 - prob_betaD
    - prob_kappa1 - prob_kappa2 - prob_kappa3
    - prob_alpha1 - prob_alpha2 - prob_alpha3
    - prob_frail;

    int move; //index for which parameter to update

    //RUN MCMC
    newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    Rcpp::Rcout << "Begin MCMC: " << ctime(&newt) << "\n";

    // Rcpp::Rcout << "begin mcmc" << "\n";
    arma::uvec move_vec = arma::uvec(n_iter,arma::fill::zeros);
    for(int M = 0; M < n_iter; ++M){
      // Rcpp::Rcout << "iteration: " << M << "\n";

      move = (int) R::runif(0, numUpdate);
      move_vec(M) = move;

      //KAPPA updates
      if(move==0){ //kappa1
        Bweib_update_kappa(eta1, alpha1, kappa1, frail, y1,
                           delta1, delta1_sum, kappa1_a, kappa1_b,
                           accept_kappa1, curr_loglik1);

        // Rcpp::Rcout << "updated kappa1: " << kappa1 << "\n";
      } else if(move==1){ //kappa2
        Bweib_update_kappa(eta2, alpha2, kappa2, frail, y1,
                           delta_cr, delta_cr_sum, kappa2_a, kappa2_b,
                           accept_kappa2, curr_loglik_cr);

        // Rcpp::Rcout << "updated kappa2: " << kappa2 << "\n";
      } else if(move==2){ //kappa3
        Bweib_update_kappa(eta3, alpha3, kappa3, frail_sm, y_sm,
                           delta_sm, delta_sm_sum, kappa3_a, kappa3_b,
                           accept_kappa3, curr_loglik_sm);
        // Rcpp::Rcout << "updated kappa3: " << kappa3 << "\n";
      } else if(move==3){ //alpha1
        //ALPHA updates (in-place)
        Bweib_update_alpha(eta1, alpha1, kappa1, frail,    y1,   delta1,
                           mhProp_alpha1_var, alpha1_a, alpha1_b,
                           accept_alpha1, curr_loglik1);
        // Rcpp::Rcout << "updated alpha1: " << alpha1 << "\n";
      } else if(move==4){ //alpha2
        Bweib_update_alpha(eta2, alpha2, kappa2, frail,    y1,   delta_cr,
                           mhProp_alpha2_var, alpha2_a, alpha2_b,
                           accept_alpha2, curr_loglik_cr);
        // Rcpp::Rcout << "updated alpha2: " << alpha2 << "\n";
      } else if(move==5){ //alpha3
        Bweib_update_alpha(eta3, alpha3, kappa3, frail_sm, y_sm, delta_sm,
                           mhProp_alpha3_var, alpha3_a, alpha3_b,
                           accept_alpha3, curr_loglik_sm);
        // Rcpp::Rcout << "updated alpha3: " << alpha3 << "\n";
      } else if(move==6){ //frailties
        // Rcpp::Rcout << "updating frailties... " << "\n";

        //This is the update for gamma-distributed frailties
        // BweibScrSM_update_frail_gamma(frail, frail_sm, theta, eta1, eta2, eta3,
        //                               kappa1, alpha1, kappa2, alpha2, kappa3, alpha3,
        //                               y1, y_sm, sm_ind, sm_ind_long,
        //                               delta1, delta_cr, delta_sm, accept_frail,
        //                               curr_loglik1, curr_loglik_cr, curr_loglik_sm);

        //This is the update for lognormal-distributed frailties
        BweibScrSMlogit_update_frail_ln_coef(frail, frail_sm, theta, eta1, eta2, eta3, etaD, betaD,
                                        kappa1, alpha1, kappa2, alpha2, kappa3, alpha3,
                                        y1, y_sm, delta1, delta1noD, delta_cr, delta_sm,
                                        sm_ind, sm_ind_long, delta1_ind_long, XmatD,
                                        frail_prop_var, accept_frail,
                                        curr_loglik1, curr_loglik_cr, curr_loglik_sm);

        // Rcpp::Rcout << "updated frail: " << frail(arma::span(0,10)) << "\n";
      } else if(move==7){ //theta (frailty variance)
        //This is update for gamma-distributed frailty variance
        //BweibScrSM_update_theta_gamma(theta, frail, mhProp_theta_var, theta_a, theta_b, accept_theta);

        //This is update for lognormal distributed frailty variance
        BweibScrSM_update_theta_ln(theta, frail, theta_a, theta_b, accept_theta);

        // Rcpp::Rcout << "updated theta: " << theta << "\n";
      } else if(move==8){//beta1
        Bweib_update_betaj(beta1, eta1, alpha1, kappa1, frail, y1, delta1,
                           Xmat1, accept_beta1, curr_loglik1);
        // Rcpp::Rcout << "updated beta1: " << beta1.t() << "\n";
      } else if(move==9){//beta2
        Bweib_update_betaj(beta2, eta2, alpha2, kappa2, frail, y1, delta_cr,
                           Xmat2, accept_beta2, curr_loglik_cr);
        // Rcpp::Rcout << "updated beta2: " << beta2.t() << "\n";
      } else if(move==10){//beta3
        Bweib_update_betaj(beta3, eta3, alpha3, kappa3, frail_sm, y_sm, delta_sm,
                           Xmat3, accept_beta3, curr_loglik_sm);
        // Rcpp::Rcout << "updated beta3: " << beta3.t() << "\n";
      } else if(move==11){//betaD
        // Rcpp::Rcout << "updating betaD... " << "\n";

        Logit_update_beta_frail_coef(betaD, etaD, XmatD, delta1D_sub, accept_betaD);

        // Rcpp::Rcout << "updated betaD: " << betaD.t() << "\n";
      }

      // Storing posterior samples
      if( ( (M+1) % thin ) == 0 && (M+1) > n_burnin) {
        StoreInx = (M+1 - n_burnin)/thin;

        sample_alpha1(StoreInx - 1) = alpha1;
        sample_alpha2(StoreInx - 1) = alpha2;
        sample_alpha3(StoreInx - 1) = alpha3;
        sample_kappa1(StoreInx - 1) = kappa1;
        sample_kappa2(StoreInx - 1) = kappa2;
        sample_kappa3(StoreInx - 1) = kappa3;
        sample_theta(StoreInx - 1) = theta;
        sample_beta1.col(StoreInx - 1) = beta1;
        sample_beta2.col(StoreInx - 1) = beta2;
        sample_beta3.col(StoreInx - 1) = beta3;
        sample_betaD.col(StoreInx - 1) = betaD;
        sample_frail.col(StoreInx - 1) = frail;
      }

      if( ( (M+1) % 10000 ) == 0){
        newt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        Rcpp::Rcout << "iteration: " << M+1 << ": " << ctime(&newt) << "\n";
        Rcpp::checkUserInterrupt(); //checks if the user hit the "stop" icon to cancel running sampler.
      }
    }

    return Rcpp::List::create(
      Rcpp::Named("samples") = Rcpp::List::create(
        Rcpp::Named("alpha1") = sample_alpha1,
        Rcpp::Named("alpha2") = sample_alpha2,
        Rcpp::Named("alpha3") = sample_alpha3,
        Rcpp::Named("kappa1") = sample_kappa1,
        Rcpp::Named("kappa2") = sample_kappa2,
        Rcpp::Named("kappa3") = sample_kappa3,
        Rcpp::Named("theta") = sample_theta,
        Rcpp::Named("beta1") = sample_beta1.t(),
        Rcpp::Named("beta2") = sample_beta2.t(),
        Rcpp::Named("beta3") = sample_beta3.t(),
        Rcpp::Named("betaD") = sample_betaD.t(),
        Rcpp::Named("gamma") = sample_frail.t()),
        Rcpp::Named("accept") = Rcpp::List::create(
          Rcpp::Named("move") = move_vec,
          Rcpp::Named("alpha1") = accept_alpha1,
          Rcpp::Named("alpha2") = accept_alpha2,
          Rcpp::Named("alpha3") = accept_alpha3,
          Rcpp::Named("kappa1") = accept_kappa1,
          Rcpp::Named("kappa2") = accept_kappa2,
          Rcpp::Named("kappa3") = accept_kappa3,
          Rcpp::Named("theta") = accept_theta,
          Rcpp::Named("beta1") = accept_beta1,
          Rcpp::Named("beta2") = accept_beta2,
          Rcpp::Named("beta3") = accept_beta3,
          Rcpp::Named("betaD") = accept_betaD,
          Rcpp::Named("gamma") = accept_frail));

}



