#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//These computations are adapted from the paper Polson, Scott & Windle (2013)

const double __PI = 3.141592653589793238462643383279502884197;

double rtigauss(double z, double t){
  // drawing from inverse-gaussian random variable truncated to be below t
  double X = t + 1.0;
  if (1/t > z) { // mu > t
    double alpha = 0.0;
    while (R::unif_rand() > alpha) {
      double E1 = R::exp_rand(); double E2 = R::exp_rand();
      while ( E1*E1 > 2 * E2 / t) {
        E1 = R::exp_rand(); E2 = R::exp_rand();
      }
      X = t / ( (1 + E1 * t) * (1 + E1 * t) );
      alpha = exp(-0.5 * z*z * X);
    }
  } else {
    double mu = 1.0 / z;
    while (X > t) {
      double Y = R::norm_rand(); Y *= Y;
      X = mu + 0.5*mu*mu*Y - 0.5*mu * sqrt(4 * mu * Y + mu*mu*Y*Y);
      if (R::unif_rand() > mu/(mu + X)){
        X = mu*mu / X;
      }
    }
  }
  return X;
}

double a_coef(int n, double x, double t) {
  //helper function to compute piecewise coefficients for PG generation
  double temp = (n + 0.5) * __PI;
  double out = 0;
  if (x > t) {
    out = temp * exp( -0.5 * temp*temp * x );
  } else if (x > 0) {
    out = exp(log(temp) - 1.5 * (log(0.5 * __PI)  + log(x))
                - 2.0 * (n+0.5)*(n+0.5) / x);
  }
  return out;
}

// [[Rcpp::export]]
double rpg1z_devroye(double z){
  // drawing from a PG(1,z) polya gamma random variable
  z = z * 0.5;
  double t = 0.64;
  double K = __PI*__PI*0.125 + z*z*0.5;
  double X = 0.0;
  double S = 1.0;
  double Y = 0.0;

  double x0 = log(K) + K * t;
  double xb = x0 - z + R::pnorm(sqrt(1.0/t)*(t*z-1), 0,1,true,true);
  double xa = x0 + z + R::pnorm(sqrt(1.0/t)*(t*z+1)*-1.0, 0,1,true,true);
  double qdivp = 4 / __PI * ( exp(xb) + exp(xa) );



  while (true) {
    // if (unif() < p/(p+q))
    if ( R::unif_rand() < (1.0 / (1.0 + qdivp)) ){
      X = t + R::exp_rand() / K; //exponential r.v. truncated at t
    } else{
      //generate truncated inverse gaussian
      X = rtigauss(z,t);
    }

    S = a_coef(0, X, t);
    Y = R::unif_rand() * S;
    int n = 0;
    bool go = true;

    // Cap the number of iterations?
    while (go) {
      ++n;
      if(n == 1000){
        return -999; //well this is not a good solution but I'll let it be for now.
      } else if (n%2==1) {
        S = S - a_coef(n, X, t);
        if ( Y<=S ) return 0.25 * X;
      } else {
        S = S + a_coef(n, X, t);
        if ( Y>S ) go = false;
      }
    }
    // Need Y <= S in event that Y = S, e.g. when X = 0.
  }
}

// generate multivariate normal (mu, sigma) and replacethe
void mvrnorm_inplace(arma::vec& x, const arma::vec& mu, const arma::mat& sigma) {
  int p = x.n_rows;
  for (int i=0; i<p; i++){
    x(i) = R::norm_rand();
  }
  x = mu + (((arma::chol(sigma)).t()) * x);
  return;
}

