#ifndef PG_utilities_H
#define PG_utilities_H

double rpg1z_devroye(double z);
void mvrnorm_inplace(arma::vec& x, const arma::vec& mu, const arma::mat& sigma);

#endif
