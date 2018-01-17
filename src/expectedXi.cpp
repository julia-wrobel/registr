// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;
//' Estimate variational parameter for the current subject.
//'
//' Function calculates value of variational parameter using maximum likelihood.
//'
//' @param theta spline basis functions for the current subject.
//' @param mu vector of spline coefficients for the population mean.
//' @param mi vector of expected mean scores for the current subject.
//' @param psi matrix of spline coefficients for the principal component basis functions.
//' @param Ci expected covariance matrix of scores for the current subject.
//' @return A vector of variational parameters for the current subject.
// [[Rcpp::export]]
std::vector<double> expectedXi(arma::mat theta, arma::vec mu, arma::vec mi, arma::mat psi, arma::mat Ci){
  int Di = theta.n_rows;
  std::vector<double> xi;
  arma::mat theta_sq, theta_sq_psi;
  double xi_sq;
  
  for(int t = 0; t < Di; t++){
    theta_sq = trans(theta.row(t)) * theta.row(t);
    theta_sq_psi = trans(psi) * theta_sq * psi;
    xi_sq = as_scalar(trans(mu) * theta_sq * mu) + 
      2.0 * as_scalar(trans(mu) * theta_sq * psi * mi) +
      trace(theta_sq_psi * Ci) + as_scalar(trans(mi) * theta_sq_psi * mi );
    
    xi.push_back(sqrt(xi_sq));
  } // end for loop
  
  return(xi);
} 
