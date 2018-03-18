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
std::vector<double> expectedXi(arma::vec mu, arma::vec mi, arma::mat psi, arma::mat Ci){
	int Kt = psi.n_rows;
  std::vector<double> xi;
  arma::mat psi_sq;
  double xi_sq;
  
  // psi_sq = trans(psi) * psi;
  // xi_sq = as_scalar(trans(mu) * mu) + 
  // 	2.0 * as_scalar(trans(mu) * psi * mi) +
  // 	trace(psi_sq * Ci) + as_scalar(trans(mi) * psi_sq * mi );
  // 
  // xi.push_back(sqrt(xi_sq));
  
  
  for(int k = 0; k < Kt; k++){
  	psi_sq = trans(psi.row(k)) * psi.row(k);
  	//mu_k = as_scalar(mu[k]);
  	xi_sq = mu[k] * mu[k] + 
  		2.0 * as_scalar(mu[k] * psi.row(k) * mi) +
  		trace(psi_sq * Ci) + as_scalar(trans(mi) * psi_sq * mi );
  	
  	xi.push_back(sqrt(xi_sq));
  } // end for loop
  
  return(xi);
} 
