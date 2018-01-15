// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Calculate expected score and score variance for the current subject.
//'
//' Calculations derived using maximum likelihood estimation.
//'
//' @param Y vector of observations for the current subject.
//' @param mu vector of spline coefficients for the population mean.
//' @param psi matrix of spline coefficients for the principal component basis functions.
//' @param theta spline basis functions for the current subject.
//' @param theta_quad quadratic form of theta for the current subject.
//' @return A list with expected score mean and variance for the current subject.
// [[Rcpp::export]]
List expectedScores(arma::vec Y, arma::vec mu, arma::mat psi, arma::mat theta, arma::mat theta_quad){
  int npc = psi.n_cols;
  arma::mat Ci_inner = arma::eye(npc, npc) - 2 * trans(psi) * theta_quad * psi;
  arma::mat Ci = inv(Ci_inner);
  
  arma::mat mi_inner = trans(Y - 0.5) * theta * psi + 2 * trans(mu) * theta_quad * psi;
  arma::mat mi = Ci * trans(mi_inner);
  
  List result;
  result["mi"] = mi;
  result["Ci"] = Ci;
  return result;
}
