#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
//' Apply lambda transformation of variational parameter.
//'
//' Simple function for use within other C++ functions.
//'
//' @param x The value to which you apply the function
//' @return A numeric value that has been transformed.
// [[Rcpp::export]]
double lambdaF(double x){
  double y = (0.5- pow(1.0+exp(-x), -1))/2.0/x;
  return y;
}


//' Calculate quadratic form of spline basis functions for the current subject.
//'
//' Calculations quadratic form of theta with diagonalized variational parameter in the center.
//'
//' @param xi vector of variational parameters for the current subject.
//' @param theta spline basis functions for the current subject.
//' @return A matrix of the quadratic form of theta for the current subject.
// [[Rcpp::export]]
arma::mat squareTheta(arma::vec xi, arma::mat theta){
  xi.transform(lambdaF); 
  arma::mat temp_xi_theta = diagmat(xi) * theta;
  arma::mat theta_quad =  trans(theta) * temp_xi_theta;
  
  return(theta_quad);
}
