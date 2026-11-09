// mm_spline.cpp - Spline values and derivatives for MM registration.
#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

inline double inv_logit_cpp(double x) {
  if (x >= 0.0) {
    double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  } else {
    double z = std::exp(x);
    return z / (1.0 + z);
  }
}

int find_span_cpp(double u, const arma::vec& U, int p) {
  int n = static_cast<int>(U.n_elem) - p - 2;
  if (u <= U[p]) {
    return p;
  }
  if (u >= U[n + 1]) {
    return n;
  }

  int low = p;
  int high = n + 1;
  int mid = (low + high) / 2;

  while (u < U[mid] || u >= U[mid + 1]) {
    if (u < U[mid]) {
      high = mid;
    } else {
      low = mid;
    }
    mid = (low + high) / 2;
  }
  return mid;
}

arma::mat ders_basis_funs_cpp(int i, double u, int p, int nDer, const arma::vec& U) {
  arma::mat ndu(p + 1, p + 1, fill::zeros);
  arma::mat ders(nDer + 1, p + 1, fill::zeros);
  arma::vec left(p + 1, fill::zeros);
  arma::vec right(p + 1, fill::zeros);
  ndu(0, 0) = 1.0;
  for (int j = 1; j <= p; ++j) {
    left(j) = u - U[i + 1 - j];
    right(j) = U[i + j] - u;
    double saved = 0.0;
    for (int r = 0; r < j; ++r) {
      ndu(j, r) = right(r + 1) + left(j - r);
      double temp = ndu(r, j - 1) / ndu(j, r);
      ndu(r, j) = saved + right(r + 1) * temp;
      saved = left(j - r) * temp;
    }
    ndu(j, j) = saved;
  }
  for (int j = 0; j <= p; ++j) {
    ders(0, j) = ndu(j, p);
  }

  arma::mat a(2, p + 1, fill::zeros);

  for (int r = 0; r <= p; ++r) {
    int s1 = 0;
    int s2 = 1;
    a(0, 0) = 1.0;

    for (int k = 1; k <= nDer; ++k) {
      double d = 0.0;
      int rk = r - k;
      int pk = p - k;
      int j1, j2;

      if (r >= k) {
        a(s2, 0) = a(s1, 0) / ndu(pk + 1, rk);
        d = a(s2, 0) * ndu(rk, pk);
      }

      if (rk >= -1) {
        j1 = 1;
      } else {
        j1 = -rk;
      }

      if (r - 1 <= pk) {
        j2 = k - 1;
      } else {
        j2 = p - r;
      }

      for (int j = j1; j <= j2; ++j) {
        a(s2, j) = (a(s1, j) - a(s1, j - 1)) / ndu(pk + 1, rk + j);
        d += a(s2, j) * ndu(rk + j, pk);
      }

      if (r <= pk) {
        a(s2, k) = -a(s1, k - 1) / ndu(pk + 1, r);
        d += a(s2, k) * ndu(r, pk);
      }

      ders(k, r) = d;
      std::swap(s1, s2);
    }
  }

  int r = p;
  for (int k = 1; k <= nDer; ++k) {
    for (int j = 0; j <= p; ++j) {
      ders(k, j) *= r;
    }
    r *= (p - k);
  }

  return ders;
}

//' Evaluate template quantities for MM registration.
//'
//' Computes the linear predictor, its derivative, and the quantities
//' needed for the quadratic MM update.
//'
//' @param t_vec Candidate registered time points.
//' @param knots_full Complete B-spline knot sequence for the template.
//' @param alpha Template spline coefficients.
//' @param ord B-spline order; default is 4 for cubic splines.
//' @param eps Lower bound for nonzero derivative magnitudes.
//' @return A list containing eta, deta, mu, w, and W2.
//' @author Linlin Wu (\email{linlin.wu@emory.edu}),
//' Jeff Goldsmith (\email{ajg2202@cumc.columbia.edu}),
//' Julia Wrobel (\email{julia.wrobel@emory.edu})
// [[Rcpp::export]]
 Rcpp::List eval_mm_quantities_cpp(const arma::vec& t_vec,
                                   const arma::vec& knots_full,
                                   const arma::vec& alpha,
                                   const int ord = 4,
                                   const double eps = 1e-8) {
   int p = ord - 1;
   int n_basis = static_cast<int>(knots_full.n_elem) - ord;

   if (alpha.n_elem != static_cast<unsigned>(n_basis)) {
     stop("length(alpha) must equal length(knots_full) - ord");
   }
   arma::vec eta(t_vec.n_elem, fill::zeros);
   arma::vec deta(t_vec.n_elem, fill::zeros);
   arma::vec mu(t_vec.n_elem, fill::zeros);
   arma::vec w(t_vec.n_elem, fill::zeros);
   arma::vec W2(t_vec.n_elem, fill::zeros);

   double left_bd = knots_full[p];
   double right_bd = knots_full[knots_full.n_elem - p - 1];

   for (arma::uword idx = 0; idx < t_vec.n_elem; ++idx) {
     double u = t_vec[idx];
     if (u < left_bd) u = left_bd;
     if (u > right_bd) u = right_bd;
     int span = find_span_cpp(u, knots_full, p);
     arma::mat ders = ders_basis_funs_cpp(span, u, p, 1, knots_full);
     int first_basis = span - p;
     double eta_i = 0.0;
     double deta_i = 0.0;
     for (int j = 0; j <= p; ++j) {
       double a_j = alpha[first_basis + j];
       eta_i += ders(0, j) * a_j;
       deta_i += ders(1, j) * a_j;
     }
     double mu_i = inv_logit_cpp(eta_i);
     double abs_d = std::abs(deta_i);
     double w_i = 0.0;
     if (deta_i > 0.0) {
       w_i = std::max(abs_d, eps);
     } else if (deta_i < 0.0) {
       w_i = -std::max(abs_d, eps);
     } else {
       w_i = 0.0;
     }

     eta[idx] = eta_i;
     deta[idx] = deta_i;
     mu[idx] = mu_i;
     w[idx] = w_i;
     W2[idx] = w_i * w_i;
   }

   return Rcpp::List::create(
     Rcpp::Named("eta") = eta,
     Rcpp::Named("deta") = deta,
     Rcpp::Named("mu") = mu,
     Rcpp::Named("w") = w,
     Rcpp::Named("W2") = W2
   );
 }
