// mm_loss.cpp - Bernoulli loss for MM registration.
#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

static inline int find_span_local(double u, const arma::vec& U, int p) {
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

static inline arma::vec basis_funs_local(int i, double u, int p, const arma::vec& U) {
  arma::vec N(p + 1, fill::zeros);
  arma::vec left(p + 1, fill::zeros);
  arma::vec right(p + 1, fill::zeros);

  N[0] = 1.0;

  for (int j = 1; j <= p; ++j) {
    left[j] = u - U[i + 1 - j];
    right[j] = U[i + j] - u;

    double saved = 0.0;
    for (int r = 0; r < j; ++r) {
      double denom = right[r + 1] + left[j - r];
      double temp = 0.0;
      if (denom != 0.0) {
        temp = N[r] / denom;
      }
      N[r] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    N[j] = saved;
  }

  return N;
}

static inline double bernoulli_nll_from_eta(double y, double eta) {
  return std::max(eta, 0.0) - y * eta + std::log1p(std::exp(-std::abs(eta)));
}

//' Calculate Bernoulli loss for a candidate inverse warping function.
//'
//' Evaluates the subject-specific template at the warped times, with
//' an optional penalty for incomplete endpoint constraints.
//'
//' @param Y Binary observations for one subject.
//' @param Theta_h Warping B-spline design matrix.
//' @param mean_coefs Subject-specific template spline coefficients.
//' @param template_knots_full Complete knot sequence for the template.
//' @param beta_inner Candidate warping coefficients, with endpoint entries
//'   determined by inc_mode.
//' @param t_min_curve Lower endpoint of the subject's time range.
//' @param t_max_curve Upper endpoint of the subject's time range.
//' @param inc_mode Endpoint mode: 0 = both fixed, 1 = leading free,
//'   2 = trailing free, 3 = both free.
//' @param lambda_inc Penalty weight for free endpoints.
//' @return Bernoulli negative log-likelihood plus any endpoint penalty.
//' @author Linlin Wu (\email{linlin.wu@emory.edu}),
//' Jeff Goldsmith (\email{ajg2202@cumc.columbia.edu}),
//' Julia Wrobel (\email{julia.wrobel@emory.edu})
// [[Rcpp::export]]
 double loss_h_binomial_nonparametric_cpp(
     const arma::vec& Y,
     const arma::mat& Theta_h,
     const arma::vec& mean_coefs,
     const arma::vec& template_knots_full,
     const arma::vec& beta_inner,
     const double t_min_curve,
     const double t_max_curve,
     const int inc_mode = 0,
     const double lambda_inc = 0.0
 ) {
   const int n = static_cast<int>(Y.n_elem);

   if (static_cast<int>(Theta_h.n_rows) != n) {
     stop("nrow(Theta_h) must equal length(Y).");
   }

   const int Kh = static_cast<int>(Theta_h.n_cols);

   arma::vec beta(Kh, fill::zeros);

   if (inc_mode == 0) {
     if (static_cast<int>(beta_inner.n_elem) != Kh - 2) {
       stop("For inc_mode = 0, length(beta_inner) must be Kh - 2.");
     }
     beta[0] = t_min_curve;
     beta.subvec(1, Kh - 2) = beta_inner;
     beta[Kh - 1] = t_max_curve;

   } else if (inc_mode == 1) {
     if (static_cast<int>(beta_inner.n_elem) != Kh - 1) {
       stop("For inc_mode = 1, length(beta_inner) must be Kh - 1.");
     }
     beta.subvec(0, Kh - 2) = beta_inner;
     beta[Kh - 1] = t_max_curve;

   } else if (inc_mode == 2) {
     if (static_cast<int>(beta_inner.n_elem) != Kh - 1) {
       stop("For inc_mode = 2, length(beta_inner) must be Kh - 1.");
     }
     beta[0] = t_min_curve;
     beta.subvec(1, Kh - 1) = beta_inner;

   } else if (inc_mode == 3) {
     if (static_cast<int>(beta_inner.n_elem) != Kh) {
       stop("For inc_mode = 3, length(beta_inner) must be Kh.");
     }
     beta = beta_inner;

   } else {
     stop("inc_mode must be one of 0,1,2,3.");
   }

   const int ord = static_cast<int>(template_knots_full.n_elem) -
     static_cast<int>(mean_coefs.n_elem);

   if (ord <= 0) {
     stop("Invalid template_knots_full / mean_coefs dimensions.");
   }

   const int p = ord - 1;

   const double left_bd  = template_knots_full[p];
   const double right_bd = template_knots_full[template_knots_full.n_elem - p - 1];

   double loss = 0.0;

   double first_hinv = NA_REAL;
   double last_hinv  = NA_REAL;

   for (int i = 0; i < n; ++i) {

     double u = arma::as_scalar(Theta_h.row(i) * beta);

     if (i == 0) first_hinv = u;
     if (i == n - 1) last_hinv = u;

     if (u < left_bd)  u = left_bd;
     if (u > right_bd) u = right_bd;

     int span = find_span_local(u, template_knots_full, p);
     arma::vec N = basis_funs_local(span, u, p, template_knots_full);

     int first_basis = span - p;

     double eta = 0.0;
     for (int j = 0; j <= p; ++j) {
       eta += N[j] * mean_coefs[first_basis + j];
     }

     loss += bernoulli_nll_from_eta(Y[i], eta);
   }

   double pen_term = 0.0;

   if ((inc_mode != 0) && (lambda_inc != 0.0)) {
     double pen_term_raw = 0.0;
     if (inc_mode == 1) {
       pen_term_raw = std::pow(first_hinv - t_min_curve, 2.0);
     } else if (inc_mode == 2) {
       pen_term_raw = std::pow(last_hinv - t_max_curve, 2.0);
     } else if (inc_mode == 3) {
       pen_term_raw = std::pow((last_hinv - first_hinv) -
         (t_max_curve - t_min_curve), 2.0);
     }
     pen_term = lambda_inc * pen_term_raw;
   }
   return loss + static_cast<double>(n) * pen_term;
 }
