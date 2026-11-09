// fosr_em.cpp - VL-FoSR EM kernels; optional OpenMP.
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_DONT_USE_OPENMP

#include <RcppArmadillo.h>
using namespace Rcpp;

inline void clamp_min_inplace(arma::vec& x, const double eps = 1e-3) {
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    if (x[j] < eps) x[j] = eps;
  }
}
inline arma::vec make_two_lambda(const arma::vec& zeta) {
  arma::vec sigma = 1.0 / (1.0 + arma::exp(-zeta));
  return (sigma - 0.5) / zeta;
}

//' Compute subject-level response crossproducts.
//'
//' Calculates the time-basis crossproducts with centered binary responses.
//'
//' @param B_full Observation-level B-spline design matrix.
//' @param Y Binary observations.
//' @param first0 Zero-based first row index for each subject.
//' @param last0 Zero-based last row index for each subject.
//' @param cores Number of threads when OpenMP is available.
//' @return A matrix of subject-specific crossproducts.
//' @author Linlin Wu (\email{linlin.wu@emory.edu}),
//'   Jeff Goldsmith (\email{ajg2202@cumc.columbia.edu}),
//'   Julia Wrobel (\email{julia.wrobel@emory.edu})
// [[Rcpp::export]]
 arma::mat subject_BtYhalf_cpp(
     const arma::mat& B_full,
     const arma::vec& Y,
     const arma::uvec& first0,
     const arma::uvec& last0,
     const int cores = 1
 ) {
   if (cores < 1) Rcpp::stop("cores must be >= 1.");
   const arma::uword I = first0.n_elem;
   const arma::uword M = B_full.n_cols;
   if (last0.n_elem != I) {
     Rcpp::stop("first0 and last0 must have the same length.");
   }
   if (B_full.n_rows != Y.n_elem) {
     Rcpp::stop("nrow(B_full) must equal length(Y).");
   }
   arma::mat out(I, M, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(cores)
#endif
   for (long long ii = 0; ii < static_cast<long long>(I); ++ii) {
     const arma::uword i = static_cast<arma::uword>(ii);
     const arma::uword a = first0[i];
     const arma::uword b = last0[i];
     arma::vec yc = Y.subvec(a, b) - 0.5;
     out.row(i) = (B_full.rows(a, b).t() * yc).t();
   }
   return out;
 }

//' Compute subject-specific fitted linear predictors.
//'
//' @param B_full Observation-level B-spline design matrix.
//' @param subject_coefs Matrix of subject-specific time-basis coefficients.
//' @param first0 Zero-based first row index for each subject.
//' @param last0 Zero-based last row index for each subject.
//' @param cores Number of threads when OpenMP is available.
//' @return Vector of fitted linear predictors at the observed times.
//' @author Linlin Wu (\email{linlin.wu@emory.edu}),
//' Jeff Goldsmith (\email{ajg2202@cumc.columbia.edu}),
//' Julia Wrobel (\email{julia.wrobel@emory.edu})
// [[Rcpp::export]]
 arma::vec fitted_values_by_subject_cpp(
     const arma::mat& B_full,
     const arma::mat& subject_coefs,
     const arma::uvec& first0,
     const arma::uvec& last0,
     const int cores = 1
 ) {
   if (cores < 1) Rcpp::stop("cores must be >= 1.");
   const arma::uword I = first0.n_elem;

   if (last0.n_elem != I) {
     Rcpp::stop("first0 and last0 must have the same length.");
   }
   if (subject_coefs.n_cols != I) {
     Rcpp::stop("ncol(subject_coefs) must equal the number of subjects.");
   }
   if (subject_coefs.n_rows != B_full.n_cols) {
     Rcpp::stop("nrow(subject_coefs) must equal ncol(B_full).");
   }
   arma::vec fits(B_full.n_rows, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(cores)
#endif
   for (long long ii = 0; ii < static_cast<long long>(I); ++ii) {
     const arma::uword i = static_cast<arma::uword>(ii);
     const arma::uword a = first0[i];
     const arma::uword b = last0[i];
     fits.subvec(a, b) = B_full.rows(a, b) * subject_coefs.col(i);
   }
   return fits;
 }

//' Perform one variational EM update for VL-FoSR.
//'
//' Updates the subject scores, variational parameters, and spline coefficients.
//'
//' @param B_full Observation-level B-spline design matrix.
//' @param Y Binary observations.
//' @param zeta_in Current variational parameters.
//' @param Q_by_id Subject-level fixed-effect design matrix.
//' @param alpha_mat Current fixed-effect spline coefficient matrix.
//' @param psi_coefs Current random-effect spline coefficient matrix.
//' @param BtYhalf Precomputed subject-level response crossproducts.
//' @param first0 Zero-based first row index for each subject.
//' @param last0 Zero-based last row index for each subject.
//' @param cores Number of threads when OpenMP is available.
//' @return A list with updated coefficients (phi_vec), subject scores,
//'   variational parameters (zeta), and accumulated EM matrices (H_all, g_all).
//' @author Linlin Wu (\email{linlin.wu@emory.edu}),
//' Jeff Goldsmith (\email{ajg2202@cumc.columbia.edu}),
//' Julia Wrobel (\email{julia.wrobel@emory.edu})
// [[Rcpp::export]]
 Rcpp::List em_iteration_structured_cpp(
     const arma::mat& B_full,
     const arma::vec& Y,
     const arma::vec& zeta_in,
     const arma::mat& Q_by_id,
     const arma::mat& alpha_mat,
     const arma::mat& psi_coefs,
     const arma::mat& BtYhalf,
     const arma::uvec& first0,
     const arma::uvec& last0,
     const int cores = 1
 ) {
   if (cores < 1) Rcpp::stop("cores must be >= 1.");
   const arma::uword I   = first0.n_elem;
   const arma::uword N   = B_full.n_rows;
   const arma::uword M   = B_full.n_cols;
   const arma::uword Q   = Q_by_id.n_cols;
   const arma::uword npc = psi_coefs.n_cols;

   if (last0.n_elem != I) {
     Rcpp::stop("first0 and last0 must have the same length.");
   }
   if (Y.n_elem != N || zeta_in.n_elem != N) {
     Rcpp::stop("B_full, Y, and zeta_in must have compatible row counts.");
   }
   if (Q_by_id.n_rows != I) {
     Rcpp::stop("nrow(Q_by_id) must equal the number of subjects.");
   }
   if (alpha_mat.n_rows != M || alpha_mat.n_cols != Q) {
     Rcpp::stop("alpha_mat must be M x ncol(Q_by_id).");
   }
   if (psi_coefs.n_rows != M) {
     Rcpp::stop("nrow(psi_coefs) must equal ncol(B_full).");
   }
   if (BtYhalf.n_rows != I || BtYhalf.n_cols != M) {
     Rcpp::stop("BtYhalf must be I x M.");
   }
   if (Q < 1) {
     Rcpp::stop("Q_by_id must contain at least the intercept column.");
   }
   if (npc < 1) {
     Rcpp::stop("npc must be at least 1.");
   }

   // Fixed-effect blocks precede random-effect blocks in the coefficient vector.
   const arma::uword G = Q + npc;
   const arma::uword n_param = M * G;
   arma::mat scores_out(I, npc, arma::fill::zeros);
   arma::vec zeta_out = zeta_in;
   arma::mat H_all(n_param, n_param, arma::fill::zeros);
   arma::vec g_all(n_param, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp parallel num_threads(cores)
#endif
{
  // Thread-local accumulators avoid contention inside the subject loop.
  arma::mat H_local(n_param, n_param, arma::fill::zeros);
  arma::vec g_local(n_param, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (long long ii = 0; ii < static_cast<long long>(I); ++ii) {
    const arma::uword i = static_cast<arma::uword>(ii);
    const arma::uword a = first0[i];
    const arma::uword b = last0[i];
    const auto B_i = B_full.rows(a, b);
    arma::vec Y_i = Y.subvec(a, b);
    arma::vec z_old = zeta_in.subvec(a, b);
    clamp_min_inplace(z_old, 1e-3);
    arma::vec q_i = Q_by_id.row(i).t();
    arma::vec fixed_coef_i = alpha_mat * q_i;
    arma::vec a_i = B_i * fixed_coef_i;
    arma::vec two_lambda_old = make_two_lambda(z_old);
    arma::vec Y_centered = Y_i - 0.5;
    arma::vec mi(npc, arma::fill::zeros);
    arma::mat Ci(npc, npc, arma::fill::zeros);
    arma::vec zeta_new;
    if (npc == 1) {
      arma::vec w = B_i * psi_coefs.col(0);
      const double h =
        1.0 + arma::dot(w, two_lambda_old % w);
      const double ci = 1.0 / ((h > 1e-12) ? h : (h + 1e-8));
      arma::vec g_vec =
        Y_centered - two_lambda_old % a_i;
      const double g = arma::dot(w, g_vec);
      const double m = ci * g;
      mi[0] = m;
      Ci(0, 0) = ci;
      arma::vec mu_star = a_i + w * m;
      arma::vec var_star = ci * arma::square(w);
      zeta_new = arma::sqrt(
        arma::square(mu_star) + var_star
      );
    } else {
      arma::mat W_i = B_i * psi_coefs;
      arma::mat WL = W_i.each_col() % two_lambda_old;
      arma::mat H = arma::eye<arma::mat>(npc, npc) +
        W_i.t() * WL;
      arma::vec g_vec = Y_centered - two_lambda_old % a_i;
      arma::vec g = W_i.t() * g_vec;
      bool ok = arma::inv_sympd(Ci, H);
      if (!ok) {
        arma::mat H_stab =
          H + 1e-8 * arma::eye<arma::mat>(npc, npc);
        ok = arma::inv_sympd(Ci, H_stab);
      }
      if (!ok) {
        Ci = arma::pinv(H);
      }
      mi = Ci * g;
      arma::vec mu_star = a_i + W_i * mi;
      arma::mat tmp = W_i * Ci;
      arma::vec var_star = arma::sum(tmp % W_i, 1);
      zeta_new = arma::sqrt(
        arma::square(mu_star) + var_star
      );
    }

    clamp_min_inplace(zeta_new, 1e-3);
    scores_out.row(i) = mi.t();
    zeta_out.subvec(a, b) = zeta_new;
    arma::vec two_lambda_new = make_two_lambda(zeta_new);
    arma::mat B_weighted = B_i;
    B_weighted.each_col() %= two_lambda_new;
    arma::mat S_i = B_i.t() * B_weighted;
    arma::vec r_i = BtYhalf.row(i).t();
    arma::mat Exx = Ci + mi * mi.t();
    arma::vec c_i = arma::join_cols(q_i, mi);
    arma::mat T_i(G, G, arma::fill::zeros);
    T_i.submat(0, 0, Q - 1, Q - 1) = q_i * q_i.t();
    T_i.submat(0, Q, Q - 1, G - 1) = q_i * mi.t();
    T_i.submat(Q, 0, G - 1, Q - 1) = mi * q_i.t();
    T_i.submat(Q, Q, G - 1, G - 1) = Exx;
    for (arma::uword u = 0; u < G; ++u) {
      const arma::uword us = u * M;
      g_local.subvec(us, us + M - 1) += c_i[u] * r_i;
      for (arma::uword v = 0; v < G; ++v) {
        const arma::uword vs = v * M;
        H_local.submat(
          us, vs,
          us + M - 1,
          vs + M - 1
        ) += T_i(u, v) * S_i;
      }
    }
  }

#ifdef _OPENMP
#pragma omp critical
#endif
{
  H_all += H_local;
  g_all += g_local;
}
}

arma::vec phi_vec;

try {
  phi_vec = arma::solve(
    H_all,
    g_all,
    arma::solve_opts::likely_sympd
  );
} catch (...) {
  phi_vec = arma::solve(
    H_all +
      1e-8 * arma::eye<arma::mat>(H_all.n_rows, H_all.n_cols),
      g_all
  );
}

return Rcpp::List::create(
  Rcpp::Named("phi_vec") = phi_vec,
  Rcpp::Named("scores")  = scores_out,
  Rcpp::Named("zeta")    = zeta_out,
  Rcpp::Named("H_all")   = H_all,
  Rcpp::Named("g_all")   = g_all
);
 }
