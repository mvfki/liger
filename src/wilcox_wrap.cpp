#include <RcppArmadillo.h>

extern "C" {
#include "wilcox.h"
}

//[[Rcpp::export()]]
double wrap_dwilcox(double x, double m, double n, int give_log){
  return Rf_dwilcox(x, m, n, give_log);
};

//[[Rcpp::export()]]
double wrap_pwilcox(double q, double m, double n, int lower_tail, int log_p){
  return Rf_pwilcox(q, m, n, lower_tail, log_p);
};

//[[Rcpp::export()]]
double wrap_qwilcox(double x, double m, double n, int lower_tail, int log_p){
  return Rf_qwilcox(x, m, n, lower_tail, log_p);
};
//[[Rcpp::export()]]
double wrap_rwilcox(double m, double n){
  return Rf_rwilcox(m, n);
}

// [[Rcpp::export]]
arma::mat colAggregateSum_sparse(arma::sp_mat& X,
                                 const arma::uvec& groups,
                                 unsigned ngroups) {
  arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
  for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
    res(groups[it.col()], it.row()) += *it;
  }
  return res;
}

// Non-zero counting aggregate %%%%
// NNZ - Number of Non-Zero

// [[Rcpp::export]]
arma::mat colNNZAggr_sparse(arma::sp_mat& X,
                            const arma::uvec& groups,
                            unsigned ngroups) {
  arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
  for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
    if (*it > 0) res(groups[it.col()], it.row())++;
  }
  return res;
}