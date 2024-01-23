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