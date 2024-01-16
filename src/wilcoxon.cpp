#include <RcppArmadillo.h>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Rcpp;
using namespace arma;

// Algorithm adopted from https://github.com/immunogenomics/presto
// With corrections and removed supports for dense matrices which do not seem
// necessary for LIGER

// X - feature x cell
// output - cell x featureRank
vector<double> cpp_in_place_rank_mean(arma::vec& v_temp, const int idx_begin,
                                      const int idx_end) {
  vector<double> ties;
  
  if (idx_begin > idx_end) return ties;
  vector<pair<double, size_t> > v_sort(idx_end - idx_begin + 1);
  for (int i = idx_begin; i <= idx_end; i++) {
    v_sort[static_cast<size_t>(i) - idx_begin] = make_pair(v_temp[i], i - idx_begin);
  }
  
  
  ranges::sort(v_sort);
  
  double rank_sum = 0;
  double n = 1;
  size_t i;
  for (i = 1U; i < v_sort.size(); i++) {
    if (v_sort[i].first != v_sort[i - 1].first) { // NOLINT(clang-diagnostic-float-equal)
      // if current val != prev val
      // set prev val to something
      for (unsigned j = 0; static_cast<double>(j) < n; j++) {
        v_temp[static_cast<uword>(v_sort[i - 1 - j].second) + idx_begin] =
          rank_sum / n + 1;
      }
      // restart count ranks
      rank_sum = static_cast<double>(i);
      if (n > 1) ties.push_back(n);
      n = 1;
    }
    else {
      // if curr val is a tie, 
      // don't set anything yet, start computing mean rank
      rank_sum += static_cast<float>(i);
      n++;
    }
  }
  // set the last element(s)
  for (unsigned j = 0; static_cast<float>(j) < n; j++)
    v_temp[static_cast<uword>(v_sort[i - 1 - j].second) + idx_begin] = (rank_sum / n) + 1;
  
  return ties;
}


dmat cpp_rank_matrix_dgc(
    SpMat& X) {
  dvec x = nonzeros(X);
  uvec p = vector(X.col_ptrs, X.col_ptrs + X.n_cols + 1);
  dmat ties(X.n_rows, X.n_cols, fill::zeros);
  uvec append_nzero(X.n_rows);
  for (uword i = 0; i < X.n_cols; i++)
  {
    if (p[i + 1] == p[i]) continue;
    const uword n_zero = X.n_rows - p[i + 1] - p[i];
    ties.col(i) = dvec(cpp_in_place_rank_mean(x, static_cast<int>(p[i]),
                       static_cast<int>(p[i + 1]) - 1));
    append_nzero(i) = n_zero;
    X(p(i + 1) - 1) += n_zero;
    
  }
  return ties;
}


arma::mat rowAggregateSum_sparse(arma::sp_mat& X,
                                 const arma::uvec& groups,
                                 unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        res(groups[it.row()], it.col()) += *it;
    }
    return res;
}


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


arma::mat colNNZAggr_sparse(arma::sp_mat& X,
                            const arma::uvec& groups,
                            unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        if (*it > 0) res(groups[it.col()], it.row())++;
    }
    return res;
}

arma::mat rowNNZAggr_sparse(arma::sp_mat& X,
                            const arma::uvec& groups,
                            unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        if (*it > 0) res(groups[it.row()], it.col())++;
    }
    return res;
}

// [[Rcpp::export]]
Rcpp::List cpp_wilcoxauc(const arma::sp_mat &X, std::vector<std::string> y)
{
    if (!ranges::is_sorted(y))
        auto enditer = ranges::sort(y);
    vector<string> uniques;
    ranges::unique_copy(y, uniques.begin());
    uword sized = static_cast<uword>(size(uniques));
    uvec counts(sized);
    uvec countsplus(sized);
    uvec countsminus(sized);
    auto bound_count = bind(count<vector<string>, string>, uniques.begin(), uniques.end(), std::placeholders::_1);
    transform(uniques.begin(), uniques.end(), counts.begin(), bound_count);
    ranges::transform(counts, countsplus.begin(), [](const uword x)
                      { return x + 1; });
    ranges::transform(counts, countsminus.begin(), [](const uword x)
                      { return x - 1; });
    uvec n1n2 = counts * (X.n_cols - counts);
    sp_mat Xt = X.t();
    dmat ties = cpp_rank_matrix_dgc(Xt);
    umat gnz = counts - rowNNZAggr_sparse(Xt,
                                             countsminus, sized);
    uvec r = vector(Xt.col_ptrs, Xt.col_ptrs + Xt.n_cols + 1);
    dvec zeroRanks = conv_to<dvec>::from((Xt.n_rows - diff(r) + 1)) / 2;
    dmat ustat = (gnz.t() * zeroRanks).t() - counts * (counts + 1) / 2;
    dmat auc = (ustat / n1n2).t();
    dmat z = ustat - .5 * conv_to<dvec>::from(n1n2);
    z = z - sign(z) * .5;
    uword x1 = Xt.n_cols ^ 3 - Xt.n_cols;
    double x2 = 1 / (static_cast<double>(12) * (Xt.n_cols ^ 2 - Xt.n_cols));
    fmat rhs(ties.n_cols, ties.n_rows);
    rhs.each_col([x1, x2](const fvec &tvals)
                 { return x1 - (sum(((tvals % tvals) % tvals) - tvals)) * x2; });
    Row flattened = vectorise(rhs);
    fvec usigma = sqrt(n1n2.as_col() * flattened);
    z = (z / usigma).t();
    dvec pvals = 2 * vecnorm(-arma::abs(z));
    Rcpp::Environment stats("package:stats");
    Rcpp::Function pAdjust = stats["p.adjust"];
    dvec fdr = pvals.for_each([pAdjust](const double x)
                              { return pAdjust(x, "BH"); }); // Eventually remove R calls
                                                             // AuxStats
    mat group_sums = colAggregateSum_sparse(X, countsminus, sized);
    umat group_nnz = colNNZAggr_sparse(X, countsminus, sized);
    int cs = 0; // TODO IMPLEMENT LOGFC
    int gs = 0;
    int lfc = 0;
    return Rcpp::List::create(Rcpp::Named("auc") = auc, Rcpp::Named("pval") = pvals, Rcpp::Named("padj") = fdr,
                              Rcpp::Named("pct_in") = group_nnz * 100,
                              Rcpp::Named("pct_out") = 0, // todo fixme
                              Rcpp::Named("avgExpr") = group_sums,
                              Rcpp::Named("statistic") = ustat.t(),
                              Rcpp::Named("logFC") = lfc);
}
