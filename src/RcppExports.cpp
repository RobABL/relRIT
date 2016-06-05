// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_sim
double compute_sim(List const& interaction, DataFrame const& data, LogicalVector const& isFactor);
RcppExport SEXP RIT_compute_sim(SEXP interactionSEXP, SEXP dataSEXP, SEXP isFactorSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List const& >::type interaction(interactionSEXP);
    Rcpp::traits::input_parameter< DataFrame const& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< LogicalVector const& >::type isFactor(isFactorSEXP);
    __result = Rcpp::wrap(compute_sim(interaction, data, isFactor));
    return __result;
END_RCPP
}
// cpp_Relaxed_RIT
Rcpp::List cpp_Relaxed_RIT(Rcpp::List const& datas, Rcpp::NumericVector const& theta, Rcpp::LogicalVector const& factor, double epsilon_cont, double epsilon_cat, int n_trees, int depth, int branch, int min_inter_sz, bool es);
RcppExport SEXP RIT_cpp_Relaxed_RIT(SEXP datasSEXP, SEXP thetaSEXP, SEXP factorSEXP, SEXP epsilon_contSEXP, SEXP epsilon_catSEXP, SEXP n_treesSEXP, SEXP depthSEXP, SEXP branchSEXP, SEXP min_inter_szSEXP, SEXP esSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List const& >::type datas(datasSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector const& >::type factor(factorSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_cont(epsilon_contSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon_cat(epsilon_catSEXP);
    Rcpp::traits::input_parameter< int >::type n_trees(n_treesSEXP);
    Rcpp::traits::input_parameter< int >::type depth(depthSEXP);
    Rcpp::traits::input_parameter< int >::type branch(branchSEXP);
    Rcpp::traits::input_parameter< int >::type min_inter_sz(min_inter_szSEXP);
    Rcpp::traits::input_parameter< bool >::type es(esSEXP);
    __result = Rcpp::wrap(cpp_Relaxed_RIT(datas, theta, factor, epsilon_cont, epsilon_cat, n_trees, depth, branch, min_inter_sz, es));
    return __result;
END_RCPP
}