// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rollmean
NumericVector rollmean(NumericVector x, int k, bool na_pad, bool na_rm, std::string align);
RcppExport SEXP _cragr_rollmean(SEXP xSEXP, SEXP kSEXP, SEXP na_padSEXP, SEXP na_rmSEXP, SEXP alignSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type na_pad(na_padSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< std::string >::type align(alignSEXP);
    rcpp_result_gen = Rcpp::wrap(rollmean(x, k, na_pad, na_rm, align));
    return rcpp_result_gen;
END_RCPP
}
// rollsum
NumericVector rollsum(NumericVector x, int k, bool na_pad, bool na_rm, std::string align);
RcppExport SEXP _cragr_rollsum(SEXP xSEXP, SEXP kSEXP, SEXP na_padSEXP, SEXP na_rmSEXP, SEXP alignSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type na_pad(na_padSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< std::string >::type align(alignSEXP);
    rcpp_result_gen = Rcpp::wrap(rollsum(x, k, na_pad, na_rm, align));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cragr_rollmean", (DL_FUNC) &_cragr_rollmean, 5},
    {"_cragr_rollsum", (DL_FUNC) &_cragr_rollsum, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_cragr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
