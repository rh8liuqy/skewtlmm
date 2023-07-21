// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fastLikelihood
double fastLikelihood(const NumericVector& covar_params, const double& nu, const arma::mat& Ga, const List& Z, const double& delta, const List& cent, const double& s2, const List& T_i, const List& dd, const NumericVector& ni, const double& N);
RcppExport SEXP _skewtlmm_fastLikelihood(SEXP covar_paramsSEXP, SEXP nuSEXP, SEXP GaSEXP, SEXP ZSEXP, SEXP deltaSEXP, SEXP centSEXP, SEXP s2SEXP, SEXP T_iSEXP, SEXP ddSEXP, SEXP niSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type covar_params(covar_paramsSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Ga(GaSEXP);
    Rcpp::traits::input_parameter< const List& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const List& >::type cent(centSEXP);
    Rcpp::traits::input_parameter< const double& >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< const List& >::type T_i(T_iSEXP);
    Rcpp::traits::input_parameter< const List& >::type dd(ddSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type ni(niSEXP);
    Rcpp::traits::input_parameter< const double& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(fastLikelihood(covar_params, nu, Ga, Z, delta, cent, s2, T_i, dd, ni, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_skewtlmm_fastLikelihood", (DL_FUNC) &_skewtlmm_fastLikelihood, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_skewtlmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
