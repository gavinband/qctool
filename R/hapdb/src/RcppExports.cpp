// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_uncompress_bitpack_haplotypes
IntegerMatrix rcpp_uncompress_bitpack_haplotypes(List rawData, int N);
RcppExport SEXP hapdb_rcpp_uncompress_bitpack_haplotypes(SEXP rawDataSEXP, SEXP NSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type rawData(rawDataSEXP );
        Rcpp::traits::input_parameter< int >::type N(NSEXP );
        IntegerMatrix __result = rcpp_uncompress_bitpack_haplotypes(rawData, N);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

NumericMatrix rcpp_uncompress_floatarray_genotypes(List rawData, int N);
RcppExport SEXP rcpp_uncompress_floatarray_genotypes(SEXP rawDataSEXP, SEXP NSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type rawData(rawDataSEXP );
        Rcpp::traits::input_parameter< int >::type N(NSEXP );
        NumericMatrix __result = rcpp_uncompress_floatarray_genotypes(rawData, N);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
