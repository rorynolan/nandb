// [[Rcpp::export]]
int SplitEverykGives(IntegerVector vec, int k) {
  return(sum(ceiling(as<NumericVector>(vec) / k)));
}
