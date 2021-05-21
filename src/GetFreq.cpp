#include <Rcpp.h>
using namespace Rcpp;

//' Compute .. frequency values
//'
//' Function GetFreq description
//'
//' @param Freq Frequency matrix
//' @export
// [[Rcpp::export]]
CharacterVector GetFreq(NumericMatrix Freq ) {

  int Nrow = Freq.nrow();
  std::vector<std::string> resv(Nrow);
  const double epsilon=0.0000001;
  double f_ref;
  int i;

  for( i = 0; i < Nrow; ++i ){


    if ( Freq(i,0) < epsilon && Freq(i,1) < epsilon && Freq(i,2) < epsilon){
      resv[i] = "1";
      continue;
    }

    f_ref =  1-( ( Freq(i,0) + Freq(i,1) + Freq(i,2) ) / 2.);

    resv[i] == "";
    if (f_ref >= epsilon ) {
      resv[i] = std::to_string( f_ref );
      if (Freq(i, 0) >= epsilon ) resv[i] = resv[i] + "/" + std::to_string( Freq(i, 0)/2. );
    } else {
      if (Freq(i, 0) >= epsilon ) resv[i] = std::to_string( Freq(i, 0)/2. );
    }

    if (Freq(i, 1) >= epsilon ) resv[i] = resv[i] + "/" + std::to_string( Freq(i, 1)/2. );
    if (Freq(i, 2) >= epsilon ) resv[i] = resv[i] + "/" + std::to_string( Freq(i, 2)/2. );

    continue;
  }

  return(wrap(resv));


}

