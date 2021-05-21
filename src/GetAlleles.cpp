#include <Rcpp.h>
using namespace Rcpp;

//' Compute .. allels
//'
//' Function GetAlleles description
//'
//' @param Both A character vector vector
//' @param Allele An integer matrix
//' @param n_sample Number of elements
//' @export
// [[Rcpp::export]]
CharacterVector GetAlleles(std::vector<std::string> Both,
                           IntegerMatrix Allele, int n_sample ) {
  int Nrow = Allele.nrow();
  int Nboth = Nrow/n_sample;
  std::vector<std::string> resv(Nrow);

  int i, j;

  for( i = 0; i < Allele.nrow(); ++i ){

    if(Allele(i, 2)==-1){

      if ( Allele(i, 0) == Allele(i, 1) ){
        resv[i] = Both[ Nboth * Allele(i, 0) + (i % Nboth)  ] ;
        continue;
      } else {
        resv[i] = Both[ Nboth * Allele(i, 0) + (i % Nboth)  ] + '/'+
                  Both[ Nboth * Allele(i, 1) + (i % Nboth)  ];
        continue;
      }

    }

    if(Allele(i, 3) == -1){
      resv[i] = Both[ Nboth * Allele(i, 0) + (i % Nboth)  ] + '/'+
                Both[ Nboth * Allele(i, 1) + (i % Nboth)  ] + '/'+
                Both[ Nboth  *Allele(i, 2) + (i % Nboth)  ];
      continue;

    }
    else{
      resv[i] = Both[ Nboth * Allele(i, 0) + (i % Nboth)  ] + '/'+
                Both[ Nboth * Allele(i, 1) + (i % Nboth)  ] + '/'+
                Both[ Nboth * Allele(i, 2) + (i % Nboth)  ] + '/'+
                Both[ Nboth * Allele(i, 3) + (i % Nboth)  ];
      continue;

    }
  }

  return(wrap(resv));


}

