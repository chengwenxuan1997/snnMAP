#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// generic function for kl_divergence
template <typename InputIterator1, typename InputIterator2>
inline double distance(InputIterator1 begin1, 
                       InputIterator1 end1, 
                       InputIterator2 begin2) {
  
  // value to return
  double rval = 0;
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    
    // accumulate if appropirate
    rval += (d1-d2)*(d1-d2);
  }
  return rval;
}


// [[Rcpp::export]]
NumericVector rcpp_distance(NumericMatrix mat, 
                            NumericVector X, 
                            NumericVector Y) {
  
  // allocate the matrix we will return
  NumericVector D(X.size());
  
  for (int i = 0; i < D.size(); i++) {

      // rows we will operate on
      NumericMatrix::Row row1 = mat.row(X[i]-1);
      NumericMatrix::Row row2 = mat.row(Y[i]-1);
      
      // write to output matrix
      D[i] = sqrt(distance(row1.begin(), row1.end(), row2.begin()));
      // D = row2[2];
    }
  
  return D;
}
  
