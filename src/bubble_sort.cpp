#include <Rcpp.h>
using namespace Rcpp;

//----------------------------------------------------------------------------80
//
//----------------------------------------------------------------------------80
//' Perform one-shot adjacent swapping for each element.
//'
//' @param listdata A list of vector and integer.
//'   
//' @return A List.
//' @export
//'
//' @examples
//' swap_pass(list(vec = c(1, 1, 0), cnt = 0))
//'
// [[Rcpp::export]]
List swap_pass(List listdata){
  int           i, tmp;
  IntegerVector vec = listdata[0];
  int           cnt = listdata[1];
  int           m = vec.length() - 1, n = vec.length();
  IntegerVector newvec(n);
  int           newcnt = cnt;

  for(i=0; i<n; ++i){
    newvec[i] = vec[i];
  }

  if(newcnt % 1000 == 0){
    Rcpp::checkUserInterrupt();
  }

  for(i=0; i<m; ++i){
    if(newvec[i] > newvec[i+1]){
      tmp = newvec[i];
      newvec[i] = newvec[i+1];
      newvec[i+1] = tmp;
      newcnt++;
    }
  }

  List newdata = List::create(newvec, newcnt);
  return newdata;
}
//----------------------------------------------------------------------------80
//
//----------------------------------------------------------------------------80
//' Perform bubble sorting, counting the number of steps.
//'
//' @param listdata A list of vector and integer.
//'   For example, in R code, listdata = list(vec = c(1, 0, 1, ...), cnt = 0).
//'   The integer (cnt = 0) is the initial number of steps for bubble sorting.
//'
//' @return A List.
//' @export
//'
//' @examples
//' bubble_sort(list(vec = c(1, 1, 0), cnt = 0))
//'
// [[Rcpp::export]]
List bubble_sort(List listdata){
  List          newdata = swap_pass(listdata);
  IntegerVector vec = listdata[0], newvec = newdata[0];
  int           i = 0, n = vec.length(), diff = 0;

  for(i=0; i<n; ++i){
    diff += abs(newvec[i] - vec[i]);
  }

  if(diff == 0){
    return newdata;
  }else{
    return bubble_sort(newdata);
  }
}
