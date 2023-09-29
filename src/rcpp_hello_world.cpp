#include <Rcpp.h>
using namespace Rcpp;

//' rcpp_hello_world
//' 
//' This function returns a list with character and numeric elements,  demonstrating the use of lists in C++, and testing that Rcpp is integrated properly into the package build process.  
//' @param no no params
//' @export
// [[Rcpp::export]]
List rcpp_hello_world() {
	CharacterVector x =
		CharacterVector::create("foo", "bar");
	NumericVector y =
		NumericVector::create(0.0, 1.0);
	List z = List::create(x, y);
	return z;
}