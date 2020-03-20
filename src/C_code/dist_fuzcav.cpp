#include <Rcpp.h>

using namespace Rcpp;

double min (double a , double b){
	if (a<b) return a;
	return b;
}

// [[Rcpp::export]]
double dist_fuzcav_ph(NumericVector a , NumericVector b){
	double res = 0;
	double countA = 0;
	double countB = 0;
	for (int i = 0 ; i< 4834; i++){
    if (b[i] !=0 && a[i] != 0) res++;
		if (a[i]!=0) countA++;
		if (b[i]!=0) countB++;
		
	}
	double mymin = min(countA,countB);
	if (mymin == 0) return 0;
	res = res /mymin ;
	return res;
}

