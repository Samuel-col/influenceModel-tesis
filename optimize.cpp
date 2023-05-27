// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <math.h>
#include <iomanip>
// [[Rcpp::depends(RcppArmadillo)]]

double cpp_f(Rcpp::Function f,double x,double target_value){
  return Rcpp::as<double>(f(x)) - target_value;
}

// [[Rcpp::export]]
Rcpp::List cpp_optimize(Rcpp::Function f,double target_value){
  int maxit = 10000;
  double tol = 1e-7;
  double p0 = -7, p1 = -11, p;
  
  while(std::fabs(p1 - p0) > tol && maxit > 0 &&
        fabs(cpp_f(f,p1,target_value) - cpp_f(f,p0,target_value)) > tol){
    p  = p1;
    p1 = p0 - cpp_f(f,p0,target_value)/(cpp_f(f,p1,target_value) - cpp_f(f,p0,target_value))*(p1 - p0);
    p0 = p;
    maxit--;
  }
  double meanDeg = cpp_f(f,p1,0);
  std::cout << "shift : " << std::setprecision(6) <<
    p1 << " | MeanDeg : " << std::setprecision(3) <<
    meanDeg << std::endl; 
  
  return Rcpp::List::create(Rcpp::Named("shift",p1),
                            Rcpp::Named("meanDeg",meanDeg));
}