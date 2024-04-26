#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List EstNull_cpp(NumericVector x, double gamma = 0.1) {
  int n = x.length();
  double gan = std::pow(n, -gamma);
  
  double s;
  double phiplus;
  double phiminus;
  double phi;
  double dphiplus;
  double dphiminus;
  double shat = 0;
  double uhat = 0;
  
  for (int t = 1; t <= 1000; t++){
    s = t * 1.0 / 200;
    phiplus = mean(cos(s * x));
    phiminus = mean(sin(s * x));
    phi = std::sqrt(std::pow(phiplus, 2) + std::pow(phiminus, 2));
    if (phi <= gan){
      dphiplus = -mean(x*sin(s*x));
      dphiminus = mean(x*cos(s*x));
      shat = std::sqrt(-(phiplus * dphiplus + phiminus * dphiminus) / (s * phi * phi));
      uhat = -(dphiplus * phiminus - dphiminus * phiplus) / (phi * phi);
      break;
    }
  }
  
  return List::create(Rcpp::Named("mean") = uhat,
                      Rcpp::Named("std") = shat);
}

