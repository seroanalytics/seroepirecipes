#include <Rcpp.h>
using namespace Rcpp;

//' Fast implementation of the power-law waning antibody kinetics function
//' @export
//' @inheritParams kinetics_power_function
// [[Rcpp::export]]
NumericVector kinetics_power_function_cpp(NumericVector times, NumericVector pars) {
  int n_times = times.size();
  NumericVector titers(n_times);
  double y0 = pars["y0"];
  double peak = pars["peak"];
  double k = pars["k"];
  double r = pars["r"];
  double v = pars["v"];
  double t_peak = pars["t_peak"];
  
  peak = (peak - 1.0) * exp(-y0*k) + 1.0;
  double mu = (1.0/t_peak) * log(peak);
    
  double peak_reached = peak*y0;
  
  for(int t = 0; t < n_times; t++){
    if(times[t] <= t_peak){
      titers[t] =  y0*exp(mu*times[t]);
    } else {
      titers[t] = peak_reached * pow((1.0 + v*(r - 1.0)*(pow(peak_reached,(r-1.0)))*(times[t]-t_peak)),(-1/(r-1)));
    }
  }
  return titers;
}

//' Fast observation error function continuous
//' Calculate the probability of a set of observed titres given a corresponding set of predicted titres assuming continuous, bounded observations. FAST IMPLEMENTATION
//' @name Fast observation error function continuous
//' @param theta NumericVector, a named parameter vector giving the normal distribution standard deviation and the max observable titre
//' @param obs NumericVector, the vector of observed log titres
//' @param predicted_titres NumericVector, the vector of predicted log titres
//' @param a vector of same length as the input data giving the probability of observing each observation given the predictions
//' @return a likelihood for each observed titre
//' @export
// [[Rcpp::export(rng = false)]]
NumericVector likelihood_func_fast_continuous(const NumericVector &theta, const NumericVector &obs, const NumericVector &predicted_titres){
 int total_titres = predicted_titres.size();
 NumericVector ret(total_titres);
 const double sd = theta["error"];
 const double den = sd*M_SQRT2;
 const double den2 = log(sd*2.50662827463);
 const double max_titre = theta["MAX_TITRE"];
 const double min_titre = theta["MIN_TITRE"];
 const double log_const = log(0.5);
 
 for(int i = 0; i < total_titres; ++i){
   // Most titres are between 0 and max_titre, this is the difference in normal cdfs
   if(obs[i] < max_titre && obs[i] > min_titre){
     ret[i] = -0.5*(pow((obs[i]-predicted_titres[i])/sd, 2)) - den2;
     // For titres above the maximum, 
   } else if(obs[i] >= max_titre) {
     ret[i] = log_const + log(erfc((max_titre - predicted_titres[i])/den));
   } else {
     ret[i] = log_const + log(1.0 + erf((theta["MIN_TITRE"] - predicted_titres[i])/den));
   }
 }
 return(ret);
} 

