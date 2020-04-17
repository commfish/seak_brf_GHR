#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // data -------------------------------
  DATA_VECTOR(age);
  DATA_VECTOR(length);
  int n = age.size(); 
  
  // parameters -------------------------
  PARAMETER(logLinf);
  PARAMETER(logkappa);
  PARAMETER(t0);
  PARAMETER(logSigma);
  
  // setup, procedures, inits -----------
  Type Linf = exp(logLinf);
  Type kappa = exp(logkappa);
  Type Sigma = exp(logSigma);
  
  vector<Type> fit(n);
  Type nll = 0.0; // negative log likelihood
  
  // analysis --------------------------
  
  for(int i=0; i<n; i++){
    fit(i) = Linf * (1.0 - exp(-kappa * (age(i) - t0)));
  }
  
  nll = -sum(dnorm(length, fit, Sigma, true));
  
  
  // reports ---------------------------
  REPORT(Linf);
  REPORT(kappa);
  REPORT(t0);
  REPORT(Sigma);
  REPORT(fit);
  
  return nll;
}
