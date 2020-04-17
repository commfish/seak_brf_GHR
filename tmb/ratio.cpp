#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(male);
  DATA_VECTOR(female);
  DATA_INTEGER(model);
  int n = male.size();
  
  PARAMETER(p);
  //PARAMETER(dummy);
  PARAMETER(tau);
  
  Type logit_p = log(p / (Type(1.0) - p));
  
  
  vector<Type> pfit(n);
  vector<Type> sigma_y(n);
  
  Type nll = 0.0;
  
  
  for (int i=0;i<n;i++) {
    
    // prediction
    pfit(i) = female(i) / (male(i) + female(i));    
    sigma_y(i) = sqrt(pow(tau, Type(2.0)) + pfit(i) * (Type(1.0) - pfit(i)) / (male(i) + female(i)));
    
    switch(model) {
    case 1 : // Model 1 (binomial)
      nll -= dbinom_robust(male(i), male(i) + female(i), logit_p, true);
      break;
      
    case 2 : // Model 2 (binomial with overdispersion)
      nll += log(sigma_y(i)) + Type(1.0) / (Type(2.0) * pow(sigma_y(i), Type(2.0))) * (pow(pfit(i) - p, Type(2.0)));
      break;
      
      
    }
  }
  
  
  REPORT(sigma_y);
  REPORT(pfit);
  return nll;
}
