#include <TMB.hpp>

template <class Type>
Type logitSelectivity(int a, Type mu, Type ups){
  Type tmp = Type(1) / (Type(1) + exp(Type(-1.0) * ups * (a - mu)));
  return tmp;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data inputs ---------------------------
  DATA_VECTOR(ages);
  DATA_VECTOR(paaC);      // prop at age female - commercial
  DATA_SCALAR(mu_M);      // M prior female- not log space
  DATA_SCALAR(sd_M);      // M prior female sd
  DATA_SCALAR(mu_F);      // F prior female - not log space
  DATA_SCALAR(sd_F);      // F prior female sd
  // DATA_INTEGER(slx_type); //selectivity type - 0 = logistic, 1 - dome

  
  // Parameters -----------------------------
  PARAMETER(logM);		    // log nat mort rate female
  PARAMETER(logF); 		  // log female F commercial current
  PARAMETER(logmu);		    // log age full selectivity females, comm
  PARAMETER(logupsilon);	// log scale selectivity curve, comm
  PARAMETER(logsigR);
  // PARAMETER(logmax_sel);
  // PARAMETER(logmin_sel);
  // PARAMETER(logmu2);
  // PARAMETER(logupsilon2);
  
  // Setup, procedures, init -------------
  
  int A = ages.size(); 	// ages a=1,...,A=30
  Type M = exp(logM);
  Type F = exp(logF);
  Type mu = exp(logmu);
  Type upsilon = exp(logupsilon);
  Type sigR = exp(logsigR);
  // Type max_sel = exp(logmax_sel);
  // Type min_sel = exp(logmin_sel);
  // Type mu2 = exp(logmu2);
  // Type upsilon2 = exp(logupsilon2);
  
  Type expnegM = exp(-M);
  vector<Type> Na(A); // Naa 
  vector<Type> Va(A); // unfished 
  vector<Type> Fa(A);

  vector<Type> saC(A); // selectivity females - comm
  vector<Type> Ca(A); // catch females
  vector<Type> propC(A); // prop females - comm

  Type nll = 0.0; //init neg loglik

  // Priors --------------------------
  nll -=dnorm(M, mu_M, sd_M, true);
  nll -=dnorm(F, mu_F, sd_F, true);

  // State dynamics ----------------------
  
  
  for(int a=0; a<A; a++){
    // switch (slx_type) {
    // case 0: // logistic selectivity
        saC(a) = logitSelectivity(a+1, mu, upsilon);
      // break;
    // case 1: // double normal
      //   
      // if(a < mu) saC(a) = pow(2.0, -((a - mu) / upsilon * (a - mu) / upsilon));
      // if(a >= mu) saC(a) = pow(2.0, -((a - mu) / upsilon2 * (a - mu) / upsilon2));
      // 
      // break;    
    // }
  }
  
  // Type maxsa = max(saC);
  // for(int a=0; a<A; a++){
  //   saC(a) = saC(a) / maxsa;
  // }
  
  for(int a=0; a<A; a++){
    Fa(a) = saC(a) * F;
  }
  
  Na(0) = 1000 ;
  Va(0) = 1000; 
  
  for(int a=1; a<(A-1); a++){
    Na(a) = Na(a-1) * exp(Type(-1.0) * (M + Fa(a-1))); 
    Va(a) = Va(a-1) * expnegM;
  }
  
  Na(A-1) = (Na(A-2) * exp(Type(-1.0) * (M + Fa(A-2)))) / (Type(1.0) - exp(Type(-1.0) * (M + Fa(A-1))));
  Va(A-1) = (Va(A-2) * expnegM) / (Type(1.0) - expnegM);
  
  for(int a=0; a<A; a++){
    Ca(a) = Na(a) * (Type(1.0) - exp(Type(-1.0) * M - Fa(a))) * Fa(a) / (M + Fa(a));
  }
  
  Type maxf = max(Ca);

  for(int a=0; a<A; a++){
    propC(a) = Ca(a) / maxf;
  }
  
  for(int a=0; a<A; a++){
    nll -=dnorm(propC(a), paaC(a), sigR, true);
  }
  

  // Reports -------------------------------------
  
  REPORT(M);
  REPORT(F);
  REPORT(Fa);
  REPORT(saC);
  REPORT(propC);
  REPORT(mu);
  REPORT(upsilon);
  // REPORT(mu2);
  // REPORT(upsilon2);
  // REPORT(min_sel);
  // REPORT(max_sel);
  REPORT(Ca);
  REPORT(Na);
  REPORT(Va);
  
  return nll;
  
}
