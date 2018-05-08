
// State-space Gompertz model
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() () 
{
// data:
DATA_VECTOR(log_yt);
DATA_INTEGER(model_switch);

// parameters:
PARAMETER_VECTOR(pop_par); // population growth rate parameter //density dependence parameter
PARAMETER(log_sigma); // log(process SD) = log(obs SD)
PARAMETER(log_x0);
PARAMETER_VECTOR(log_xt); // state

// transform the parameters
Type sigma = exp(log_sigma);

// reports
ADREPORT(sigma)

int n = log_yt.size(); // get time series length
Type nll = 0.0; // initialize negative log likelihood

nll -= dnorm(log_xt[0], log_x0, sigma, true) ;

Type m;

// process model:
for(int i = 1; i < n; i++){
	if(model_switch == 1){
		m = pop_par[0]+pop_par[1]*log_xt[i-1] ;    //Gompertz model
	}
    if(model_switch == 2){
    	m = log_xt[i-1] + pop_par[0]*(1-exp(log_xt[i-1])/pop_par[1]) ;    //Ricker model
    }
    nll -= dnorm(log_xt[i], m, sigma, true); //likelihood for random effects
}

// observation model:
for(int i = 0; i < n; i++){
  if(!isNA(log_yt[i])){	
  nll -= dnorm(log_yt[i], log_xt[i], sigma, true); //likelihood for observations
  }
 }

SIMULATE{
	log_xt[0] = rnorm(log_x0,sigma);
	log_yt[0] = rnorm(log_xt[0],sigma);
	for (int i = 1; i < n; i++){
		log_xt[i] = rnorm(pop_par[0]+pop_par[1]*log_xt[i-1], sigma);
		log_yt[i] = rnorm(log_xt[i], sigma);
	}
	REPORT(log_xt);
	REPORT(log_yt);
}

return nll;
}


















