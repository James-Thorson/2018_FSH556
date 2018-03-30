#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_i );

  // Parameters
  PARAMETER( log_mean );
  PARAMETER( log_CV );

  // Objective funcction
  Type jnll = 0;
  Type mean = exp(log_mean);
  Type CV = exp(log_CV);
  Type shape = pow(CV,-2);
  Type scale = mean * pow(CV,2);

  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<y_i.size(); i++){
    jnll -= dgamma( y_i(i), shape, scale, true );
  }
  
  REPORT( mean );
  REPORT( CV );

  ADREPORT( mean );
  ADREPORT( CV );

  return jnll;
}
