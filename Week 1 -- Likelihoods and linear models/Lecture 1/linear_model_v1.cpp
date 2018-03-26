// Space time
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_i );

  // Parameters
  PARAMETER( mean );
  PARAMETER( log_sd );

  // Objective funcction
  Type sd = exp(log_sd);
  Type jnll = 0;
  int n_data = y_i.size();

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n_data; i++){
    jnll -= dnorm( y_i(i), mean, sd, true );
  }
  
  // Reporting
  return jnll;
}
