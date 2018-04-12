#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_groups );
  DATA_IVECTOR( g_i );
  DATA_VECTOR( y_i );
  
  // Parameters
  PARAMETER( beta0 );
  PARAMETER( log_SD0 );
  PARAMETER( log_SDZ );
  PARAMETER_VECTOR( z_g );
  
  // Objective funcction
  Type jnll = 0;
  int n_i = y_i.size();

  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<n_i; i++){
    jnll -= dnorm( y_i(i), beta0 + z_g(g_i(i)), exp(log_SD0), true );
  }
  
  // Probability of random coefficients
  for( int g=0; g<n_groups; g++){
    jnll -= dnorm( z_g(g), Type(0.0), exp(log_SDZ), true );
  }
  
  // Reporting
  Type SDZ = exp(log_SDZ);
  Type SD0 = exp(log_SD0);

  REPORT( SDZ );
  REPORT( SD0 );
  ADREPORT( SDZ );
  ADREPORT( SD0 );

  return jnll;
}
