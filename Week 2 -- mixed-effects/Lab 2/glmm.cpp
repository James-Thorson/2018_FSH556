#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_y );
  DATA_INTEGER( n_s );
  DATA_IVECTOR( s_i );
  DATA_VECTOR( y_i );
  
  // Parameters
  PARAMETER( beta0 );
  PARAMETER( log_sdz );
  PARAMETER_VECTOR( z_s );
  
  // Objective funcction
  Type jnll = 0;
  
  // Probability of data conditional on fixed and random effect values
  vector<Type> ypred_i(n_y);
  for( int i=0; i<n_y; i++){
    ypred_i(i) = exp( beta0 + z_s(s_i(i)) );
    jnll -= dpois( y_i(i), ypred_i(i), true );
  }
  
  // Probability of random coefficients
  for( int s=0; s<n_s; s++){
    jnll -= dnorm( z_s(s), Type(0.0), exp(log_sdz), true );
  }
  
  // Reporting
  Type sdz = exp(log_sdz);

  REPORT( sdz );

  ADREPORT( sdz );

  return jnll;
}
