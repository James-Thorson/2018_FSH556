#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_i );
  DATA_INTEGER( n_s );
  DATA_IVECTOR( s_i );
  DATA_VECTOR( y_i );
  
  // Parameters
  PARAMETER( beta0 );
  PARAMETER( log_sd_epsilon );
  PARAMETER( log_sd_delta );
  PARAMETER_VECTOR( epsilon_s );
  PARAMETER_VECTOR( delta_i );

  // Objective funcction
  Type jnll = 0;
  
  // Probability of data conditional on fixed and random effect values
  vector<Type> ypred_i(n_i);
  for( int i=0; i<n_i; i++){
    ypred_i(i) = exp( beta0 + epsilon_s(s_i(i)) + delta_i(i) );
    jnll -= dpois( y_i(i), ypred_i(i), true );
  }
  
  // Probability of site-level variation
  for( int s=0; s<n_s; s++){
    jnll -= dnorm( epsilon_s(s), Type(0.0), exp(log_sd_epsilon), true );
  }
  
  // Probability of overdispersion
  for( int i=0; i<n_i; i++){
    jnll -= dnorm( delta_i(i), Type(0.0), exp(log_sd_delta), true );
  }

  // Reporting
  Type sd_epsilon = exp(log_sd_epsilon);
  Type sd_delta = exp(log_sd_delta);

  REPORT( sd_epsilon );
  REPORT( sd_delta );
  REPORT( epsilon_s );
  REPORT( delta_i );
  REPORT( beta0 );
  
  ADREPORT( sd_epsilon );
  ADREPORT( sd_delta );

  return jnll;
}
