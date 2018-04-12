#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( nt );
  DATA_VECTOR( log_b_t );
  DATA_VECTOR( SD_log_b_t );

  // Parameters
  PARAMETER( log_d0 );
  PARAMETER( log_sigmaP );
  PARAMETER( log_sigmaM );
  PARAMETER( alpha );
  PARAMETER( rho );
  PARAMETER_VECTOR( log_d_t );
  
  // Objective funcction
  Type jnll = 0;
  
  // Reporting
  Type sigmaP = exp(log_sigmaP);
  vector<Type> sigmaM_t(nt);
  for( int t=0; t<nt; t++){
    sigmaM_t(t) = pow( pow(exp(log_sigmaM),2) + pow(SD_log_b_t(t),2), 0.5 );
  }

  // Probability of random coefficients
  jnll -= dnorm( log_d_t(0), log_d0, exp(log_sigmaP), true );
  for( int t=1; t<nt; t++){
    jnll -= dnorm( log_d_t(t), alpha + rho*log_d_t(t-1), sigmaP, true );
  }

  // Probability of data conditional on fixed and random effect values
  for( int t=0; t<nt; t++){
    jnll -= dnorm( log_b_t(t), log_d_t(t), sigmaM_t(t), true );
  }
  
  REPORT( sigmaP );
  REPORT( sigmaM_t );
  REPORT( log_d_t );

  ADREPORT( sigmaP );
  ADREPORT( sigmaM_t );

  return jnll;
}
