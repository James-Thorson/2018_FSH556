#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( nt );
  DATA_VECTOR( y_t );

  // Parameters
  PARAMETER( x0 );
  PARAMETER( log_sigmaP );
  PARAMETER( log_sigmaM );
  PARAMETER( alpha );
  PARAMETER_VECTOR( x_t );
  
  // Objective funcction
  Type jnll = 0;
  
  // Probability of random coefficients
  jnll -= dnorm( x_t(0), x0, exp(log_sigmaP), true );
  for( int t=1; t<nt; t++){
    jnll -= dnorm( x_t(t), x_t(t-1) + alpha, exp(log_sigmaP), true );
  }

  // Probability of data conditional on fixed and random effect values
  for( int t=0; t<nt; t++){
    jnll -= dnorm( y_t(t), x_t(t), exp(log_sigmaM), true );
  }
  
  // Reporting
  Type sigmaP = exp(log_sigmaP);
  Type sigmaM = exp(log_sigmaM);

  REPORT( sigmaP );
  REPORT( sigmaM );
  REPORT( x_t );

  ADREPORT( sigmaP );
  ADREPORT( sigmaM );
  ADREPORT( x_t );

  return jnll;
}
