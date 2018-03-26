#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_data );
  DATA_INTEGER( n_factors );
  DATA_IVECTOR( Factor );
  DATA_VECTOR( Y );                
  
  // Parameters
  PARAMETER( X0 );
  PARAMETER( log_SD0 );
  PARAMETER( log_SDZ );
  PARAMETER_VECTOR( Z );
  
  // Objective funcction
  Type jnll = 0;
  
  // Probability of data conditional on fixed and random effect values
  for( int i=0; i<n_data; i++){
    jnll -= dnorm( Y(i), X0 + Z(Factor(i)), exp(log_SD0), true );
  }
  
  // Probability of random coefficients
  for( int i=0; i<n_factors; i++){
    jnll -= dnorm( Z(i), Type(0.0), exp(log_SDZ), true );
  }
  
  // Reporting
  Type SDZ = exp(log_SDZ);
  Type SD0 = exp(log_SD0);

  REPORT( SDZ );
  REPORT( SD0 );
  REPORT( Z );
  REPORT( X0 );

  ADREPORT( SDZ );
  ADREPORT( SD0 );
  ADREPORT( Z );
  ADREPORT( X0 );

  return jnll;
}
