// Space time
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_i );
  DATA_MATRIX( X_ij )

  // Parameters
  PARAMETER_VECTOR( b_j );
  PARAMETER( log_sd );

  // Objective funcction
  Type sd = exp(log_sd);
  Type jnll = 0;
  int n_data = y_i.size();

  // Linear predictor
  vector<Type> linpred_i( n_data );
  linpred_i = X_ij*b_j;

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n_data; i++){
    jnll -= dnorm( y_i(i), linpred_i(i), sd, true );
  }
  
  // Reporting
  return jnll;
}
