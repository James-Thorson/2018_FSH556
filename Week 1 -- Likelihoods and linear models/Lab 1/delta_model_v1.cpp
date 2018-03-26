
#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_i );
  DATA_MATRIX( X_ij );

  // Parameters
  PARAMETER_VECTOR( b_j );
  PARAMETER_VECTOR( theta_z );

  // Objective funcction
  Type zero_prob = 1 / (1 + exp(-theta_z(0)));
  Type logsd = exp(theta_z(1));
  Type jnll = 0;
  int n_data = y_i.size();

  // Linear predictor
  vector<Type> linpred_i( n_data );
  linpred_i = X_ij * b_j;

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n_data; i++){
    if(y_i(i)==0) jnll -= log( zero_prob );
    if(y_i(i)!=0) jnll -= log( 1-zero_prob ) + dlognorm( y_i(i), linpred_i(i), logsd, true );
  }
  
  // Reporting
  REPORT( zero_prob );
  REPORT( logsd );
  REPORT( linpred_i );
  return jnll;
}
