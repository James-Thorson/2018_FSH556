
#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// Main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR( y_i );
  DATA_MATRIX( X_ij )
  DATA_VECTOR( predTF_i );

  // Parameters
  PARAMETER_VECTOR( b_j );
  PARAMETER_VECTOR( theta_z );

  // Objective funcction
  Type zero_prob = 1 / (1 + exp(-theta_z(0)));
  Type logsd = exp(theta_z(1));
  int n_data = y_i.size();
  int n_j = X_ij.row(0).size();

  vector<Type> jnll_i(n_data);
  Type jnll = 0;
  Type pred_jnll = 0;

  // Linear predictor
  vector<Type> linpred_i( n_data );
  for( int i=0; i<n_data; i++){
    linpred_i(i) = 0;
    for( int j=0; j<n_j; j++){
      linpred_i(i) += X_ij(i,j) * b_j(j);
    }
  }

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n_data; i++){
    if(y_i(i)==0) jnll_i(i) -= log( zero_prob );
    if(y_i(i)!=0) jnll_i(i) -= log( 1-zero_prob ) + dlognorm( y_i(i), linpred_i(i), logsd, true );
    // Running counter
    if( predTF_i(i)==0 ) jnll += jnll_i(i);
    if( predTF_i(i)==1 ) pred_jnll += jnll_i(i);
  }

  // Reporting
  REPORT( zero_prob );
  REPORT( logsd );
  REPORT( linpred_i );
  REPORT( pred_jnll );
  REPORT( jnll_i );
  return jnll;
}
