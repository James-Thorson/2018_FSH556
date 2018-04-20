
#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// dlnorm
template<class Type>
Type dinvgauss(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.5*log(shape) - 0.5*log(2*M_PI*pow(x,3)) - (shape * pow(x-mean,2) / (2*pow(mean,2)*x));
  if(give_log) return logres; else return exp(logres);
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Options
  DATA_IVECTOR( Options_vec );

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
  vector<Type> jnll_i(n_data);
  jnll_i.setZero();
  Type jnll = 0;
  Type pred_jnll = 0;

  // Linear predictor
  vector<Type> linpred_i( n_data );
  linpred_i = X_ij*b_j;

  // Probability of data conditional on fixed effect values
  for( int i=0; i<n_data; i++){
    if(y_i(i)==0){
      jnll_i(i) -= log( zero_prob );
    }else{
      if( Options_vec(0)==0 ) jnll_i(i) -= log( 1-zero_prob ) + dlognorm( y_i(i), linpred_i(i), logsd, true );
      if( Options_vec(0)==1 ) jnll_i(i) -= log( 1-zero_prob ) + dgamma( y_i(i), 1/pow(logsd,2), exp(linpred_i(i))*pow(logsd,2), true );
      if( Options_vec(0)==2 ) jnll_i(i) -= log( 1-zero_prob ) + dinvgauss( y_i(i), exp(linpred_i(i)), logsd, true );
    }
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
