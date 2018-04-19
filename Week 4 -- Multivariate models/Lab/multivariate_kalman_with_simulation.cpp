#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_MATRIX( y_tp );

  // Calculations
  int nt = y_tp.col(0).size();
  int np = y_tp.row(0).size();

  // Parameters
  PARAMETER( log_VarM );
  PARAMETER_VECTOR( x0_p );
  PARAMETER_MATRIX( Lp );
  PARAMETER_MATRIX( delta_tp );
  
  // Simulate random effects
  // See:  http://kaskr.github.io/adcomp/Simulation.html
  SIMULATE {
    for( int t=0; t<nt; t++){
    for( int p=0; p<np; p++){
      delta_tp(t,p) = rnorm( Type(0), Type(1) );
    }}
  }

  // Derived quantities
    // Covariance matrix
  matrix<Type> CovP_pp( np, np );
  CovP_pp = Lp * Lp.transpose();
    // Correlated errors
  matrix<Type> epsilon_tp( nt, np );
  epsilon_tp.setZero();
  for( int t=0; t<nt; t++){
  for( int p1=0; p1<np; p1++){
  for( int p2=0; p2<np; p2++){
    epsilon_tp(t,p1) += Lp(p1,p2) * delta_tp(t,p2);
  }}}
  matrix<Type> x_tp( nt, np );
    // Predicted
  for( int p=0; p<np; p++){
    x_tp(0,p) = x0_p(p) + epsilon_tp(0,p);
    for( int t=1; t<nt; t++){
      x_tp(t,p) = x_tp(t-1,p) + epsilon_tp(t,p);
    }
  }

  // Objective funcction
  Type jnll = 0;
  
  // Probability of random coefficients
  for( int t=0; t<nt; t++){
  for( int p=0; p<np; p++){
    jnll -= dnorm( delta_tp(t,p), Type(0), Type(1), true );
  }}

  // Probability of data conditional on fixed and random effect values
  for( int t=0; t<nt; t++){
  for( int p=0; p<np; p++){
    if( !isNA(y_tp(t,p)) ){
      jnll -= dnorm( y_tp(t,p), x_tp(t,p), exp(log_VarM), true );
    }
  }}
  
  // Simulate data
  SIMULATE {
    for( int t=0; t<nt; t++){
    for( int p=0; p<np; p++){
      y_tp(t,p) = rnorm( x_tp(t,p), exp(log_VarM) );
    }}
    REPORT( y_tp )
  }

  // Reporting
  Type VarM = exp(log_VarM);

  REPORT( CovP_pp );
  REPORT( VarM );
  REPORT( x_tp );
  REPORT( epsilon_tp );
  REPORT( Lp );

  ADREPORT( CovP_pp );
  ADREPORT( VarM );
  ADREPORT( x_tp );

  return jnll;
}
