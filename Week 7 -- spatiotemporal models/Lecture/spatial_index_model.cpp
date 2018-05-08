
#include <TMB.hpp>
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( n_s );
  DATA_INTEGER( n_t );

  DATA_VECTOR( a_s );  // Area associated with location s
  DATA_VECTOR( c_i );  // counts for observation i
  DATA_FACTOR( s_i );  // Random effect index for observation i
  DATA_FACTOR( t_i );  // Random effect index for observation i

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_tau_O );
  PARAMETER( ln_tau_E );
  PARAMETER( ln_kappa );

  // Random effects
  PARAMETER_VECTOR( omega_s );
  PARAMETER_ARRAY( epsilon_st );

  // Objective funcction
  using namespace density;
  int n_i = c_i.size();
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaO = 1 / sqrt(4 * M_PI * exp(2*ln_tau_O) * exp(2*ln_kappa));
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau_E) * exp(2*ln_kappa));

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  jnll_comp(1) += SCALE( GMRF(Q), 1/exp(ln_tau_O) )( omega_s );
  for( int t=0; t<n_t; t++){
    jnll_comp(2) += SCALE( GMRF(Q), 1/exp(ln_tau_E) )( epsilon_st.col(t) );
  }

  // True density and abundance
  array<Type> log_d_st( n_s, n_t );
  for( int t=0; t<n_t; t++){
  for( int s=0; s<n_s; s++){
    log_d_st(s,t) = beta0 + omega_s(s) + epsilon_st(s,t);
  }}

  // Probability of data conditional on random effects
  for( int i=0; i<n_i; i++){
    if( !isNA(c_i(i)) ) jnll_comp(0) -= dpois( c_i(i), exp(log_d_st(s_i(i),t_i(i))), true );
  }

  // Objective function
  Type jnll = jnll_comp.sum();

  // Derived quantities
  vector<Type> D_t( n_t );
  D_t.setZero();
  for( int t=0; t<n_t; t++){
  for( int s=0; s<n_s; s++){
    D_t(t) += a_s(s) * exp(log_d_st(s,t));
  }}

  // Reporting
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( SigmaO );
  REPORT( log_d_st );
  REPORT( D_t );
  ADREPORT( D_t );

  return jnll;
}
