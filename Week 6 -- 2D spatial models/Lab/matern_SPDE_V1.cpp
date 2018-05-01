
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
  using namespace density;

  // Data
  DATA_VECTOR( c_i );  // counts for observation i
  DATA_FACTOR( j_i );  // Random effect index for observation i

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_tau );
  PARAMETER( ln_kappa );

  // Random effects
  PARAMETER_VECTOR( epsilon_j );

  // Objective funcction
  int n_i = c_i.size();
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();

  // Derived quantities
  Type Range = sqrt(8) / exp( ln_kappa );
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau) * exp(2*ln_kappa));

  // Probability of random effects
  Eigen::SparseMatrix<Type> Q = exp(4*ln_kappa)*M0 + Type(2.0)*exp(2*ln_kappa)*M1 + M2;
  jnll_comp(1) += SCALE( GMRF(Q), 1/exp(ln_tau) )( epsilon_j );

  // Probability of data conditional on random effects
  Type Total_Abundance = 0;
  for( int i=0; i<n_i; i++){
    if( !isNA(c_i(i)) ) jnll_comp(0) -= dpois( c_i(i), exp(beta0 + epsilon_j(j_i(i))), true );
    Total_Abundance += exp(beta0 + epsilon_j(j_i(i)));
  }

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( Total_Abundance );
  ADREPORT( Total_Abundance );

  return jnll;
}
