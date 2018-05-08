
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
  DATA_MATRIX( c_xy );

  // Sparse matrices for forming the AR precision matrix (Version #4 below)
  // Q = M0*(1+rho^2)^2 + M1*(1+rho^2)*(-rho) + M2*rho^2
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_sigma2 );
  PARAMETER( logit_rho );

  // Random effects
  PARAMETER_ARRAY( epsilon_xy );

  // Objective funcction
  int n_x = c_xy.col(0).size();
  int n_y = c_xy.row(0).size();
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();
  Type sigma2 = exp(ln_sigma2);
  Type rho = 1 / (1 + exp(-logit_rho));

  //// Calculate using built-in TMB functions
  using namespace density;
  Eigen::SparseMatrix<Type> Q_zz = ( M0*pow(1+pow(rho,2),2) + M1*(1+pow(rho,2))*(-rho) + M2*pow(rho,2) ) / sigma2;
  jnll_comp(1) += GMRF(Q_zz)( epsilon_xy );

  // Probability of data conditional on random effects
  Type Total_Abundance = 0;
  for( int x=0; x<n_x; x++){
  for( int y=0; y<n_y; y++){
    if( !isNA(c_xy(x,y)) ) jnll_comp(0) -= dpois( c_xy(x,y), exp(beta0 + epsilon_xy(x,y)), true );
    Total_Abundance += exp(beta0 + epsilon_xy(x,y));
  }}

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( sigma2 );
  REPORT( rho );
  REPORT( Total_Abundance );
  ADREPORT( Total_Abundance );

  return jnll;
}
