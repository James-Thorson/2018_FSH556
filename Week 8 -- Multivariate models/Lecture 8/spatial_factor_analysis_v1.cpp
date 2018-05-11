
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
  DATA_MATRIX(Y_sp);       	// Responses
  DATA_INTEGER(n_f);        // Number of factors
  DATA_INTEGER(n_x); // Number of vertices in mesh
  DATA_IVECTOR(x_s);	// Convert site s to vertex x
  DATA_MATRIX(X_sj);		// Covariate design matrix

  // SPDE objects
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters
  PARAMETER_MATRIX(beta_jp);
  PARAMETER_VECTOR(Loadings_vec);
  PARAMETER(log_kappa);

  // Random effect
  PARAMETER_ARRAY(Omega_xf);

  //
  using namespace density;
  int n_p = Y_sp.row(0).size();    // n_p is number of species
  int n_s = Y_sp.col(0).size();
  int n_j = X_sj.row(0).size();
  Type jnll = 0;

  // Unpack loadings matrix
  matrix<Type> Loadings_pf(n_p, n_f);
  int Count = 0;
  for(int f=0; f<n_f; f++){
    for(int p=0; p<n_p; p++){
      if(p>=f){
        Loadings_pf(p,f) = Loadings_vec(Count);
        Count++;
      }else{
        Loadings_pf(p,f) = 0.0;
      }
    }
  }

  // Spatial variables
  Type log_tau = log( 1 / (exp(log_kappa) * sqrt(4*M_PI)) );  // Ensures that MargSD = 1
  Type Range = sqrt(8) / exp( log_kappa );
  Eigen::SparseMatrix<Type> Q = exp(log_kappa*4)*M0 + Type(2.0)*exp(log_kappa*2)*M1 + M2;
  for(int f=0; f<n_f; f++){
    jnll += SCALE( GMRF(Q), 1/exp(log_tau))( Omega_xf.col(f) );
  }
  
  // Likelihood contribution from observations
  matrix<Type> ln_yexp_sp( n_s, n_p );
  matrix<Type> eta_sp( n_s, n_p );
  eta_sp = X_sj * beta_jp;
  for(int s=0; s<n_s; s++){
  for(int p=0; p<n_p; p++){
    // Predictor
    ln_yexp_sp(s,p) = eta_sp(s,p);
    for(int f=0; f<n_f; f++){
      ln_yexp_sp(s,p) += Omega_xf(x_s(s),f) * Loadings_pf(p,f);
    }
    // Likelihood
    if( !isNA(Y_sp(s,p)) ) jnll -= dpois( Y_sp(s,p), exp(ln_yexp_sp(s,p)), true );
  }}
  
  // Reporting
  REPORT( log_kappa );
  REPORT( log_tau );
  REPORT( Range );
  REPORT( Omega_xf );
  REPORT( Loadings_pf );
  REPORT( jnll );
  REPORT( ln_yexp_sp );

  return jnll;
}
