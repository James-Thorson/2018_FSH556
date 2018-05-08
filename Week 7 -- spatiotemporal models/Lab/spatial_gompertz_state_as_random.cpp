// Space time 
#include <TMB.hpp>
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Indices
  DATA_INTEGER( n_i );         // Total number of observations
  DATA_INTEGER( n_x );         // Number of vertices in SPDE mesh
  DATA_INTEGER( n_t );          // Number of years
  DATA_INTEGER( n_p );          	// number of columns in covariate matrix X

  // Data
  DATA_IVECTOR( x_s );	      // Association of each station with a given vertex in SPDE mesh
  DATA_VECTOR( c_i );       	// Count data
  DATA_IVECTOR( s_i );        // Station for each sample
  DATA_IVECTOR( t_i );        // Time for each sample
  DATA_MATRIX( X_xp );		            // Covariate design matrix

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER_VECTOR(alpha);   // Mean of Gompertz-drift field
  PARAMETER(phi);            // Offset of beginning from equilibrium
  PARAMETER(log_tau_U);      // log-inverse SD of Epsilon
  PARAMETER(log_tau_O);      // log-inverse SD of Omega
  PARAMETER(log_kappa);      // Controls range of spatial variation
  PARAMETER(rho);             // Autocorrelation (i.e. density dependence)

  // Random effects
  PARAMETER_ARRAY(log_D_xt);  // Spatial process variation
  PARAMETER_VECTOR(Omega_input);   // Spatial variation in carrying capacity

  // objective function -- joint negative log-likelihood 
  using namespace density;
  Type jnll = 0;
  vector<Type> jnll_comp(3);
  jnll_comp.setZero();

  // Spatial parameters
  Type kappa2 = exp(2.0*log_kappa);
  Type kappa4 = kappa2*kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp( log_kappa );
  Type SigmaU = 1 / sqrt(4*pi*exp(2*log_tau_U)*exp(2*log_kappa));
  Type SigmaO = 1 / sqrt(4*pi*exp(2*log_tau_O)*exp(2*log_kappa));
  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;

  // Objects for derived values
  vector<Type> eta_x(n_x);
  array<Type> log_Dpred_xt(n_x, n_t);
  vector<Type> Omega_x(n_x);
  vector<Type> Equil_x(n_x);
 
  // Transform GMRFs
  eta_x = X_xp * alpha.matrix();
  for(int x=0; x<n_x; x++){
    Omega_x(x) = Omega_input(x);
    Equil_x(x) = ( eta_x(x) + Omega_x(x) ) / (1-rho);
  }
  
  // Calculate expectation for state-vector
  for(int x=0; x<n_x; x++){
  for(int t=0; t<n_t; t++){
    if(t==0) log_Dpred_xt(x,t) = phi + Equil_x(x);
    if(t>=1) log_Dpred_xt(x,t) = rho*log_D_xt(x,t-1) + eta_x(x) + Omega_x(x);
  }}
  
  // Probability of Gaussian-Markov random fields (GMRFs)
  // Omega
  jnll_comp(0) = SCALE(GMRF(Q), exp(-log_tau_O))(Omega_x);

  // Epsilon
  for(int t=0; t<n_t; t++){
    jnll_comp(1) += SCALE( GMRF(Q), exp(-log_tau_U) )( log_D_xt.col(t) - log_Dpred_xt.col(t) );
  }

  // total_abundance contribution from observations
  vector<Type> chat_i(n_i);
  for (int i=0; i<n_i; i++){
    chat_i(i) = exp( log_D_xt(x_s(s_i(i)),t_i(i)) );
    if( !isNA(c_i(i)) ){
      jnll_comp(2) -= dpois( c_i(i), chat_i(i), true );
    }
  }
  jnll = jnll_comp.sum();

  // Diagnostics
  REPORT( jnll_comp );
  REPORT( jnll );
  // Spatial field summaries
  REPORT( Range );
  REPORT( SigmaU );
  REPORT( SigmaO );
  REPORT( rho );
  ADREPORT( Range );
  ADREPORT( SigmaU );
  ADREPORT( SigmaO );
  // Fields
  REPORT( log_D_xt );
  REPORT( log_Dpred_xt );
  REPORT( Omega_x );
  REPORT( Equil_x );

  return jnll;
}
