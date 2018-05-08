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
  PARAMETER(log_tau_E);      // log-inverse SD of Epsilon
  PARAMETER(log_tau_O);      // log-inverse SD of Omega
  PARAMETER(log_kappa);      // Controls range of spatial variation
  PARAMETER(rho);             // Autocorrelation (i.e. density dependence)

  // Random effects
  PARAMETER_ARRAY(Epsilon_input);  // Spatial process variation
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
  Type SigmaE = 1 / sqrt(4*pi*exp(2*log_tau_E)*exp(2*log_kappa));
  Type SigmaO = 1 / sqrt(4*pi*exp(2*log_tau_O)*exp(2*log_kappa));
  Eigen::SparseMatrix<Type> Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;

  // Objects for derived values
  vector<Type> eta_x(n_x);
  array<Type> log_Dpred_xt(n_x, n_t);
  vector<Type> Omega_x(n_x);
  vector<Type> Equil_x(n_x);
  matrix<Type> Epsilon_xt(n_x, n_t);

  // Probability of Gaussian-Markov random fields (GMRFs)
  jnll_comp(0) += GMRF(Q)(Omega_input);
  jnll_comp(1) = SEPARABLE(AR1(rho),GMRF(Q))(Epsilon_input);

  // Transform GMRFs
  eta_x = X_xp * alpha.matrix();
  for(int x=0; x<n_x; x++){
    Omega_x(x) = Omega_input(x) / exp(log_tau_O);
    Equil_x(x) = ( eta_x(x) + Omega_x(x) ) / (1-rho);
    for( int t=0; t<n_t; t++){
      Epsilon_xt(x,t) = Epsilon_input(x,t) / exp(log_tau_E);
    }
  }

  // Likelihood contribution from observations
  vector<Type> log_chat_i(n_i);
  for (int i=0; i<n_i; i++){
    log_chat_i(i) = phi*pow(rho,t_i(i)) + Epsilon_xt(x_s(s_i(i)),t_i(i)) + (eta_x(x_s(s_i(i))) + Omega_x(x_s(s_i(i))) ) / (1-rho);
    if( !isNA(c_i(i)) ){
      jnll_comp(2) -= dpois( c_i(i), exp(log_chat_i(i)), true );
    }
  }
  jnll = jnll_comp.sum();

  // Diagnostics
  REPORT( jnll_comp );
  REPORT( jnll );
  // Spatial field summaries
  REPORT( Range );
  REPORT( SigmaE );
  REPORT( SigmaO );
  REPORT( rho );
  ADREPORT( Range );
  ADREPORT( SigmaE );
  ADREPORT( SigmaO );
  // Fields
  REPORT( Epsilon_xt );
  REPORT( Omega_x );
  REPORT( Equil_x );

  return jnll;
}
