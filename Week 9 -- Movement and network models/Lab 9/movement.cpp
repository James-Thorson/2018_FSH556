#include <TMB.hpp>
// Space time

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Options
  DATA_FACTOR( Options_vec );
  
  // Dimensions
  DATA_INTEGER(n_i);         // Number of observations (stacked across all years)
  DATA_INTEGER(n_r);        // Number of triangles
  DATA_INTEGER(n_g);        // Number of vertices in GMRF approximation to triangle values
  DATA_INTEGER(n_t);         // Number of real "strata" (i.e., verticies containing data)
  DATA_INTEGER(n_tdiv);

  // Data vectors
  DATA_VECTOR( c_i );
  DATA_FACTOR( r_i );
  DATA_FACTOR( t_i );

  // Movement matrices
  DATA_SPARSE_MATRIX( M1 );
  DATA_SPARSE_MATRIX( M2 );
  DATA_SPARSE_MATRIX( M3 );
  DATA_SPARSE_MATRIX( M4 );

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER( log_beta );
  PARAMETER( alpha );
  PARAMETER_VECTOR( log_mvec );   // Movement parameters (east, north, west, south)
  PARAMETER( logkappa );
  PARAMETER( logSigmaU );
  PARAMETER( logSigmaO );
  PARAMETER( logmeanu0 );

  // -- Gaussian random fields
  PARAMETER_ARRAY( ln_u_gt );            // State-matrix
  PARAMETER_VECTOR( Omegainput_g );      // Spatial variation in growth rates

  // Global parameters
  using namespace density;
  Type pi = 3.141592;
  Type beta = exp( log_beta );
  Type rho = 1.0 - beta; 
  
  // Objective function
  vector<Type> NLL_c(3);    // 0:ln_u; 1:Epsilon; 2:Likelihood
  NLL_c.setZero();
  vector<Type> NLL_u_t(n_t);
  NLL_u_t.setZero();
  Type NLL = 0;                // Objective function
  vector<Type> NLL_i(n_i);
  
  // Derived parameters
  Type SigmaU = exp(logSigmaU);
  Type SigmaO = exp(logSigmaO);
  Type kappa_pow2 = exp(2.0*logkappa);
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  Type tauU = 1 / sqrt(4*pi*exp(2*logSigmaU)*exp(2*logkappa));
  Type tauO = 1 / sqrt(4*pi*exp(2*logSigmaO)*exp(2*logkappa));
  Type Range_raw = sqrt(8) / exp( logkappa );
  vector<Type> mvec(4);
  mvec = exp( log_mvec );
 
  // Movement matrix
  Eigen::SparseMatrix<Type> M_sparse(n_r,n_r);
  Eigen::SparseMatrix<Type> Mdiv_sparse(n_r,n_r);
  Eigen::SparseMatrix<Type> Identity_sparse(n_r,n_r);
  if( Options_vec(0)==0 ){
    M_sparse.setIdentity();
    Identity_sparse.setIdentity();
    Mdiv_sparse.setIdentity();
  }
  if( Options_vec(0)==1 | Options_vec(0)==2 ){
    M_sparse = mvec(0)*M1 + mvec(1)*M2 + mvec(2)*M3 + mvec(3)*M4;
    Identity_sparse.setIdentity();
    Mdiv_sparse = M_sparse / Type(n_tdiv) + Identity_sparse; // Euler approximation
  }
  
//   Random field probability
  Eigen::SparseMatrix<Type> Q;
  Q = kappa_pow4*G0 + Type(2.0)*kappa_pow2*G1 + G2;
  GMRF_t<Type> gmrf_nll = GMRF(Q);

  // Derived fields
  vector<Type> Omega_g(n_g);
  Omega_g = Omegainput_g / tauO;

  // Track estimates and predictions of effort (dimensions n_r x n_t, i.e., density for each triangle in the domain)
  array<Type> u_rt(n_r,n_t);
  array<Type> upred_rt(n_r,n_t);
  // Track estimates and predictions of effort (dimensions n_gmrf x n_t, i.e., padded with zeros)
  array<Type> ln_uhat_gt(n_g,n_t);

  // Predict dynamics of state-vector
  for(int r=0; r<n_r; r++){
    ln_uhat_gt(r,0) = logmeanu0 + Omega_g(r);
    u_rt(r,0) = exp( ln_u_gt(r,0) );
  }
  for(int g=n_r; g<n_g; g++) ln_uhat_gt(g,0) = 0.0;
  for(int t=1; t<n_t; t++){
    for(int r=0; r<n_r; r++){
      u_rt(r,t) = exp( ln_u_gt(r,t) );
      upred_rt(r,t) = u_rt(r,t-1);
    }
    for(int tdev=0; tdev<n_tdiv; tdev++){
      upred_rt.col(t) = Mdiv_sparse * upred_rt.col(t).matrix();
    }
    for(int r=0; r<n_r; r++){
      if(Options_vec(1)==0) ln_uhat_gt(r,t) = log( upred_rt(r,t) * exp(alpha - beta*log(upred_rt(r,t)) + Omega_g(r)) );
      if(Options_vec(1)==1) ln_uhat_gt(r,t) = log( upred_rt(r,t) * exp(alpha - beta*upred_rt(r,t) + Omega_g(r)) );
    }
    for(int g=n_r; g<n_g; g++) ln_uhat_gt(g,t) = 0.0;
  }

  // Probability
  for(int t=0; t<n_t; t++){
    NLL_u_t(t) += SCALE(gmrf_nll, 1.0/tauU)( ln_u_gt.col(t)-ln_uhat_gt.col(t) );
  }
  NLL_c(0) = NLL_u_t.sum();

  // Epsilon probability
  NLL_c(1) = gmrf_nll(Omegainput_g);

  // Likelihood
  vector<Type> chat_i(n_i);
  for(int i=0; i<n_i; i++){
    // Calculate deviance
    chat_i(i) = exp(ln_u_gt(r_i(i),t_i(i)));
    NLL_i(i) = -1 * dpois( c_i(i), chat_i(i), true);
  }
  NLL_c(2) = NLL_i.sum();
  NLL = NLL_c.sum();

  // Diagnostic output
  REPORT( chat_i );
  REPORT( NLL_u_t );
  REPORT( u_rt );
  REPORT( ln_u_gt );
  REPORT( upred_rt );
  REPORT( ln_uhat_gt );
  REPORT( NLL_i );
  REPORT( Omega_g );
  REPORT( Range_raw );
  REPORT( SigmaU );
  REPORT( SigmaO );
  REPORT( tauU );
  REPORT( tauO );
  REPORT( NLL_c );
  REPORT( NLL );
  REPORT( rho );
  REPORT( log_mvec );
  REPORT( mvec );

  // Movement stuff
  REPORT( M_sparse );
  REPORT( Mdiv_sparse );

  return NLL;
}
