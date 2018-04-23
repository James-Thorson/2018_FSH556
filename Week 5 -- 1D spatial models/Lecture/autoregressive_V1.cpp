
#include <TMB.hpp>

// trace of a matrix
template<class Type>
Type trace( matrix<Type> mat ){
  Type Return = 0;
  for(int i=0; i<mat.col(0).size(); i++) Return += mat(i,i);
  return Return;
}

template<class Type>
Type logdet( matrix<Type> mat ){
  int dim = mat.col(0).size();
  matrix<Type> chol(dim,dim);
  chol = mat.llt().matrixL();   /* L0 L0' = Q0 */
  Type logdet_mat = 0;
  for(int i=0; i<dim; i++ ) logdet_mat += Type(2.0) * log(chol(i,i));
  return logdet_mat;
}

template<class Type>
Type dmvnorm( vector<Type> x, matrix<Type> Q, int give_log=0 ){
  int n_x = x.size();
  Type logres = 0;
  vector<Type> Tmp_x(n_x);
  Type logdet_Q = logdet( Q );
  Tmp_x =  Q * x.matrix();
  logres += ( Type(0.5)*logdet_Q );
  logres += Type(-0.5) * (x * Tmp_x).sum();  //
  if (give_log) return logres; else return exp(logres);
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Settings
  DATA_VECTOR( Options_vec );
  // Slot 0: method for calculating probability of random effects

  // Data
  DATA_MATRIX( c_si );

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( ln_sigma2 );
  PARAMETER( logit_rho );

  // Random effects
  PARAMETER_VECTOR( epsilon_s );

  // Objective funcction
  int n_s = c_si.col(0).size();
  int n_i = c_si.row(0).size();
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();
  Type sigma2 = exp(ln_sigma2);
  Type rho = 1 / (1 + exp(-logit_rho));

  // Probability of random effects
  using namespace density;
  //// Going downstream in the random effect vector
  if( Options_vec(0)==0 ){
    jnll_comp(1) -= dnorm( epsilon_s(0), Type(0.0), pow(sigma2,0.5), true );
    for(int s=1; s<n_s; s++) jnll_comp(1) -= dnorm( epsilon_s(s), rho*epsilon_s(s-1), pow(sigma2,0.5), true );
  }
  //// Calculate using precision matrix
  if( Options_vec(0)==1 ){
    matrix<Type> Q_ss(n_s,n_s);
    Q_ss.setZero();
    for(int s=0; s<n_s; s++) Q_ss(s,s) = (1+pow(rho,2))/sigma2;
    for(int s=1; s<n_s; s++){
      Q_ss(s-1,s) = -rho/sigma2;
      Q_ss(s,s-1) = -rho/sigma2;
    }
    REPORT( Q_ss )
    jnll_comp(1) -= dmvnorm( epsilon_s, Q_ss, true );
  }
  //// Calculate using covariance matrix
  if( Options_vec(0)==2 ){
    matrix<Type> Cov_ss(n_s,n_s);
    for(int s1=0; s1<n_s; s1++){
    for(int s2=s1; s2<n_s; s2++){
      Cov_ss(s1,s2) = sigma2 / (1-pow(rho,2)) * pow( rho, double(s2-s1) );
      if(s1!=s2) Cov_ss(s2,s1) = Cov_ss(s1,s2);
    }}
    REPORT( Cov_ss );
    jnll_comp(1) += MVNORM(Cov_ss)( epsilon_s );
  }
  //// Calculate using built-in TMB functions
  if( Options_vec(0)==3 ){
    jnll_comp(1) += SCALE( AR1(rho), pow(sigma2 / (1-pow(rho,2)),0.5))( epsilon_s );
  }

  // Probability of data conditional on random effects
  for( int s=0; s<n_s; s++){
  for( int i=0; i<n_i; i++){
    jnll_comp(0) -= dpois( c_si(s,i), exp(beta0 + epsilon_s(s)), true );
  }}

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( sigma2 );
  REPORT( rho );

  return jnll;
}
