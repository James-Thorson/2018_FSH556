
#include <TMB.hpp>
// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
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

  // Probability of random effects
  using namespace density;
  ///// Make precision matrix for single axis
  matrix<Type> Q_yy(n_y,n_y);
  Q_yy.setZero();
  for(int y=0; y<n_y; y++) Q_yy(y,y) = (1+pow(rho,2));
  for(int y=1; y<n_y; y++){
    Q_yy(y-1,y) = -rho;
    Q_yy(y,y-1) = -rho;
  }
  REPORT( Q_yy )
  //// Going downstream in the random effect vector
  if( Options_vec(0)==0 ){
    vector<Type> Tmp_y(n_y);
    matrix<Type> Q0_yy(n_y,n_y);
    Q0_yy = Q_yy * ( 1-pow(rho,2) ) / sigma2;
    matrix<Type> Q1_yy(n_y,n_y);
    Q1_yy = Q_yy / sigma2;
    for(int x=0; x<n_x; x++){
      for(int y=0; y<n_y; y++){
        if(x==0) Tmp_y(y) = epsilon_xy(0,y);
        if(x>=1) Tmp_y(y) = epsilon_xy(x,y) - rho*epsilon_xy(x-1,y);
      }
      if(x==0) jnll_comp(1) -= dmvnorm( Tmp_y, Q0_yy, true );
      if(x>=1) jnll_comp(1) -= dmvnorm( Tmp_y, Q1_yy, true );
    }
  }
  //// Calculate using precision matrix
  if( Options_vec(0)==1 ){
    int n_z = n_x * n_y;
    matrix<Type> Q_zz(n_z, n_z);
    Q_zz = kronecker( Q_yy, Q_yy );
    Q_zz = Q_zz / sigma2;
    REPORT( Q_zz );
    vector<Type> epsilon_z(n_z);
    int Count = 0;
    for( int x=0; x<n_x; x++){
    for( int y=0; y<n_y; y++){
      epsilon_z(Count) = epsilon_xy(x,y);
      Count++;
    }}
    jnll_comp(1) -= dmvnorm( epsilon_z, Q_zz, true );
    //jnll_comp(1) += GMRF(Q_zz)( epsilon_z );
  }
  //// Calculate using built-in TMB functions
  if( Options_vec(0)==3 ){
    jnll_comp(1) += SCALE( SEPARABLE(AR1(rho),AR1(rho)), pow(sigma2,0.5) / pow(1-pow(rho,2),0.5) / pow(1-pow(rho,2),0.5) )( epsilon_xy );      // Include "pow(1-pow(rho,2),0.5)" twice for 2D unit variance
  }
  // Sparse precision matrix using external computation of constituent matrices
  if( Options_vec(0)==4 ){
    Eigen::SparseMatrix<Type> Q_zz = ( M0*pow(1+pow(rho,2),2) + M1*(1+pow(rho,2))*(-rho) + M2*pow(rho,2) ) / sigma2;
    jnll_comp(1) += GMRF(Q_zz)( epsilon_xy );
    //Eigen::SparseMatrix<Type> Q_zz = ( M0*pow(1+pow(rho,2),2) + M1*(1+pow(rho,2))*(-rho) + M2*pow(rho,2) );
    //jnll_comp(1) += SCALE( GMRF(Q_zz), pow(sigma2,-0.5) )( epsilon_xy );
  }

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
