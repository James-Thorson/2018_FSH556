// Space time 
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Settings
  DATA_IVECTOR( Options_vec );  // Slot 0: method;  0=conditiona;  1=joint

  // Data
  DATA_VECTOR(c_i);         // Count data
  DATA_FACTOR(b_i);         // Branch number b for each observationp

  // Inputs regarding branched network
  DATA_FACTOR(parent_b); // index of parent
  DATA_FACTOR(child_b); // index of child for branch b
  DATA_VECTOR(dist_b);  // distance to parent
  
  // Fixed effects
  PARAMETER(log_theta);             // Autocorrelation (i.e. density dependence)
  PARAMETER(log_SD);
  PARAMETER(log_mean);

  // Random effects
  PARAMETER_VECTOR(Epsiloninput_b);  // Spatial process variation

  // objective function -- joint negative log-likelihood
  int n_i = c_i.size();
  int n_b = Epsiloninput_b.size();
  Type jnll = 0; 

  // Derived parameters
  Type SDinput = exp(log_SD);
  Type theta = exp(log_theta);
  Type mu = exp(log_mean);
  Type SDmarg = SDinput / sqrt(2*theta);

  // Probability of GRF on network -- SPATIAL
  // Version 0: Sweep upstream to downstream through network
  if( Options_vec(0)==0 ){
    vector<Type> rho_b(n_b);
    vector<Type> SD_b(n_b);
    vector<Type> SDinput_b(n_b);
    for(int b=0; b<n_b; b++){
      if( isNA(dist_b(b)) ){
        // Correlation between i and parent(i) as distance -> INF
        rho_b(b) = 0;
        // SD of Ornstein-Uhlenbeck process as distance -> INF
        SDinput_b(b) = SDinput / pow(2*theta, 0.5);
        // conditional probability
        jnll -= dnorm(Epsiloninput_b(child_b(b)), Type(0.0), SDinput_b(b), true);
      }
      if( !isNA(dist_b(b)) ){
        // Correlation between i and parent(i)
        rho_b(b) = exp(-theta * dist_b(b));
        // SD of O-U process
        SDinput_b(b) = pow( pow(SDinput,2)/(2*theta) * (1-exp(-2*theta*dist_b(b))), 0.5 );
        // conditional probability
        jnll -= dnorm(Epsiloninput_b(child_b(b)), rho_b(b)*Epsiloninput_b(parent_b(b)), SDinput_b(b), true);
      }
    }
  }
    // Sweep upstream to downstream through network
  if( Options_vec(0)==1 ){
    Eigen::SparseMatrix<Type> Q( n_b, n_b );;
    for(int b=0; b<n_b; b++){
      Q.coeffRef( b, b ) = pow(SDmarg, -2);
    }
    for(int b=1; b<n_b; b++){
      Q.coeffRef( parent_b(b), child_b(b) ) = -exp(-theta*dist_b(b)) / (pow(SDmarg,2) * (1-exp(-2*theta*dist_b(b))));
      Q.coeffRef( child_b(b), parent_b(b) ) = Q.coeffRef( parent_b(b), child_b(b) );
      Q.coeffRef( parent_b(b), parent_b(b) ) += exp(-2*theta*dist_b(b)) / (pow(SDmarg,2) * (1-exp(-2*theta*dist_b(b))));
      Q.coeffRef( child_b(b), child_b(b) ) += exp(-2*theta*dist_b(b)) / (pow(SDmarg,2) * (1-exp(-2*theta*dist_b(b))));
    }
    REPORT( Q );
    jnll += density::GMRF(Q)(Epsiloninput_b);
  }

  // Probability of data given GRF
  vector<Type> lambda_i( n_i );
  for (int i=0; i<n_i; i++){
    lambda_i(i) = exp( log_mean + Epsiloninput_b(b_i(i)) );
    if( !isNA(c_i(i)) ){
      jnll -= dpois(c_i(i), lambda_i(i), true);
    }
  }

  REPORT( SDmarg );

  return jnll;
}
