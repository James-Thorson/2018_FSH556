// Space time
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_SCALAR( dummy );

  // Parameters
  PARAMETER_VECTOR( Params );

  // Objective funcction
  Type jnll = pow(1-Params(0),2) + 100*pow(Params(1)-pow(Params(0),2),2);

  // Reporting
  return jnll;
}
