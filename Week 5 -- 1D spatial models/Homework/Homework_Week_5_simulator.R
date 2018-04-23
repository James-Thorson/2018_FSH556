
Sim_Fn = function( n_i=1000, Scale=2, Sigma2=1, logSD_spatial=0.1, L0=10, Linf_0=100, beta_y=0.02, growth_rate=0.1, mortality_rate=growth_rate*1.6, logSD_resid=0.05){
  require( RandomFields )

  # Sample locations
  n_i = 1000
  y_i = runif( n=n_i, min=32, max=49 )

  # Simulate spatial process
  RMmodel = RMgauss(var=logSD_spatial^2, scale=Scale)
  Linf_i = Linf_0 * exp( RFsimulate(model=RMmodel, x=rep(0,n_i), y=y_i)@data[,1] - Sigma2/2 ) * exp( beta_y*(y_i-40.5) )
  plot( y=Linf_i, x=y_i )

  # Simulate ages of samples
  Survival_a = exp( -mortality_rate * 1:100 )
  a_i = sample( size=n_i, x=1:100, prob=Survival_a, replace=TRUE)
  l_i = Linf_i - (Linf_i-L0) * exp( -growth_rate * a_i ) * exp( rnorm(n_i, mean=-logSD_resid^2/2, sd=logSD_resid) )
  plot( x=a_i, y=l_i )

  # Bundle and return stuff
  DF = data.frame( "y_i"=y_i, "a_i"=a_i, "l_i"=l_i, "Linf_i"=Linf_i)
  return(DF)
}

