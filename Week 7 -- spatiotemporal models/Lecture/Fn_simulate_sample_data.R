
  # logmean=1; Scale=0.2; SD_omega=1; SD_epsilon=1; n_per_year=100; n_years=10
Sim_Fn = function( logmean=1, Scale=0.2, SD_omega=1, SD_epsilon=1, n_per_year=100, n_years=10 ){

  loc_xy = cbind( "x"=runif(n_per_year), "y"=runif(n_per_year))

  RF_omega = RMgauss(var=SD_omega^2, scale=Scale)
  RF_epsilon = RMgauss(var=SD_epsilon^2, scale=Scale)

  log_d_rt = Epsilon_rt = matrix(NA, ncol=n_years, nrow=n_per_year)
  Omega_r = RFsimulate(model=RF_omega, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1]
  for(t in 1:n_years){
    Epsilon_rt[,t] = RFsimulate(model=RF_epsilon, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1]
    log_d_rt[,t] = logmean + Epsilon_rt[,t] + Omega_r
  }
  b_t = colSums( exp(log_d_rt) )

  # Simulate samples for each site and year
  DF = expand.grid("s_i"=1:n_per_year, "t_i"=1:10)
  DF = cbind( DF, "log_d_rt"=log_d_rt[ as.matrix(DF[,c('s_i','t_i')]) ] )
  DF = cbind( DF, "c_i"=rpois(n_per_year*n_years, lambda=exp(DF[,'log_d_rt'])) )

  # Return stuff
  Return = list("DF"=DF, "Epsilon_rt"=Epsilon_rt, "Omega_r"=Omega_r, "loc_xy"=loc_xy, "b_t"=b_t, "n_per_year"=n_per_year, "n_years"=n_years)
  return( Return )
}
