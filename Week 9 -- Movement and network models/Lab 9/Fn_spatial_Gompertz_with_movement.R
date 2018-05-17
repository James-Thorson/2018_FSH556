
Sim_Fn = function( MoveMat, SD_omega=1, SD_epsilon=1, SD_effort=1, Scale, Dynamical_Model, n_s, n_t, r_s, n_r, loc_r, alpha, beta ){

  # Simulate
  RF_omega = RMgauss(var=SD_omega^2, scale=Scale)
  RF_epsilon = RMgauss(var=SD_epsilon^2, scale=Scale)

  # Simulate density for each triangle
  # Gompertz: u(t+1) = u(t) * exp( alpha - beta*log(u(t)) )
  # Moran-Ricker: u(t+1) = u(t) * exp( alpha - beta*u(t) )
  upred_rt = u_rt = Epsilon_rt = matrix(NA, ncol=n_t, nrow=n_r)
  Omega_r = RFsimulate(model=RF_omega, x=loc_r[,1], y=loc_r[,2])@data[,1]
  for(t in 1:n_t){
    Epsilon_rt[,t] = RFsimulate(model=RF_epsilon, x=loc_r[,1], y=loc_r[,2])@data[,1]
    if(t==1) u_rt[,t] = exp( logmeanu0 + Epsilon_rt[,t] + Omega_r )
    if(t>=2){
      upred_rt[,t] = as.vector( MoveMat %*% u_rt[,t-1] )
      if( Dynamical_Model=="Gompertz" ) u_rt[,t] = upred_rt[,t] * exp(alpha - beta*log(upred_rt[,t]) + Omega_r + Epsilon_rt[,t])
      if( Dynamical_Model=="Ricker" ) u_rt[,t] = upred_rt[,t] * exp(alpha - beta*log(upred_rt[,t]) + Omega_r + Epsilon_rt[,t])
    }
  }

  # Simulate samples for each site and year
  DF = expand.grid("s_i"=1:n_s, "t_i"=1:n_t)
  DF = cbind( DF, "r_i"=r_s[DF[,'s_i']] )
  DF = cbind( DF, "cexp_i"=u_rt[ as.matrix(DF[,c('r_i','t_i')]) ] )
  DF = cbind( DF, "c_i"=rpois(n_s*n_t, lambda=DF[,'cexp_i']) )

  # Return stuff
  Return = list("DF"=DF, "MoveMat"=MoveMat, "upred_rt"=upred_rt, "u_rt"=u_rt, "Epsilon_rt"=Epsilon_rt, "Omega_r"=Omega_r)
  return( Return )
}

