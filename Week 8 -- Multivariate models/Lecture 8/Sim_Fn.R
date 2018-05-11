Sim_Fn <-
function( n_p, n_s=200, n_f=2, SpatialScale=0.1, SD_O=1.0, logMeanDens=3.0, Loadings_pf=NULL, Loc=NULL ){
  # n_p: number of species
  # n_f: number of factors
  # n_s: Number of sites

  # Loadings matrix
  if( is.null(Loadings_pf) ){
    Loadings_pf = matrix( rnorm(n_f*n_p), nrow=n_p, ncol=n_f)
    for(fI in 1:n_f) Loadings_pf[seq(from=1,to=fI-1,length=fI-1),fI] = 0
  }

  # Species intercepts
  Beta_p = rep(logMeanDens, n_p)

  # Locations for sampling
  if( is.null(Loc) ) Loc_xy = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )

  # Initialize spatial model
  model_O <- RMgauss(var=SD_O^2, scale=SpatialScale)

  # Simulate fields
  Omega_sf = matrix(NA, ncol=n_f, nrow=n_s)
  for(fI in 1:n_f){
    Omega_sf[,fI] = RFsimulate(model=model_O, x=Loc_xy[,'x'], y=Loc_xy[,'y'])@data[,1]
  }

  # Calculate linear predictor
  ln_Yexp_sp = Omega_sf%*%t(Loadings_pf) + outer(rep(1,n_s),Beta_p)

  # Simulate data
  Y_sp = matrix(rpois( n_s*n_p, lambda=exp(ln_Yexp_sp) ), ncol=n_p, byrow=FALSE)

  # Create dummy matrix of covariates
  X_sj = cbind( rep(1,nrow(Y_sp)) )

  # Return stuff
  Sim_List = list("Y_sp"=Y_sp, "X_sj"=X_sj, "Loadings_pf"=Loadings_pf, "Loc_xy"=Loc_xy, "Omega_sf"=Omega_sf, "ln_Yexp_sp"=ln_Yexp_sp)
  return(Sim_List)
}
