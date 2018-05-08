Sim_Gompertz_Fn <-
function( n_years, n_stations=100, phi=NULL, SpatialScale=0.1, SD_O=0.5, SD_E=0.2, SD_extra=0.1, rho=0.8, logMeanDens=1, Loc=NULL ){
  # Parameters
  if( is.null(phi) ) phi = rnorm(1, mean=0, sd=1)
  alpha = logMeanDens * (1-rho)

  # Spatial model
  if( is.null(Loc) ) Loc = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )
  model_O <- RMgauss(var=SD_O^2, scale=SpatialScale)
  model_E <- RMgauss(var=SD_E^2, scale=SpatialScale)

  # Simulate Omega
  Omega = RFsimulate(model = model_O, x=Loc[,'x'], y=Loc[,'y'])@data[,1]

  # Simulate Epsilon
  Epsilon = array(NA, dim=c(n_stations,n_years))
  for(t in 1:n_years){
    Epsilon[,t] = RFsimulate(model=model_E, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
  }

  # Calculate Psi
  Theta = array(NA, dim=c(n_stations,n_years))
  for(t in 1:n_years){
    if(t==1) Theta[,t] = phi + Epsilon[,t] + (alpha + Omega)/(1-rho)
    if(t>=2) Theta[,t] = rho * Theta[,t-1] + alpha + Omega + Epsilon[,t]
  }

  # Simulate data
  DF = NULL
  for(s in 1:n_stations){
  for(t in 1:n_years){
    Tmp = c("Site"=s, "Year"=t, "Simulated_example"=rpois(1,lambda=exp(Theta[s,t]+SD_extra*rnorm(1))))#, 'Lon..DDD.DDDDD.'=Loc[s,1], 'Lat..DD.DDDDD.'=Loc[s,2] )
    DF = rbind(DF, Tmp)
  }}
  DF = cbind( DF, 'Longitude'=Loc[DF[,'Site'],1], 'Latitude'=Loc[DF[,'Site'],2] )
  DF = data.frame(DF, row.names=NULL)

  # Return stuff
  Sim_List = list("DF"=DF, "phi"=phi, "Loc"=Loc, "Omega"=Omega, "Epsilon"=Epsilon, "Theta"=Theta)
  return(Sim_List)
}
