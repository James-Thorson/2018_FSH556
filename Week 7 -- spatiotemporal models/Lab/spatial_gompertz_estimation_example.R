

setwd( "C:/Users/James.Thorson/Desktop/Project_git/2016_classes_private/Spatio-temporal models/Week 7 -- spatiotemporal models/Lab" )

#########################
# Spatial Gompertz model
# SEE: James T. Thorson, Hans Skaug, Kasper Kristensen, Andrew O. Shelton, Eric J. Ward, John Harms, Jim Benante. In press. The importance of spatial models for estimating the strength of density dependence. Ecology.
########################

# load libraries
library(INLA)
library(TMB)
library(RandomFields)
library(TMBdebug)

source( "Sim_Gompertz_Fn.R" )

# Read data
# n_years=10; n_stations=100; SpatialScale=0.1; SD_O=0.5; SD_E=0.2; SD_extra=0; rho=0.8; logMeanDens=1; phi=NULL; Loc=NULL
Sim_List = Sim_Gompertz_Fn( n_years=10, n_stations=100, SpatialScale=0.1, SD_O=0.4, SD_E=0.2, SD_extra=0, rho=0.5, logMeanDens=1, phi=0.0, Loc=NULL )
DF = Sim_List[["DF"]]
loc_xy_orig = loc_xy = Sim_List[["Loc"]]

# Reduce number of stations -- OPTIONAL
n_knots = 50
if( n_knots < nrow(loc_xy) ){
  library(RANN)
  knots_xy = kmeans( x=loc_xy, centers=n_knots )
  # Modify data
  loc_xy = knots_xy$centers
  DF[,'Site'] = knots_xy$cluster[DF[,'Site']]
}

plot( loc_xy_orig, cex=2, pch=20 )
points( loc_xy, cex=2, pch=20, col="red")

# Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
mesh = inla.mesh.create( loc_xy )
spde = inla.spde2.matern( mesh )

# display stations
#plot( x=loc_xy[,'x'], y=loc_xy[,'y'])


###################
# Parameter estimation
###################

#####  Version 0 -- Sweep upstream to downstream through time
Version = "spatial_gompertz_state_as_random"

# Compile
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

# Build inputs
X_xp = matrix( 1, ncol=1, nrow=mesh$n)
Data = list( n_i=nrow(DF), n_x=mesh$n, n_t=max(DF$Year), n_p=ncol(X_xp), x_s=mesh$idx$loc-1, c_i=DF[,'Simulated_example'], s_i=DF[,'Site']-1, t_i=DF[,'Year']-1, X_xp=X_xp, G0=spde$param.inla$M0, G1=spde$param.inla$M1, G2=spde$param.inla$M2)
Parameters = list(alpha=c(0.0), phi=0.0, log_tau_U=1.0, log_tau_O=1.0, log_kappa=0.0,	rho=0.5, log_D_xt=matrix(rnorm(mesh$n*Data$n_t),nrow=mesh$n,ncol=Data$n_t), Omega_input=rnorm(mesh$n))
Random = c("log_D_xt","Omega_input")

# Make object
obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=FALSE, DLL=Version)

# Run optimizer
start_time = Sys.time()
opt0 = nlminb(obj$par, objective=obj$fn, gradient=obj$gr, lower=c(rep(-20,2),rep(-10,3),-0.999), upper=c(rep(20,2),rep(10,3),0.999), control=list(eval.max=1e4, iter.max=1e4, trace=1))
opt0[["final_gradient"]] = obj$gr( opt0$par )
opt0[["total_time"]] = Sys.time() - start_time

# Get standard errors
Report0 = obj$report()
SD0 = try( sdreport(obj) )

#####  Version 3 -- Joint analysis using TMB functions
Version = "spatial_gompertz"

# Compile
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

# Build inputs
X_xp = matrix( 1, ncol=1, nrow=mesh$n)
Data = list( n_i=nrow(DF), n_x=mesh$n, n_t=max(DF$Year), n_p=ncol(X_xp), x_s=mesh$idx$loc-1, c_i=DF[,'Simulated_example'], s_i=DF[,'Site']-1, t_i=DF[,'Year']-1, X_xp=X_xp, G0=spde$param.inla$M0, G1=spde$param.inla$M1, G2=spde$param.inla$M2)
Parameters = list(alpha=c(0.0), phi=0.0, log_tau_E=1.0, log_tau_O=1.0, log_kappa=0.0,	rho=0.5, Epsilon_input=matrix(rnorm(mesh$n*Data$n_t),nrow=mesh$n,ncol=Data$n_t), Omega_input=rnorm(mesh$n))
Random = c("Epsilon_input","Omega_input")

# Make object
obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=FALSE, DLL=Version)

# Run optimizer
start_time = Sys.time()
opt3 = nlminb(obj$par, objective=obj$fn, gradient=obj$gr, lower=c(rep(-20,2),rep(-10,3),-0.999), upper=c(rep(20,2),rep(10,3),0.999), control=list(eval.max=1e4, iter.max=1e4, trace=1))
opt3[["final_gradient"]] = obj$gr( opt3$par )
opt3[["total_time"]] = Sys.time() - start_time

# Get standard errors
Report3 = obj$report()
SD3 = try( sdreport(obj) )

######## Compare results

# Report
unlist( Report0[c('Range','SigmaO','SigmaU','SigmaE','rho')] )
unlist( Report3[c('Range','SigmaO','SigmaU','SigmaE','rho')] )

