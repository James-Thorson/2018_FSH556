

setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 7 -- spatiotemporal models/Lab" )

#########################
# Spatial Gompertz model
# SEE: James T. Thorson, Hans Skaug, Kasper Kristensen, Andrew O. Shelton, Eric J. Ward, John Harms, Jim Benante. In press. The importance of spatial models for estimating the strength of density dependence. Ecology.
########################

# load libraries
library(INLA)
library(TMB)
library(RandomFields)
library(raster)
library(RANN)

source( "Sim_Gompertz_Fn.R" )

# Read data
set.seed( 2 )
Sim_List = Sim_Gompertz_Fn( n_years=10, n_stations=1000, SpatialScale=0.1, SD_O=0.4, SD_E=1, SD_extra=0, rho=0.5, logMeanDens=1, phi=-2, Loc=NULL )
DF = Sim_List[["DF"]]
loc_xy_orig = loc_xy = Sim_List[["Loc"]]

# Reduce sample sizes to 100 per year
Which2Keep = sample(1:nrow(DF), size=100*Sim_List$n_years, replace=FALSE)
Which2Drop = setdiff(1:nrow(DF),Which2Keep)
DF[Which2Drop,'Simulated_example'] = NA

# Reduce number of stations -- OPTIONAL
n_knots = 50
if( n_knots < nrow(loc_xy) ){
  knots_xy = kmeans( x=loc_xy_orig, centers=n_knots )
  # Modify data
  loc_xy = knots_xy$centers
  DF[,'Site'] = knots_xy$cluster[DF[,'Site']]
}

# Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
mesh = inla.mesh.create( loc_xy, refine=TRUE, extend=-0.5 )
spde = inla.spde2.matern( mesh )

# Visualize mesh and predictive process
plot(mesh)
points( loc_xy_orig, cex=1.5, pch=20 )
points( loc_xy, cex=2, pch=3, col="green", lwd=5)

# Generate grid to visualize density
vizloc_xy = expand.grid( x=seq(0,1,by=0.001), y=seq(0,1,by=0.001) )
knots_xy = nn2( data=loc_xy_orig, query=vizloc_xy, k=1 )

# Plot densities
par( mfrow=c(2,5), mar=c(2,2,2,0), mgp=c(1.5,0.25,0) )
for( tI in 1:Sim_List$n_years ){
  vizTheta_xy = array(Sim_List$Theta[ cbind(knots_xy$nn.idx,tI) ], dim=c(1001,1001) )
  rasterTheta_xy = raster( vizTheta_xy )
  plot( rasterTheta_xy, xlim=c(0,1), ylim=c(0,1), main=paste0("Year ",tI) )
}


###################
#
# Parameter estimation
#
###################

#####################
#  Version 0 -- Sweep upstream to downstream through time "State-space parameterization"
#####################

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
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=FALSE, DLL=Version)

# Run optimizer
Opt0 = TMBhelper::Optimize( obj=Obj, lower=c(rep(-Inf,5),-0.999), upper=c(rep(Inf,5),0.999), getsd=TRUE, newtonsteps=1 )

# Get standard errors
Report0 = Obj$report()
H0 = Obj$env$spHess()

##################
#  Version 3 -- Joint analysis using TMB functions  "Innovations parameterization"
##################

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
Obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=FALSE, DLL=Version)

# Run optimizer
Opt3 = TMBhelper::Optimize( obj=Obj, lower=c(rep(-Inf,5),-0.999), upper=c(rep(Inf,5),0.999), getsd=TRUE, newtonsteps=1 )

# Get standard errors
Report3 = Obj$report()
H3 = Obj$env$spHess()

######## Compare results

# Check parameter estimates
unlist(Report0[c('Range','SigmaO','SigmaU','rho','phi')])
unlist(Report3[c('Range','SigmaO','SigmaE','rho','phi')])
Sim_List[["Parameters"]][c('SpatialScale','SigmaO','SigmaE','rho','phi')]

# Compare sparseness
sum( H0!=0 ) / prod(dim(H0))
sum( H3!=0 ) / prod(dim(H0))

# Show inner hessian
image( H0, main="Version0: Sweep through time" );
dev.new(); image(H3, main="Version3: TMB functions")

# Run times
Opt0$run_time
Opt3$run_time

