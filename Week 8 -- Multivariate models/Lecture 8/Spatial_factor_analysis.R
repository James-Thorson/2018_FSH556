
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 8 -- Multivariate models/Lecture 8/" )

# Load libraries
library(TMB)
library(INLA)
library(RandomFields)

#########################
# Simulation example
#########################

# Source
source( "Sim_Fn.R")

# Simulation settings
n_factors_true = 6
n_species = 10
n_stations = 100

# Estimation settings
n_factors_estimation = 2

# Simulate data
Data_List = Sim_Fn( n_p=n_species, n_s=n_stations, n_f=n_factors_true )
Y_sp = Data_List[["Y_sp"]]
X_sj = Data_List[["X_sj"]]
Loc_xy = Data_List[["Loc_xy"]]

# Create SPDE mesh
mesh = inla.mesh.create( Loc_xy )
spde = inla.spde2.matern( mesh )

# Data
Data = list("Y_sp"=Y_sp, "n_f"=n_factors_estimation, "n_x"=mesh$n, "x_s"=mesh$idx$loc-1, "X_sj"=X_sj, "M0"=spde$param.inla$M0, "M1"=spde$param.inla$M1, "M2"=spde$param.inla$M2)
# Parameters
Params = list( "beta_jp"=matrix(0,nrow=ncol(Data$X_sj),ncol=ncol(Data$Y_sp)), "Loadings_vec"=rep(1,Data$n_f*ncol(Data$Y_sp)-Data$n_f*(Data$n_f-1)/2), "log_kappa"=log(1), "Omega_xf"=matrix(0,nrow=mesh$n,ncol=Data$n_f) )
# Declare random
Random = c("Omega_xf")

# Compile if necessary
Version = "spatial_factor_analysis_v1"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )                                                         # log_tau=0.0,

# Initialization
Obj <- MakeADFun(data=Data, parameters=Params, random=Random, hessian=FALSE, inner.control=list(maxit=1000) )
table(names(Obj$env$last.par))

# Run model
Opt = TMBhelper::Optimize( obj=Obj, getsd=TRUE, newtonsteps=1 )

# Summarize
Report = Obj$report()

# Compare with simulated values
# Psi (loadings matrix)
Report$Loadings_pf
Data_List$Loadings_pf

# Compare covariances
plot( x=Report$Loadings_pf%*%t(Report$Loadings_pf), y=Data_List$Loadings_pf%*%t(Data_List$Loadings_pf), xlab="Estimated covariance", ylab="True covariance" )
abline( a=0, b=1, lty="dotted", lwd=2 )

# Omega (prediction for each species)
par( mfrow=c(1,min(n_factors_true,n_factors_estimation)) )
for(fI in 1:min(n_factors_true,n_factors_estimation)){
  plot(y=Report$Omega_xf[1:n_stations,fI], x=Data_List$Omega_sf[,fI], xlab="True", ylab="Estimated", main=paste("Factor",fI))
}

##########################
# Real dataset
##########################

# Load and format data
load( "EBS_Nspecies=20.RData" )
DF = DF[which( DF$year == 1990 ),]
SpeciesSet = unique( DF$spp )

# Identify unique locations
Match = match(unique(DF$TowID), DF$TowID)
Loc_xy = DF[Match,c('long','lat')]

# Make mesh
mesh = inla.mesh.create( Loc_xy )
spde = inla.spde2.matern( mesh )

# Format data frame
Y_sp = matrix( NA, nrow=length(Match), ncol=length(SpeciesSet), dimnames=list(NULL,SpeciesSet))
for(pI in 1:ncol(Y_sp)){
  DF_subset = DF[which(DF$spp==SpeciesSet[pI]),]
  Match = match( unique(DF$TowID), DF_subset[,'TowID'] )
  Y_sp[,pI] = round( DF_subset[Match,'catch'] )
}

# Data
Data = list("Y_sp"=Y_sp, "n_f"=3, "n_x"=mesh$n, "x_s"=mesh$idx$loc-1, "X_sj"=matrix(1,nrow=nrow(Y_sp),ncol=1), "M0"=spde$param.inla$M0, "M1"=spde$param.inla$M1, "M2"=spde$param.inla$M2)
# Parameters
Params = list( "beta_jp"=matrix(0,nrow=ncol(Data$X_sj),ncol=ncol(Data$Y_sp)), "Loadings_vec"=rep(1,Data$n_f*ncol(Data$Y_sp)-Data$n_f*(Data$n_f-1)/2), "log_kappa"=log(1), "Omega_xf"=matrix(0,nrow=mesh$n,ncol=Data$n_f) )
# Declare random
Random = c("Omega_xf")

# Initialization
Obj <- MakeADFun(data=Data, parameters=Params, random=Random, hessian=FALSE, inner.control=list(maxit=1000) )
table(names(Obj$env$last.par))

# Run model
Opt = TMBhelper::Optimize( obj=Obj, getsd=TRUE, newtonsteps=1 )

# Summarize
Report = Obj$report()
Report$Loadings_pf
Cov_pp = Report$Loadings_pf %*% t(Report$Loadings_pf)
