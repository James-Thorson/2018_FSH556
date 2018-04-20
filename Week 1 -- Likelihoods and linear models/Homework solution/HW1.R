
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 1 -- Likelihoods and linear models/Homework solution" )

library(TMB)
library( TMBdebug )  # Used to debug TMB models in windows (doesn't hurt for other operating systems)
library( SpatialDeltaGLMM )
library( statmod ) # Inverse-gaussian distribution pinvgauss

# compile template file
compile( "HW1.cpp" )
dyn.load( dynlib("HW1") )

###########
# Part 1 -- analysis of real data for Alaska pollock
###########

# Settings
set.seed(1)
K = 10

# Load data and divide into partitions
data( EBS_pollock_data )
CPUE = EBS_pollock_data$catch
X = cbind( "Intercept"=rep(1,length(CPUE)) )
Partition_i = sample( x=1:K, size=length(CPUE), replace=TRUE )

# Stuff to record
PredNLL_kj = matrix(NA, ncol=3, nrow=K)
Table1 = array(NA, dim=c(3,3), dimnames=list(c("lognormal","gamma","invgauss"),c("negloglike","num_params","pred_negloglike_per_datum")) )

# Cross-validation
for(j in 1:3){
for(k in 1:K){
  Params = list("b_j"=rep(0,ncol(X)), "theta_z"=c(0,0))
  Data = list( "Options_vec"=j-1, "y_i"=CPUE, "X_ij"=X, predTF_i=ifelse(Partition_i==k,1,0) )
  Obj = MakeADFun( data=Data, parameters=Params, DLL="HW1")
  Obj$env$beSilent()

  # Step 3 -- test and optimize
  Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )

  # Extract stuff
  Report = Obj$report()
  PredNLL_kj[k,j] = Report$pred_jnll / sum(Partition_i==k)
  if( PredNLL_kj[k,j]>10 ) stop("Check results")
}}

# Build table of results
for(j in 1:3){
  Params = list("b_j"=rep(0,ncol(X)), "theta_z"=c(0,0))
  Data = list( "Options_vec"=j-1, "y_i"=CPUE, "X_ij"=X, predTF_i=rep(0,length(CPUE)) )
  Obj = MakeADFun( data=Data, parameters=Params, DLL="HW1")

  # Step 3 -- test and optimize
  for(i in 1:3) Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
  Table1[j,] = c(Opt$objective, length(Opt$par), mean(PredNLL_kj[,j]))
}

# Look at results
print(Table1)

###########
# Part 2 -- Simulation design
###########

# Settings
set.seed(1)
n_rep = 100

# Stuff to record
Results = array(NA, dim=c(n_rep,3,3,2), dimnames=list(NULL,paste("True=",c("lognormal","gamma","invgauss"),sep=""),paste("Est=",c("lognormal","gamma","invgauss"),sep=""),c("True","Est")) )

# Loop through simulation scenarios
for(simI in 1:3){

  # Get parameter values from fitting to real data
  Params = list("b_j"=rep(0,ncol(X)), "theta_z"=c(0,0))
  Data = list( "Options_vec"=simI-1, "y_i"=CPUE, "X_ij"=X, predTF_i=rep(0,length(CPUE)) )
  Obj = MakeADFun( data=Data, parameters=Params, DLL="HW1")
  Opt_real = nlminb( start=Obj$env$last.par.best, objective=Obj$fn, gradient=Obj$gr )

  # Loop through replicates and estimation models
  for(repI in 1:n_rep){

    # Simulate data
    if(simI==1) ysim_i = rbinom(n=length(CPUE), size=1, prob=1-plogis(Opt_real$par)[2]) * rlnorm(n=length(CPUE), meanlog=Opt_real$par[1], sdlog=exp(Opt_real$par[3]))
    if(simI==2) ysim_i = rbinom(n=length(CPUE), size=1, prob=1-plogis(Opt_real$par)[2]) * rgamma(n=length(CPUE), shape=1/(exp(Opt_real$par[3])^2), scale=exp(Opt_real$par[1])*(exp(Opt_real$par[3])^2))
    if(simI==3) ysim_i = rbinom(n=length(CPUE), size=1, prob=1-plogis(Opt_real$par)[2]) * rinvgauss(n=length(CPUE), mean=exp(Opt_real$par[1]), shape=exp(Opt_real$par[3]))

    # Loop through estimation models
    for(estI in 1:3){
      # Get parameter values
      Params = list("b_j"=rep(0,ncol(X)), "theta_z"=c(0,0))
      Data = list( "Options_vec"=estI-1, "y_i"=ysim_i, "X_ij"=X, predTF_i=rep(0,length(CPUE)) )
      Obj = MakeADFun( data=Data, parameters=Params, DLL="HW1")
      Opt = nlminb( start=Obj$env$last.par.best, objective=Obj$fn, gradient=Obj$gr )

      # Record estimated and true intercept
      Results[repI,simI,estI,1] = Opt_real$par[1]
      Results[repI,simI,estI,2] = Opt$par[1]
    }
  }
}

# Plot
par( mfrow=c(3,3), mar=c(3,3,0,0), mgp=c(2,0.5,0), tck=-0.02, yaxs="i")
for( simI in 1:3 ){
for( estI in 1:3 ){
  hist( Results[,simI,estI,2], xlim=range(Results,na.rm=TRUE) )
  abline( v=Results[1,simI,estI,1], lwd=3, lty="dotted" )
}}



