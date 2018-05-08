
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2016_Spatio-temporal_models/Week 2 -- mixed-effects/Homework solution" )

# Simulation settings
Nrep = 100
Target_Coverage = 0.5 # Target nominal CI coverage

# Parameter settings
Nsite = 10
Nobs_per_site = 10
Site_logMean = 2
Site_logSd = 1
Overdispersion_logSd = 0.5

# Libraries
library(TMB)
#library(TMBdebug)

# Compile model
Version = "HW2"
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )

# Recording object
Results = array(NA, dim=c(Nrep,4,3), dimnames=list(NULL,paste0("Random:",c("None","Site","Obs","Both")),c("mu_hat","SE_mu_hat","CoveredTF")))

# Helper functions
CoveredTF = function( x, CI ) x>CI[1] & x<CI[2]

# Simulation loop
for(repI in 1:Nrep){

  # Simulation
  s_i = rep( 1:Nsite, each=Nobs_per_site)
  epsilon_s = rnorm( Nsite, mean=0, sd=Site_logSd )
  delta_i = rnorm(Nobs_per_site*Nsite, mean=0, sd=Overdispersion_logSd)
  yhat_i = exp( Site_logMean + epsilon_s[s_i] + delta_i )
  y_i = rpois( Nsite*Nobs_per_site, lambda=yhat_i )

  # Loop through estimation models
  for(estI in 1:4){

    # Build inputs
    Data = list( "n_i"=length(y_i), "n_s"=length(unique(s_i)), "s_i"=s_i-1, "y_i"=y_i)
    Parameters = list( "beta0"=-10, "log_sd_epsilon"=2, "log_sd_delta"=2, "epsilon_s"=rep(0,Data$n_s), "delta_i"=rep(0,Data$n_i) )
    Random = c("epsilon_s", "delta_i")

    # Turn off parameters for some estimation models
    Map = list()
    # Turn off overdispersion
    if(estI %in% c(1,2)){
      Map[["delta_i"]] = factor(rep(NA,Data$n_i))
      Map[["log_sd_delta"]] = factor(NA)
    }
    # Turn off site-level variation
    if(estI %in% c(1,3)){
      Map[["epsilon_s"]] = factor(rep(NA,Data$n_s))
      Map[["log_sd_epsilon"]] = factor(NA)
    }
    Random = setdiff( Random, names(Map))
    if( length(Random)==0 ) Random = NULL

    # Build object
    Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, map=Map)
    Obj$env$beSilent()

    # Optimize
    Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
    SD = sdreport( Obj )

    # Save results
    Results[repI,estI,c("mu_hat","SE_mu_hat")] = summary(SD)['beta0',]
    Results[repI,estI,'CoveredTF'] = CoveredTF(CI=Results[repI,estI,'mu_hat']+c(-1,1)*qnorm(1-Target_Coverage/2)*Results[repI,estI,'SE_mu_hat'], x=Site_logMean)
  }
}

#################
# Analysis
#################

# Plot estimates
par( mfrow=c(1,4), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02, yaxs="i")
for(estI in 1:4){
  hist( Results[,estI,'mu_hat'], breaks=seq(1,5,length=25), main=dimnames(Results)[[2]][estI], xlab="", ylab="")
  abline( v=Site_logMean, lwd=3, lty="dotted" )
}

# Summarize coverage results
apply( Results[,,'CoveredTF'], MARGIN=2, FUN=mean, na.rm=TRUE )

# Summarize bias
apply( Results[,,'mu_hat'], MARGIN=2, FUN=mean, na.rm=TRUE )



