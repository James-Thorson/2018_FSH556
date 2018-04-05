
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 2 -- mixed-effects/Lab 2" )
Use_REML = TRUE

#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

############
# Generalized linear mixed model
############
library(lme4)

###### Simulate data
# Parameters
Nsite = 10
Nobs_per_site = 10
Site_logMean = log(10)
Site_logSd = 1

# Bookkeeping
s_i = rep( 1:Nsite, each=Nobs_per_site)

# Simulation
z_s = rnorm( Nsite, mean=0, sd=Site_logSd )
Mean_s = exp( Site_logMean + z_s )
y_i = rpois( Nsite*Nobs_per_site, lambda=Mean_s[s_i] )

# Plot data
library(lattice)
histogram( ~ y_i | factor(s_i), breaks=seq( min(y_i), max(y_i), length=10), type="density", panel=function(x,...){ panel.histogram(x, ...); panel.mathdensity(dmath=dnorm, col="black", args = list(mean=mean(x),sd=sd(x))) } )      #

###### Fit using R
# No site level (Not recommended)
GLM = glm( y_i ~ 1, family="poisson" )
print( summary(GLM) )

# Using fixed effects (Not recommended)
GLM = glm( y_i ~ 0 + factor(s_i), family="poisson" )
print( summary(GLM) )

# Using mixed effects (Recommended) -- doesn't appear to use REML
library(lme4)
GLMM = glmer( y_i ~ 1 + (1 | factor(s_i)), family="poisson" )
print( summary(GLMM) )

####################
# Fit using TMB
####################
library(TMB)

# Compile model
Version = "glmm"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "n_y"=length(y_i), "n_s"=length(unique(s_i)), "s_i"=s_i-1, "y_i"=y_i)
Parameters = list( "beta0"=-10, "log_sdz"=2, "z_s"=rep(0,Data$n_s) )
Random = c("z_s")
if( Use_REML==TRUE ) Random = union( Random, "beta0")

# Build object
dyn.load( dynlib(Version) )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  #

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )

# Get reporting and SEs
Report = Obj$report()
ParHat = as.list( Opt$SD, "Estimate" )

######################
# Compare estimates
######################

# Global mean
c( "Lme4"=fixef(GLMM), "TMB"=ParHat$beta0 )

# Random effects
cbind( "True"=z_s, "Lme4"=ranef(GLMM)[['factor(s_i)']], "TMB"=ParHat$z_s )


