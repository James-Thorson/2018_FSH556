
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 2 -- mixed-effects/Lecture 2" )
Use_REML = FALSE

devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

######################
# Simulate data
######################

# Simulate predictors
group_i = rep( 1:10, each=10)
z_g = rnorm( length(unique(group_i)), mean=0, sd=1)
beta0 = 0

# Simulate response
y_i = z_g[group_i] + beta0 + rnorm( length(group_i), mean=0, sd=1)

######################
# Run in R
######################

library(lme4)
Lme = lmer( y_i ~ 1|factor(group_i), REML=Use_REML)

######################
# Run in TMB
######################

library(TMB)

# Compile model
Version = "linear_mixed_model"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "n_groups"=length(unique(group_i)), "g_i"=group_i-1, "y_i"=y_i)
Parameters = list( "beta0"=-10, "log_SD0"=2, "log_SDZ"=2, "z_g"=rep(0,Data$n_groups) )
Random = c("z_g")
if( Use_REML==TRUE ) Random = union( Random, "beta0")

# Build object
dyn.load( dynlib("linear_mixed_model") )
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
# Shrinkage estimator
######################

# jim's attempt at replicating this from first principles
Mu = mean(y_i)
  Mu_s = tapply( y_i, INDEX=group_i, FUN=mean)
Sigma = sd( Mu_s )
  Sigma_s = sd( y_i - Mu_s[group_i] )
Weights_hat = c( 1/Sigma^2, length(y_i)/length(unique(group_i))/Sigma_s^2 )
  Weights_hat = Weights_hat / sum(Weights_hat)

# Predictions
Mu_s_hat = ( Mu*Weights_hat[1] + Mu_s*Weights_hat[2] )

######################
# Compare estimates
######################

# Global mean
c( fixef(Lme), ParHat$beta0, Mu )

# Random effects
cbind( "True"=z_g, "Lme4"=ranef(Lme)[['factor(group_i)']], "TMB"=ParHat$z_g, "Shrinkage_estimator"=Mu_s-Mu )

# Variances
summary(Lme)
unlist( Report[c("SDZ","SD0")] )

