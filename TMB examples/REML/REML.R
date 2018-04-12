
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/TMB examples/REML" )
#devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")

library(lme4)
library(TMB)

# Compile model
Version = "REML"
compile( paste0(Version,".cpp") )
dyn.load( dynlib("REML") )

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
# Run WITHOUT REML in R and TMB
######################

# In R
Lme_noREML = lmer( y_i ~ 1|factor(group_i), REML=FALSE)

# In TMB
Data = list( "n_groups"=length(unique(group_i)), "g_i"=group_i-1, "y_i"=y_i)
Parameters = list( "beta0"=-10, "log_SD0"=2, "log_SDZ"=2, "z_g"=rep(0,Data$n_groups) )
Random = c("z_g")
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  #
Opt_noREML = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )

# TMB
exp(Opt_noREML$par['log_SDZ'])
# R
summary(Lme_noREML)

######################
# Run WITH REML in R and TMB
######################

# In R
Lme_REML = lmer( y_i ~ 1|factor(group_i), REML=TRUE)

# In TMB
Data = list( "n_groups"=length(unique(group_i)), "g_i"=group_i-1, "y_i"=y_i)
Parameters = list( "beta0"=-10, "log_SD0"=2, "log_SDZ"=2, "z_g"=rep(0,Data$n_groups) )
Random = c("z_g","beta0")
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  #
Opt_REML = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )

# TMB
exp(Opt_REML$par['log_SDZ'])
# R
summary(Lme_REML)

