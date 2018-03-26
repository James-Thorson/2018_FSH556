
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2016_Spatio-temporal_models/TMB examples/map argument" )

######################
# Simulate data
######################

# Simulate predictors
Factor = rep( 1:10, each=10)
Z = rnorm( length(unique(Factor)), mean=0, sd=1)
X0 = 0

# Simulate response
Y = Z[Factor] + X0 + rnorm( length(Factor), mean=0, sd=1)

######################
# Run in TMB
######################

library(TMB)

# Compile model
Version = "linear_mixed_model"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "n_data"=length(Y), "n_factors"=length(unique(Factor)), "Factor"=Factor-1, "Y"=Y)
Parameters = list( "X0"=-10, "log_SD0"=2, "log_SDZ"=2, "Z"=rep(0,Data$n_factor) )

# Turn off random effects
Map = list()
Map[["log_SDZ"]] = factor(NA)
Map[["Z"]] = factor( rep(NA,Data$n_factor) )

# Build object
dyn.load( dynlib("linear_mixed_model") )
Obj = MakeADFun(data=Data, parameters=Parameters, map=Map)  #

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
start_time = Sys.time()
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1) )

# Get reporting and SEs
Report = Obj$report()
  SD = sdreport( Obj )

