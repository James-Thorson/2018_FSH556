
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 6 -- 2D spatial models/Lecture/" )

library(TMB)
library(RandomFields)
library(Matrix)

###################
# Equal distance 2D autoregressive
###################

Dim = c("n_x"=10, "n_y"=10)
loc_xy = expand.grid("x"=1:Dim['n_x'], "y"=1:Dim['n_y'])
Scale = 5
Sigma2 = (0.5) ^ 2
beta0 = 3
prob_missing = 0.2

# Simulate spatial process
RMmodel = RMexp(var=Sigma2, scale=Scale)
epsilon_xy = array(RFsimulate(model=RMmodel, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1], dim=Dim)
image( z=epsilon_xy )

# SImulate counts
c_xy = array(NA, dim=dim(epsilon_xy))
for(x in 1:nrow(c_xy)){
for(y in 1:ncol(c_xy)){
  c_xy[x,y] = rpois(1, exp(beta0 + epsilon_xy[x,y]) )
  if( rbinom(n=1, size=1, prob=prob_missing)==1) c_xy[x,y] = NA
}}
true_abundance =  sum( exp(beta0 + epsilon_xy) )

# Generate sparse matrices for precision matrix of 2D AR1 process
M0 = as( ifelse(as.matrix(dist(loc_xy,diag=TRUE,upper=TRUE))==0,1,0), "dgTMatrix" )
M1 = as( ifelse(as.matrix(dist(loc_xy,diag=TRUE,upper=TRUE))==1,1,0), "dgTMatrix" )
M2 = as( ifelse(as.matrix(dist(loc_xy,diag=TRUE,upper=TRUE))==sqrt(2),1,0), "dgTMatrix" )

# Compile
Params = list( "beta0"=0, "ln_sigma2"=0, "logit_rho"=0, "epsilon_xy"=array(rnorm(prod(dim(loc_xy)))-100,dim=dim(epsilon_xy)) )
compile( "autoregressive_grid_V1.cpp" )
dyn.load( dynlib("autoregressive_grid_V1") )

######## Version 0 -- Sweep downstream
# Build object
Data = list("Options_vec"=c(0), "c_xy"=c_xy, "M0"=M0, "M1"=M1, "M2"=M2 )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid_V1" )
# Optimize
Opt0 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
h0 = Obj$env$spHess(random=TRUE)
report0 = Obj$report()

######## Version 1 -- Analytic precision matrix
# Build object
Data = list("Options_vec"=c(1), "c_xy"=c_xy, "M0"=M0, "M1"=M1, "M2"=M2 )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid_V1" )
# Optimize
Opt1 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
h1 = Obj$env$spHess(random=TRUE)
report1 = Obj$report()

######## Version 3 -- Built-in function for AR process
# Build object
Data = list("Options_vec"=c(3), "c_xy"=c_xy, "M0"=M0, "M1"=M1, "M2"=M2 )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid_V1" )
# Optimize
Opt3 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
h3 = Obj$env$spHess(random=TRUE)
report3 = Obj$report()

######## Version 4 -- Assemble sparse precision matrix using external computation (probably the easiest to generalize)
# Build object
Data = list("Options_vec"=c(3), "c_xy"=c_xy, "M0"=M0, "M1"=M1, "M2"=M2 )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid_V1" )
# Optimize
Opt4 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
h4 = Obj$env$spHess(random=TRUE)
report4 = Obj$report()

# Compare results
cbind( par0, par1, par3, par4 )
c( Opt0$run_time, Opt1$run_time, Opt3$run_time, Opt4$run_time )

# Compare hessian sparseness
image( h0, main="Version 0" ); dev.new()
image( h1, main="Version 1" ); dev.new()
image( h3, main="Version 3" ); dev.new()
image( h4, main="Version 4" )
