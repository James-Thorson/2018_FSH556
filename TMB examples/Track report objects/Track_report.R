
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/TMB examples/Track report objects/" )

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
compile( "autoregressive_grid.cpp" )
dyn.load( dynlib("autoregressive_grid") )

# Build object
Data = list("c_xy"=c_xy, "M0"=M0, "M1"=M1, "M2"=M2 )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_xy", DLL="autoregressive_grid" )

# Overload marginal log-likelihood function
Obj$fn_orig = Obj$fn
Obj$fn = function( par ){
  MargNLL = Obj$fn_orig( par )
  Report = Obj$report()
  Total_Abundance = rbind(Total_Abundance, c(nrow(Total_Abundance)+1, Report$Total_Abundance) )
  colnames(Total_Abundance) = c("Function_call","Total_Abundance")
  assign( "Total_Abundance", Total_Abundance, pos=.GlobalEnv )
  return( MargNLL )
}

# Optimize
Total_Abundance = data.frame()
Opt = TMBhelper::Optimize( obj=Obj, fn=Obj$fn_trace, newtonsteps=1 )

