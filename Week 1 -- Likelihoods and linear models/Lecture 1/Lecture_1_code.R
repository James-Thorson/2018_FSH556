
# Install packages
install.packages("devtools")
devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")

# Load TMB
library( TMB )

# Set working directory (change on other machine!)
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 1 -- Likelihoods and linear models/Lecture 1" )


############
# Example 1 -- average CPUE for canary rockfish
############

# Data
data(WCGBTS_Canary_example, package="SpatialDeltaGLMM")
CPUE = WCGBTS_Canary_example$HAUL_WT_KG
# Visualize
png( file="Canary_histogram.png", width=4, height=4, res=200, units="in")
  par( mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
  hist( log(1+CPUE) )
dev.off()
# Probability for example parameters
dnorm( CPUE, mean=20, sd=2 )           # Likelihood for each datum
sum(dnorm( CPUE, mean=20, sd=2, log=TRUE ))     # Log-likelihood for all data

######### Method 1 -- Pre-made functions in R
Lm = lm( CPUE ~ 1 )
summary(Lm)

######### Method 2 -- Optimize using R
# Step 1 -- define function
NegLogLike_Fn = function(Par, Data){
  # Parameters
  Mean_hat = Par[1]
  SD_hat = Par[2]
  # Log-likelihood
  LogLike_i = dnorm( Data$y, mean=Mean_hat, sd=SD_hat, log=TRUE )
  NegLogLike = -1 * sum(LogLike_i)
  return( NegLogLike )
}
# step 2 -- minimize negative log-likelihood to estimate parameters
Data = list( 'y'=CPUE )
Start = c(1,1)
NegLogLike_Fn(Par=Start, Data=Data)
Opt = optim( par=Start, fn=NegLogLike_Fn, Data=Data, lower=c(0.01,0.01), upper=Inf, method="L-BFGS-B", hessian=TRUE )
print( Opt$par ) # Estimated parameters
print(sqrt( diag( solve(Opt$hessian) )) ) # standard errors

###### Method 3 -- Optimize using TMB
# Step 1 -- make and compile template file
compile( "linear_model_v1.cpp" )

# Step 2 -- build inputs and object
dyn.load( dynlib("linear_model_v1") )
Params = list("mean"=0, "log_sd"=0)
Data = list( "y_i"=CPUE )
Obj = MakeADFun( data=Data, parameters=Params, DLL="linear_model_v1")

# Step 3 -- test and optimize
Obj$fn( Obj$par )
Obj$gr( Obj$par )
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Opt$diagnostics = data.frame( "name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
Opt$par # estimated parameters
SD = sdreport( Obj ) # standard errors

############
# Example 2 -- spatially stratified CPUE for canary rockfish
############

# Define a new design matrix
X = cbind( "CA"=ifelse(WCGBTS_Canary_example$BEST_LAT_DD<42,1,0), "OR"=ifelse(WCGBTS_Canary_example$BEST_LAT_DD>=42&WCGBTS_Canary_example$BEST_LAT_DD<46,1,0), "WA"=ifelse(WCGBTS_Canary_example$BEST_LAT_DD>=46,1,0), "Pass"=WCGBTS_Canary_example$PASS-1.5 )

# Step 1 -- make and compile template file
compile( "linear_model_v2.cpp" )

# Step 2 -- build inputs and object
dyn.load( dynlib("linear_model_v2") )
Params = list("b_j"=rep(0,ncol(X)), "log_sd"=0)
Data = list( "y_i"=CPUE, "X_ij"=X )
Obj = MakeADFun( data=Data, parameters=Params, DLL="linear_model_v2")

# Step 3 -- test and optimize
Obj$fn( Obj$par )
Obj$gr( Obj$par )
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Opt$diagnostics = data.frame( "name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
Opt$par # estimated parameters
SD = sdreport( Obj ) # standard errors




