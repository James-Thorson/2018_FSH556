

library(TMB)
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 5 -- 1D spatial models/Lecture/" )

############################
# Visualize Gaussian process
############################

library(RandomFields)

png( "Gaussian_process.png", width=8, height=5, res=200, units="in")
  par( mfrow=c(2,2), mar=c(2,2,2,0), mgp=c(2,0.5,0) )
  for(i in 1:4){
    Scale = c(10, 10, 25, 25)[i]
    SD = c(0.5, 1, 0.5, 1)[i]
    x_i = seq(0,100,length=100)
    RMmodel = RMgauss(var=SD^2, scale=Scale)
    plot( 1, type="n", xlim=range(x_i), ylim=c(-3,3), main=paste0("SD=",SD," Scale=",Scale), xlab="", ylab="")
    for(j in 1:20){
      gp_i = RFsimulate(model=RMmodel, x=x_i, y=rep(0,length(x_i)))@data[,1]
      lines( x=x_i, y=gp_i )
    }
  }
dev.off()

###################
# Visualize autocorrelation
###################

Rho = 0.8
X = 0:10
Corr = Rho^X

png( file="Autocorrelation.png", width=4, height=4, res=200, units='in')
  plot( x=X, y=Corr, xlab="distance", ylab="Correlation", ylim=c(0,1) )
dev.off()

###################
# Visualize autocorrelation vs. random-walk
###################

Rho = c(-0.5,0,0.5,0.99)
Marginal_Sigma = 1
Conditional_Sigma = sqrt( 1-Rho^2 ) * Marginal_Sigma

Y = X = 1:1e2
png( file="Autocorrelated_series.png", width=4, height=6, res=200, units='in')
  par( mfrow=c(4,1), mgp=c(2,0.5,0), mar=c(2,2,1,0), xaxs="i")
  for(i in 1:4){
    Y[1] = 0
    for(s in 2:length(X)) Y[s] = Y[s-1]*Rho[i] + rnorm(1, mean=0, sd=Conditional_Sigma[i])
    plot( x=X, y=Y, xlab="location", ylab="value", ylim=c(-3,3), type="l", main=paste0("Rho = ",Rho[i]) )
  }
dev.off()

###################
# Math check on inverse-covariance matrix
###################
Rho = 0.5
Sigma2 = 2 ^ 2
x = 1:10
Dist = outer(x, x, FUN=function(a,b){abs(a-b)})

# Calculate Q directly
Cov = Sigma2 / (1-Rho^2) * Rho^Dist
Q = solve(Cov)

# Calculate Q analytically
M = diag(length(x)) * (1+Rho^2)
M[ cbind(1:9,2:10) ] = -Rho
M[ cbind(2:10,1:9) ] = -Rho
Q2 = 1/Sigma2 * M

###################
# Equal distance autoregressive
###################

x = 1:100
Rho = 0.8
Sigma2 = (0.5) ^ 2
n_rep = 3
beta0 = 3

# Simulate spatial process
epsilon_s = rep(NA, length(x))
epsilon_s[1] = rnorm(1, mean=0, sd=sqrt(Sigma2))
for(s in 2:length(x)) epsilon_s[s] = Rho*epsilon_s[s-1] + rnorm(1, mean=0, sd=sqrt(Sigma2))

# SImulate counts
c_si = matrix( nrow=length(x), ncol=n_rep)
for(s in 1:nrow(c_si)){
for(i in 1:ncol(c_si)){
  c_si[s,i] = rpois(1, exp(beta0 + epsilon_s[s]) )
}}

# Compile
Params = list( "beta0"=0, "ln_sigma2"=0, "logit_rho"=0, "epsilon_s"=rnorm(length(x)) )
compile( "autoregressive_V1.cpp" )
dyn.load( dynlib("autoregressive_V1") )

######## Version 0 -- Stochastic process with automatic sparseness detection
# Build object
Data = list("Options_vec"=c(0), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL="autoregressive_V1" )
# Optimize
Opt0 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
par0 = Opt0$par
h0 = Obj$env$spHess( random=TRUE )

######## Version 1 -- Analytic precision matrix
# Build object
Data = list("Options_vec"=c(1), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL="autoregressive_V1" )
# Optimize
Opt1 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
par1 = Opt1$par
h1 = Obj$env$spHess(random=TRUE)

######## Version 2 -- Covariance and built-in function
# Build object
Data = list("Options_vec"=c(2), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL="autoregressive_V1" )
# Optimize
Opt2 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
par2 = Opt2$par
h2 = Obj$env$spHess(random=TRUE)

######## Version 3 -- Built-in function for AR process
# Build object
Data = list("Options_vec"=c(3), "c_si"=c_si )
Obj = MakeADFun( data=Data, parameters=Params, random="epsilon_s", DLL="autoregressive_V1" )
# Optimize
Opt3 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
par3 = Opt3$par
h3 = Obj$env$spHess(random=TRUE)

# Compare timing
Opt0$run_time
Opt1$run_time
Opt2$run_time
Opt3$run_time

# Compare separability
library(INLA)
image(h0, main="Version 0"); dev.new()
image(h1, main="Version 1"); dev.new()
image(h2, main="Version 2"); dev.new()
image(h3, main="Version 3")

