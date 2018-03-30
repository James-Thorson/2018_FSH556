
# Parameters
Mean = 3
CV = 1

# Simulate data
Shape = 1/CV^2
Scale = Mean*CV^2
y_i = rgamma( n=1000, shape=Shape, scale=Scale )

# Check implementation
mean( y_i ) # Should be close to "Mean"
sd(y_i)/mean(y_i) # Should be close to "CV"

# Run in TMB
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/TMB examples/gamma distribution/" )
library(TMB)

compile("dgamma_example.cpp")
dyn.load( dynlib("dgamma_example") )

Params = list( "log_mean"=log(1), "log_CV"=5 )
Data = list( "y_i"=y_i )
Obj = MakeADFun( parameters=Params, data=Data )

Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )

# Compare true and estimated values
cbind( "True"=c(Mean,CV), summary(Opt$SD,"report") )
