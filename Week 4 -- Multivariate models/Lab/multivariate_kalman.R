
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 4 -- Multivariate models/Lab" )
set.seed(3)

######################
# Simulate data
######################

np = 6
nt = 120
years_without_data = 101:120 # years with no data

# Assemble variance
L = matrix(0, nrow=np, ncol=np)
L[lower.tri(L,diag=TRUE)] = rnorm( np*(np+1)/2, mean=0.5, sd=1 )

# Other parameters
alpha = 0.04
CovP = L %*% t(L)
VarM = mean(diag(CovP))
x0 = rnorm(np, mean=3, sd=sqrt(diag(CovP)))

# Calculate Cholesky (R uses an upper-triangle Cholesky, so I transpose it!)
Lp = t(chol(CovP))

# Simulate predictors
x_tp = y_tp = matrix(NA, nrow=nt, ncol=np)
x_tp[1,] = x0
for( t in 2:nt ){
  x_tp[t,] = x_tp[t-1,] + (Lp%*%rnorm(np, mean=alpha, sd=1))[,1]
}
for( t in 1:nt ){
for( p in 1:np ){
  y_tp[t,p] = x_tp[t,p] + rnorm(1, mean=0, sd=sqrt(VarM))
}}

# Exclude years in forecast period
y_tp[years_without_data,] = NA

png( file="simulation.png", width=8, height=5, res=200, units="in")
  par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  matplot( x=1:nt, y=y_tp, type="l", lwd=2, lty="solid", col="blue", ylab="Value", xlab="Year" )
  #matplot( x=1:nt, y=x_tp, type="l", lty="solid", col="black", add=TRUE )
dev.off()

#######################
# Fit in TMB
#######################

library(TMB)

# Compile model
Version = "multivariate_kalman"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "y_tp"=y_tp )
Parameters = list( "log_VarM"=log(1), "x0_p"=rep(0,np), "Lp"=matrix(0,nrow=np,ncol=np), "delta_tp"=matrix(rnorm(nt*np),ncol=np) )
Parameters$Lp[lower.tri(Parameters$Lp,diag=TRUE)] = rnorm( np*(np+1)/2, mean=0.5, sd=1 )
Random = c("delta_tp")

# Build Map
Map = NULL
Map[["Lp"]] = matrix( 1:(np*np), nrow=np, ncol=np )
Map[["Lp"]][upper.tri(Map[["Lp"]])] = NA
Map[["Lp"]] = factor(Map[["Lp"]])

# Build object
dyn.load( dynlib("multivariate_kalman") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, map=Map)  #

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
Opt = TMBhelper::Optimize( obj=Obj, getsd=TRUE, newtonsteps=1 )

# Get reporting and SEs
Report = Obj$report()
X_tp = array( summary(Opt$SD,"report")[rownames(summary(Opt$SD,"report"))=="x_tp",], dim=c(dim(Report$x_tp),2), dimnames=list(NULL,NULL,c("Estimate","Std. Error")) )

# Plot fit
for( z in 1:3 ){
  png( file=paste0("fit_",z,".png"), width=8, height=5, res=200, units="in")
    par( mfrow=c(2,3), mar=c(1,1,2,1), mgp=c(2,0.5,0), tck=-0.02, oma=c(2,2,0,0) )
    for( p in 1:ncol(X_tp) ){
      matplot( x=1:nt, y=y_tp[,p], col="blue", cex=1.2, type="l", lty="solid", xlab="", ylab="", pch=20, main=paste0("Population ",p), ylim=range(x_tp) )
      if(z>=2){
        matplot( x=1:nt, y=X_tp[,p,'Estimate'], col="red", type="p", pch=21, add=TRUE )
        for( t in 1:nrow(X_tp) ){
          lines( x=rep(t,2), y=X_tp[t,p,'Estimate']+c(-1.96,1.96)*X_tp[t,p,'Std. Error'], col="red" )
        }
      }
      if(z>=3) matplot( x=1:nt, y=x_tp[,p], col="black", lwd=2, add=TRUE, lty="solid", type="l" )
    }
    mtext( side=1:2, outer=TRUE, text=c("Year","Value"), line=c(1,1) )
  dev.off()
}

# Compare covariance
plot( x=CovP, y=Report$CovP_pp )
write.csv( Opt$diagnostics, file="parameter_estimates.csv" )

###########################
# Conduct simulation experiment
###########################

# Compile model
Version = "multivariate_kalman_with_simulation"
compile( paste0(Version,".cpp") )

# Build object
dyn.load( dynlib("multivariate_kalman_with_simulation") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, map=Map, DLL="multivariate_kalman_with_simulation")  #

# Optimize
Opt = TMBhelper::Optimize( obj=Obj, getsd=TRUE, newtonsteps=1 )

# Don't run -- very slow !
if(FALSE){
  # Do simulation
  SimPar = matrix( NA, nrow=100, ncol=length(Opt$par), dimnames=list(NULL,names(Opt$par)) )
  for( i in 1:nrow(SimPar) ){
    SimData <- Obj$simulate( par=Obj$env$last.par, complete=TRUE )
    Obj2 <- MakeADFun(data=SimData['y_tp'], parameters=Parameters, random=Random, map=Map, silent=TRUE, DLL="multivariate_kalman_with_simulation")
    SimPar[i,] = TMBhelper::Optimize(obj=Obj2, getsd=FALSE, newtonsteps=1)$par
    if( i%%10 == 0 ) message(paste0("Finished replicate ",i))
  }

  # Plot distribution against true value
  png( file="parametric_bootstrap.png", width=8, height=5, res=200, units="in")
    par( mfrow=c(4,7), mar=c(1,1,2,1), mgp=c(2,0.5,0), tck=-0.02, oma=c(1,1,0,0) )
    for( z in 1:length(Obj$par)){
      hist( SimPar[,z], breaks=10, xlab="", ylab="", main=names(Obj$par)[z], xlim=range(c(SimPar[,z],Opt$par[z]),na.rm=TRUE), col=rgb(1,0,0,0.2) )
      abline( v=Opt$par[z], lwd=2 )
    }
  dev.off()
}
