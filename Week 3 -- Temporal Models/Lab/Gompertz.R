
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 3 -- Temporal Models/Lab" )
Use_REML = FALSE
set.seed(2)

#####################
# Explore Gompertz model
#####################

beta= 0.2
alpha = 1
d_equil = exp(alpha/beta)

d_1 = seq(0,d_equil*2,length=1e4)
d_2 = d_1 * exp(alpha) * exp( - beta*log(d_1))

# Dynamics
png( file="Gompertz_dynamics.png", width=8, height=4, res=200, units="in")
  par( mfrow=c(1,2), mar=c(3,3.5,2,0), mgp=c(1.75,0.25,0), tck=-0.02)
  plot( x=d_1, y=d_2, type="l", lwd=3, xlab=expression(Biomass[t]), ylab=expression(Biomass[t+1]), main="Production")
  abline( a=0, b=1, lty="dotted")
  # Log-dynamics
  plot( x=log(d_1[-1]), y=log(d_2[-1]/d_1[-1]), type="l", lwd=3, xlab=expression(log(Biomass[t])), ylab=expression(log(Biomass[t+1]/Biomass[t])), main="log-Biomass ratio" )
  abline( a=1, b=0, lty="dotted")
dev.off()

######################
# Simulate data
######################

nt = 100
log_d0 = 3
sigmaP = 0.5
sigmaM = 0.5
alpha = 1
beta = 0.1

# Simulate predictors
log_d_t = log_b_t = rep(NA, nt)
log_d_t[1] = log_d0
for( t in 2:nt ){
  log_d_t[t] = alpha + (1-beta)*log_d_t[t-1] + rnorm( 1, mean=0, sd=sigmaP )
}
for( t in 1:nt ){
  log_b_t[t] = log_d_t[t] + rnorm( 1, mean=0, sd=sigmaM )
}



######################
# Run in TMB
######################

library(TMB)

# Compile model
Version = "gompertz"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "nt"=nt, "log_b_t"=log_b_t )
Parameters = list( "log_d0"=0, "log_sigmaP"=1, "log_sigmaM"=1, "alpha"=0, "rho"=0, "log_d_t"=rep(0,nt) )
Random = c("log_d_t")
if( Use_REML==TRUE ) Random = union( Random, c("log_d0","alpha","rho") )

# Build object
dyn.load( dynlib("gompertz") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  #

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
start_time = Sys.time()
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1) )
  Opt[["final_gradient"]] = Obj$gr( Opt$par )
  Opt[["total_time"]] = Sys.time() - start_time

# Get reporting and SEs
Report = Obj$report()
  SD = sdreport( Obj )

for( z in 1:3 ){
  png( file=paste0("gompertz_",z,".png"), width=8, height=5, res=200, units="in" )
    par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    plot( x=1:nt, y=log_b_t, col="blue", cex=1.2, xlab="Time", ylab="Value", pch=20 )
    if(z>=2){
      points( x=1:nt, y=as.list(SD,"Estimate")$log_d_t, col="red" )
      for( t in 1:nt) lines( x=rep(t,2), y=as.list(SD,"Estimate")$log_d_t[t]+c(-1.96,1.96)*as.list(SD,"Std. Error")$log_d_t[t], col="red" )
    }
    if(z>=3) lines( x=1:nt, y=log_d_t, col="black", lwd=2 )
  dev.off()
}

