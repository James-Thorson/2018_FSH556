
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 4 -- Multivariate models/Lecture" )
set.seed(2)

######################
# Check Cholesky method
######################

Rho = 0.2
Var = 0.7
np = 6
ni = 1e3

# Assemble variance
Cov = (Var-Rho)*diag(np) + Rho*matrix(1,nrow=np)%*%matrix(1,ncol=np)

# Calculate Cholesky (R uses an upper-triangle Cholesky, so I transpose it!)
L = t(chol(Cov))
L %*% t(L) # Check math

# Generate many draws
Delta_ip = matrix( rnorm(ni*np), ncol=np)

# Transform
Epsilon_ip = L %*% t(Delta_ip)
cov(t(Epsilon_ip))

######################
# Simulate data
######################

nt = 100
x0 = rnorm(np, mean=3, sd=sqrt(diag(Cov)))
alpha = 0.04

# Calculate Cholesky (R uses an upper-triangle Cholesky, so I transpose it!)
Lp = t(chol(Cov))

# Simulate predictors
x_tp = y_tp = matrix(NA, nrow=nt, ncol=np)
x_tp[1,] = x0
for( t in 2:nt ){
  #x_tp[t,] = x_tp[t-1,] + rmvnorm( 1, mean=rep(alpha,np), sigma=CovP )
  x_tp[t,] = x_tp[t-1,] + (Lp%*%rnorm(np, mean=alpha, sd=1))[,1]
}

png( file="simulation.png", width=8, height=5, res=200, units="in")
  par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
  matplot( x=1:nt, y=x_tp, type="l", lty="solid", col="black", ylab="Value", xlab="Year" )
dev.off()

#######################
#
#######################


