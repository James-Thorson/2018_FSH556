
# Generate 10 by 10 grid
Dim = c("n_x"=10, "n_y"=10)
loc_xy = expand.grid("x"=1:Dim['n_x'], "y"=1:Dim['n_y'])

# Define parameters
Scale = 2
Sigma2 = (0.5) ^ 2
beta0 = 3
prob_missing = 0.2

# Simulate spatial process
RMmodel = RMgauss(var=Sigma2, scale=Scale)
epsilon_xy = array(RFsimulate(model=RMmodel, x=loc_xy[,'x'], y=loc_xy[,'y'])@data[,1], dim=Dim)

# SImulate counts
c_xy = array(NA, dim=dim(epsilon_xy))
for(x in 1:nrow(c_xy)){
for(y in 1:ncol(c_xy)){
  c_xy[x,y] = rpois(1, exp(beta0 + epsilon_xy[x,y]) )
  if( rbinom(n=1, size=1, prob=prob_missing)==1) c_xy[x,y] = NA
}}

# Format data frame
DF = data.frame( "c_i"=as.vector(c_xy), "x_i"=as.vector(row(c_xy)), "y_i"=as.vector(col(c_xy)) )
DF = rbind( cbind(DF, "t_i"=1), cbind( 'c_i'=DF$c_i*2, DF[,c('x_i','y_i')], "t_i"=2) )

# Run GAM with 3D smoother
Gam = gam( c_i ~ s( x_i, y_i, t_i ), data=DF )

# Run GAM with 2D smoother in each individual year
Gam = gam( c_i ~ s( x_i, y_i, by=factor(t_i) ), data=DF )
