


#########################
# Convert from mesh vertices to average triangle value, and back
#########################

devtools::install_github( "james-thorson/movement_tools" )

# File structure
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 9 -- Movement and network models/Lab 9/" )

library( MovementTools )
library( INLA )
library( RandomFields )
library( TMB )
library( expm ) # %^%
library( RANN ) # nn2()

source( "Fn_spatial_Gompertz_with_movement.R" )

# Domain
set.seed(1)
n_s = 50      # Number of samples
Cutoff = 1e-1   # 1e-1
Resolution_vec = c( "x"=50, "y"=100)

# Simulation model
SimulationDomain = c("Mesh", "Grid")[1]
Dynamical_Model = c("Gompertz", "Ricker")[1]
n_t = 20
logmeanu0 = 1
alpha = 1
beta = 0.5
mvec = rep(1, 4)   # Movement in each of four cardinal directions
SD_omega = 0.4       # SD of spatial variation
SD_epsilon = 0.6     # SD of spatio-temporal variation
Scale = 0.2          # Range/2, where Range is distance with 10% correlation

# Estimation model
n_tdiv = 3         # Number of iterations in Euler approximation
Include_movement = 1 # 0=No, 1=Stationary diffusion, 2=Independent-directions

# Derived
if( Dynamical_Model=="Gompertz" ){
  rho = 1 - beta
  phi = logmeanu0 - alpha/(1-rho)
}

# Locations on a grid
loc_g = expand.grid( "x"=seq(0,1,length=Resolution_vec["x"]), "y"=seq(0,1,length=Resolution_vec["y"]) )

# Locations for sampling
loc_s = cbind( "x"=runif(n_s), "y"=runif(n_s) )

# Make necessary mesh
Boundary = NULL   # inla.nonconvex.hull(points=cbind(x_stations, y_stations), convex=-0.1)
MeshList = Make_Movement_Mesh( loc_orig=loc_s, boundary=Boundary, Cutoff=Cutoff )

# Plot
par( mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
plot( loc_s, col="red", xlim=range(MeshList$loc_v), ylim=range(MeshList$loc_v), type="p", pch=20, cex=2)
plot( MeshList$mesh_domain, add=TRUE )

# Triangle info
TriList = TriList_Fn( MeshList$mesh_domain )
# Assign samples to triangles
r_s = Loc2Tri_Fn( locmat=loc_s, TriList=TriList )

# Movement operator on a mesh
MoveList = MovementMatrix_Fn(mesh=MeshList$mesh_domain, TriList=TriList)

# Assign samples to grid cells
NN_grid = nn2( data=loc_g, query=loc_s, k=1 )
g_s = NN_grid$nn.idx[,1]

######################
# Simulate data
######################

if( SimulationDomain=="Mesh" ){
  # Movement matrix
  MoveMat_mesh = expm( mvec[1]*MoveList[["M1"]] + mvec[2]*MoveList[["M2"]] + mvec[3]*MoveList[["M3"]] + mvec[4]*MoveList[["M4"]] )
  image( log(ifelse(MoveMat_mesh<0.0001,NA,MoveMat_mesh)) )
  # Simulate
  SimList = Sim_Fn( MoveMat=MoveMat_mesh, SD_omega=SD_omega, SD_epsilon=SD_epsilon, Scale=Scale, Dynamical_Model=Dynamical_Model, n_s=n_s, n_t=n_t, r_s=r_s, n_r=MeshList$n_r, loc_r=MeshList$loc_r, alpha=alpha, beta=beta)
}

if( SimulationDomain=="Grid" ){
  # Compute operators on a grid
  GridOperator = GridOperator_Fn( loc_g )
  # Movement matrix
  mparams = c("diffuse_x"=mean(mvec[c(1,3)]^2), "diffuse_y"=mean(mvec[c(2,4)]^2), "advect_east"=(mvec[1]-mvec[3])^2, "advect_south"=(mvec[4]-mvec[2])^2 ) / 10
  MoveMat_grid = ( mparams["diffuse_x"]*GridOperator[["nablaxx"]] + mparams["diffuse_y"]*GridOperator[["nablayy"]] + mparams["advect_east"]*GridOperator[["gradx"]] + mparams["advect_south"]*GridOperator[["grady"]] )
  MoveMat_grid = Matrix::expm( MoveMat_grid )
  image( MoveMat_grid, zlim=c(0,1) )
  # Simulate data
  SimList = Sim_Fn( MoveMat=MoveMat_grid, SD_omega=SD_omega, SD_epsilon=SD_epsilon, Scale=Scale, Dynamical_Model=Dynamical_Model, n_s=n_s, n_t=n_t, r_s=g_s, n_r=nrow(loc_g), loc_r=loc_g, alpha=alpha, beta=beta)
  SimList$DF[,"r_i"] = r_s[SimList$DF[,'s_i']]
}

# Extract count record
DF = SimList[["DF"]]


######################
# Estimate parameters 
######################

Options_vec = c("Include_movement"=Include_movement, "Dynamical_Model"=switch(Dynamical_Model,"Gompertz"=0,"Ricker"=1))
  # Include_movement: 0=No, 1=Stationary diffusion, 2=Independent-directions
  # Dynamical_Model: 0=Gompertz, 1=Ricker

Version = "movement"
dyn.unload( dynlib(Version) )
compile( paste0(Version,".cpp") )

# Data inputs
Data = list("Options_vec"=Options_vec, "n_i"=nrow(DF), "n_r"=MeshList$n_r, "n_g"=MeshList$n_g, "n_t"=n_t, "n_tdiv"=n_tdiv, "c_i"=DF[,'c_i'], "r_i"=DF[,'r_i']-1, "t_i"=DF[,'t_i']-1, "M1"=inla.as.dgTMatrix(MoveList[["M1"]]), "M2"=inla.as.dgTMatrix(MoveList[["M2"]]), "M3"=inla.as.dgTMatrix(MoveList[["M3"]]), "M4"=inla.as.dgTMatrix(MoveList[["M4"]]), "G0"=MeshList$spde_gmrf$param.inla$M0, "G1"=MeshList$spde_gmrf$param.inla$M1, "G2"=MeshList$spde_gmrf$param.inla$M2)
# Parameter inputs
Params = list("log_beta"=log(0.1), "alpha"=log(1), "log_mvec"=log(rep(0.2,length(mvec))), "logkappa"=log(1), "logSigmaU"=log(1), "logSigmaO"=log(1), "logmeanu0"=log(1), "ln_u_gt"=matrix( log(mean(Data$c_i)), nrow=MeshList$n_g, ncol=n_t), "Omegainput_g"=rnorm(MeshList$n_g))
Random = c("ln_u_gt", "Omegainput_g") # NULL #  
# Turn off params
Map = NULL
if( Options_vec["Include_movement"]==0 ){
  Map[["log_mvec"]] = factor( rep(NA,length(Params[["log_mvec"]])) )
  Params[["log_mvec"]] = rep(-1e20,length(Params[["log_mvec"]]))
}
if( Options_vec["Include_movement"]==1 ){
  Map[["log_mvec"]] = factor( rep(1,length(Params[["log_mvec"]])) )
  Params[["log_mvec"]] = rep(log(0.1),length(Params[["log_mvec"]]))
}
   
# Build object
dyn.load( dynlib(Version) )
Obj = MakeADFun( data=Data, parameters=Params, map=Map, random=Random, inner.control=list(maxit=2500))
Report = Obj$report() 
if( min(Report$u_rt)<0 ) stop("Error: Minimum u below 0!")

# Bounds
Min = function(List){ min(sapply(List,min)) }
Lwr = rep(-Inf, length(Obj$par))
  Lwr[grep("log_beta",names(Obj$par))] = log(0.01)    # Gompertz: rho = 1-beta
  Lwr[grep("log_mvec",names(Obj$par))] = log(1e-10)    
Upr = rep(Inf, length(Obj$par))
  Upr[grep("log_beta",names(Obj$par))] = log(1.99)     # Gompertz: rho = 1-beta

# Run                                                          # upper=50, lower=-50, 
Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1, getsd=TRUE, lower=Lwr, upper=Upr )
Report = Obj$report( )

# Get standard errors
ParHat = Obj$env$parList( Opt$par )

# Compare true and estimated values
cbind( "Est"=unlist(Report[c('rho','SigmaU','SigmaO','Range_raw','mvec')]), "True"=c(rho, SD_epsilon, SD_omega, Scale*2, mvec) )

# Compare movement matrices
image( log(ifelse(MoveMat_mesh<0.0001,NA,MoveMat_mesh)), main="True", zlim=c(log(0.0001),log(1)) ); dev.new()
for( tdiv in 1:Data$n_tdiv){
  if(tdiv==1) MoveMat_hat = Report$Mdiv_sparse
  if(tdiv>=2) MoveMat_hat = MoveMat_hat %*% Report$Mdiv_sparse
}
image( log(ifelse(as.matrix(MoveMat_hat)<0.0001,NA,MoveMat_mesh)), main="Estimated", zlim=c(log(0.0001),log(1)) )

# Compare movement from a given location
D0 = D0hat = rep(0,Data$n_r)
Start = sample(1:Data$n_r,size=1)
D0[Start] = D0hat[Start] = 1
D1 = MoveMat_mesh %*% D0
D1hat = MoveMat_hat %*% D0hat

# Plot movement from a given location
par(mfrow=c(1,2))
Col = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"))
ColBin = Col(50)[cut(as.vector(D1),seq(0,1,by=0.02)) ]
plot( x=MeshList$loc_r[,'x'], y=MeshList$loc_r[,'y'], col=ColBin, pch=20 )
plot( MeshList$mesh_domain, add=TRUE )
ColBin = Col(50)[cut(as.vector(D1hat),seq(0,1,by=0.02)) ]
plot( x=MeshList$loc_r[,'x'], y=MeshList$loc_r[,'y'], col=ColBin, pch=20 )
plot( MeshList$mesh_domain, add=TRUE )

