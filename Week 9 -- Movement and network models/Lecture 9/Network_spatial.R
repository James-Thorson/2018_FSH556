
library(TMB)
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 9 -- Movement and network models/Lecture 9" )

# Network parentage
set.seed(1)
parent_b = c( NA, 1, 2, 3, 3, 5, 5, 5, 8, 9, 10 )
dist_b = c( NA, rlnorm(length(parent_b)-1, meanlog=0, sdlog=1) )
DF_b = cbind( parent_b=parent_b, child_b=1:length(parent_b), dist_b=dist_b )
b_i = 1:nrow(DF_b)

# Distance matrix
D = array(NA, dim=rep(length(dist_b),2) )
D[ na.omit(DF_b)[,c("parent_b","child_b")] ] = na.omit(DF_b)[,'dist_b']
D[ na.omit(DF_b)[,c("child_b","parent_b")] ] = na.omit(DF_b)[,'dist_b']
Dsparse = as(ifelse(is.na(D),0,D),"dsCMatrix")

sparse_exp = function(A){
  B = exp(as.matrix(A))
  B = ifelse( B==1, 0, B )
  return( as(B,"dsCMatrix") )
}
exp_neg_D = sparse_exp(-D)

# Parameter values
Rho=0.3           # Correlation at distance = 1
SDmarg=0.3
theta = -log(Rho)  # Decorrelation rate per distance
log_mean = 4
detect_prob = 0.5
Npass = 3

# Detection probability for each pass
detect_per_pass = rep(NA,Npass)
for( zI in 1:Npass ){
  detect_per_pass[zI] = detect_prob * (1-detect_prob)^(zI-1)
}

# Derived
SDinput = SDmarg * sqrt(2*theta)
SDcond = SDmarg * sqrt(1-Rho^2)

# Simulate random effect
omega_b = rep(NA, length(b_i))
pow = function(a,b) a^b
Type = function(a) a
rho_b = SDinput_b = rep(NA, length(b_i))
for( b in 1:length(omega_b) ){
  if( is.na(dist_b[b]) ){
    # Correlation between i and parent(i) as distance -> INF
    rho_b[b] = 0;
    # SD of Ornstein-Uhlenbeck process as distance -> INF
    SDinput_b[b] = SDinput / pow(2*theta, 0.5);
    # conditional probability
    omega_b[b] = rnorm(n=1, Type(0.0), SDinput_b[b] );
  }
  if( !is.na(dist_b[b]) ){
    # Correlation between i and parent(i)
    rho_b[b] = exp(-theta * dist_b[b]);
    # SD of O-U process
    SDinput_b[b] = pow( pow(SDinput,2)/(2*theta) * (1-exp(-2*theta*dist_b[b])), 0.5 );
    # conditional probability
    omega_b[b] = rnorm(n=1, rho_b[b]*omega_b[parent_b[b]], SDinput_b[b] );
  }
}

# Simulate data
lambda_i = exp( omega_b + log_mean )
c_iz = matrix( rpois(n=length(lambda_i)*Npass, lambda=outer(lambda_i,detect_per_pass[1:Npass])), ncol=Npass )

# Plot detectability
if( Npass>1 ){
  prop_iz = c_iz / outer( rowSums(c_iz), rep(1,Npass) )
  png( file=paste0("Detectability.png"), width=4, height=4, res=200, units="in" )
    par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    matplot( y=t(prop_iz), type="l", col="blue", lty="dotted", log="y", xlab="Pass", ylab="Detection probability / proportion" )
    lines( detect_per_pass/sum(detect_per_pass), lwd=3 )
    legend( "topright", legend=c("True","Observed"), fill=c("black","blue"), bty="n" )
  dev.off()
}

# Format inputs for TMB
Data = list( "Options_vec"=c(0), "c_iz"=c_iz, "b_i"=b_i-1, "parent_b"=DF_b[,'parent_b']-1, "child_b"=DF_b[,'child_b']-1, "dist_b"=DF_b[,'dist_b'] )
Params = list( "log_theta"=log(theta), "log_SD"=log(SDinput), "log_mean"=log(1), "logit_detect_prob"=0, "Epsiloninput_b"=rep(0,max(b_i)) )
Random = "Epsiloninput_b"

# Turn off detect_prob if single-pass model
Map = list()
if( Npass==1 ){
  Map[["logit_detect_prob"]] = factor(NA)
  Params[["logit_detect_prob"]] = qlogis(0.9999)
}

# Exclude data from "join" nodes
c_iz[c(3,5),] = NA

# Compile TMB
compile( "network_spatial.cpp" )
dyn.load( dynlib("network_spatial") )

##################
# Check that precision matrix is correct
##################

# Extract from TMB
Data[["c_iz"]][] = NA
Obj = MakeADFun( data=Data, parameters=Params, random=Random )
Hess = Obj$env$spHess( random=TRUE )

# Calculate analytically
Dfilled = ifelse(is.na(D), Inf, D)
Q = -exp(-theta*Dfilled) / (SDmarg*sqrt(1-exp(-theta*Dfilled)^2))^2
diag(Q) = 1/SDmarg^2 + colSums(exp(-2*theta*Dfilled) / (SDmarg*sqrt(1-exp(-theta*Dfilled)^2))^2)

# Compare
summary( as.vector(Q) - as.vector(Hess) )

##############
# Compare estimates from two methods
##############

######### Version 0:  Sweep upstream to downstream
Data = list( "Options_vec"=c(0), "c_iz"=c_iz, "b_i"=b_i-1, "parent_b"=DF_b[,'parent_b']-1, "child_b"=DF_b[,'child_b']-1, "dist_b"=DF_b[,'dist_b'] )
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map )
Opt0 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
Report0 = Obj$report()

######### Version 1:  Joint probability
Data = list( "Options_vec"=c(1), "c_iz"=c_iz, "b_i"=b_i-1, "parent_b"=DF_b[,'parent_b']-1, "child_b"=DF_b[,'child_b']-1, "dist_b"=DF_b[,'dist_b'] )
Obj = MakeADFun( data=Data, parameters=Params, random=Random, map=Map )
Opt1 = TMBhelper::Optimize( obj=Obj, newtonsteps=1 )
Report1 = Obj$report()

# Compare MLE from two methods
True = c(log(theta),SDinput,log_mean)
if( Npass>1 ) True = c(True, qlogis(detect_prob))
rbind( "True"=True, "Method1"=Opt0$par, "Method2"=Opt1$par )

# Compare timings
c( Opt0$run_time, Opt1$run_time )

# Bias-correcting estimates
Opt1 = TMBhelper::Optimize( obj=Obj, newtonsteps=1, bias.correct=TRUE )
summary( Opt1$SD, "report" )


