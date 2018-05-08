library(TMB)

setwd("/Users/jiecao/Desktop/hw3_soln")
CPUE = FishData::download_catch_rates( survey="Eastern_Bering_Sea", species_set="Gadus chalcogrammus", error_tol=0.01, localdir=paste0(getwd(),"/") )
B_t = tapply( CPUE[,'Wt'], INDEX=CPUE[,'Year'], FUN=mean )
log_yt <- log(B_t)

#Run TMB model 
compile("gompertz.cpp")
dyn.load(dynlib("gompertz"))

N = length(log_yt)	
data <- list(log_yt = log_yt,model_switch=1)
parameters <- list(pop_par=c(1,0.5), log_sigma = -1,log_x0 = 4, log_xt = rep(mean(log_yt), N))

#data <- list(log_yt = log_yt,model_switch=2)
#parameters <- list(pop_par=c(0.6,100), log_sigma = -1,log_x0 = 4, log_xt = rep(mean(log_yt), N))

obj <- MakeADFun(data, parameters, random = "log_xt", DLL = "gompertz")
obj$hessian <- FALSE
opt <- do.call("optim", obj)
sd <- sdreport(obj)

# extract fixed effects:
fixed <- summary(sd, "fixed")
# extract estimated process:
log_xt <- summary(sd, "random")[, "Estimate"]
log_xt_se <- summary(sd, "random")[, "Std. Error"]

# simulation
n_sim = 100
sim <- replicate(n_sim, {
simdata <- obj$simulate()
obj2 <- MakeADFun(list(log_yt=c(simdata$log_yt[1:31],rep(NA,5)),model_switch=1), parameters, random = "log_xt", DLL = "gompertz", silent = TRUE)
TMBhelper::Optimize( obj=obj2, newtonsteps=1)
c( true = simdata$log_xt[32:36],
   est = summary(sdreport(obj2), "random")[, "Estimate"][32:36],
   sd = summary(sdreport(obj2), "random")[, "Std. Error"][32:36])
})

# calculate coverage 
cal_coverage <- function (true, est, sd, interval){
  low = est - qnorm(interval/2) * sd
  high = est + qnorm(interval/2) * sd
  return(true > low & true < high)
}

gom.matrix = matrix(NA,nrow=n_sim,ncol=5)
for(i in 1:n_sim){
  for(t in 1:5){
    gom.matrix[i,t]=cal_coverage(sim[t,i],sim[t+5,i],sim[t+10,i],0.5)
  }
}
apply(gom.matrix,2,sum)

# ricker 
parameters2 <- list(pop_par=c(0.6,100), log_sigma = -1,log_x0 = 4, log_xt = rep(mean(log_yt), N))
sim2 <- replicate(n_sim, {
simdata <- obj$simulate()
obj3 <- MakeADFun(list(log_yt=c(simdata$log_yt[1:31],rep(NA,5)),model_switch=2), parameters2, random = "log_xt", DLL = "gompertz", silent = TRUE)
TMBhelper::Optimize( obj=obj3, newtonsteps=1)
c( true = simdata$log_xt[32:36],
   est = summary(sdreport(obj3), "random")[, "Estimate"][32:36],
   sd = summary(sdreport(obj3), "random")[, "Std. Error"][32:36])
})

ricker.matrix = matrix(NA,nrow=n_sim,ncol=5)
for(i in 1:n_sim){
  for(t in 1:5){
    ricker.matrix[i,t]=cal_coverage(sim2[t,i],sim2[t+5,i],sim2[t+10,i],0.5)
  }
}
apply(ricker.matrix,2,sum)
