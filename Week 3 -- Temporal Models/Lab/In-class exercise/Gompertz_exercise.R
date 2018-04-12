
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2018_FSH556/Week 3 -- Temporal Models/Lab/In-class exercise" )

######################
# Modify to add sampling variance to measurement error
######################

# Modify to decompose measurement error into sampling and additional error
sd_B_t = sqrt( tapply( CPUE[,'Wt'], INDEX=CPUE[,'Year'], FUN=var ) / tapply( CPUE[,'Wt'], INDEX=CPUE[,'Year'], FUN=length ) )
CV_b_t = sd_B_t / B_t
SD_log_b_t = sqrt(log( CV_b_t^2 + 1 ))

# Compile model
Version = "gompertz_2"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "nt"=length(B_t), "log_b_t"=log(B_t), "SD_log_b_t"=SD_log_b_t )
Parameters = list( "log_d0"=0, "log_sigmaP"=1, "log_sigmaM"=1, "alpha"=0, "rho"=0, "log_d_t"=rep(0,Data$nt) )
Random = c("log_d_t")
if( Use_REML==TRUE ) Random = union( Random, c("log_d0","alpha","rho") )

# Build object
dyn.load( dynlib("gompertz_2") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, DLL="gompertz_2")  #
Opt2 = TMBhelper::Optimize( obj=Obj )

# Plot again
for( z in 1:2 ){
  png( file=paste0("pollock_",z,"-with_measurement_error.png"), width=8, height=5, res=200, units="in" )
    par( mar=c(3,3,1,1), mgp=c(2,0.5,0), tck=-0.02 )
    plot( x=1:Data$nt, y=Data$log_b_t, col="blue", cex=1.2, xlab="Time", ylab="Value", pch=20 )
    if(z>=2){
      points( x=1:Data$nt, y=as.list(Opt2$SD,"Estimate")$log_d_t, col="red" )
      for( t in 1:Data$nt) lines( x=rep(t,2), y=as.list(Opt2$SD,"Estimate")$log_d_t[t]+c(-1.96,1.96)*as.list(Opt2$SD,"Std. Error")$log_d_t[t], col="red" )
    }
  dev.off()
}

write.csv( cbind(Opt2$diagnostics, summary(Opt$SD,"fixed")[,'Std. Error']), file=paste0(getwd(),"/pollock_2.csv") )

###############
# Fix SigmaM at zero
###############

Map = list()
Map$log_sigmaM = factor(NA)
Parameters$log_sigmaM = log(1e-10)
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random, map=Map, DLL="gompertz_2")  #
Opt3 = TMBhelper::Optimize( obj=Obj )
Report = Obj$report()

# Compare AIC
c(Opt2$AIC, Opt3$AIC)
