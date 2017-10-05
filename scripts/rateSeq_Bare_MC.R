#Set working directory
setwd("/Users/Farah/Desktop/R_scripts")
phiML<-read.table("MC_parameter_input.txt",header=TRUE,row.names=1)

#This is the likelihood function coded by Dan
#Include this from the other R script

#Input variables
#nu is a library-specific scaling factor
phiMLmc_df<-matrix(,nrow=dim(phiML)[1],ncol=3)
nu<-c(306.449350,327.589194,132.761826,121.591145,126.523919,130.656217,54.824759,3.197605,54.849455,71.841431,28.036492,88.699189,26.008729,34.891802,14.166423,25.540578)
nMC <- 1000 #number of Monte Carlo repetitions of synthetic experimental measurements at the time points above
nuMC<-rep(nu,nMC)
a0     <- aOFmHillParamsN20["a0"] #This and the following 2 lines are parameters in the equation for the dispersion
aInf   <- aOFmHillParamsN20["aInf"]
zHalf  <- aOFmHillParamsN20["zHalf"]
z_ss   <- phiML[1]
t <- c(5,   5,  10,  10,  15,  15,  25,  25,  35,  35,  50,  50,  90,  90, 160, 160) #experimental time points
tMC <- rep(t,nMC) #Just multiplying
lambda=0.30/60 # I added this line, Will growth rate affect number of timepoints and location? YES most likely...
n=3
#for (i in 1:dim(phiML)[1]){
for (i in 1:25){       #Test
  
  zmod<-phiML[i,1]*ExpSynWithDelay(tMC,alpha=phiML[i,2]+lambda, t_d=n/phiML[i,3], n=8) #Deterministic, why is n=8?
  plot(tMC,zmod,pch=16,xlab="Time (mins)")
  muY<-nuMC*zmod #Scale back to counts
  nH<-1
  a <- a0 +( aInf - a0 )* ( z_ss[i,1]^nH/(zHalf^nH + z_ss[i,1]^nH) )
  yMC<-rnbinom(n=nMC*length(t), mu=muY, size=a) #Randomly generated counts frome statistical model
  zMC <- yMC/nuMC #Convert random counts to corresponding amol estimates
  quartz()
  plot(tMC, zMC, type="p",main=i, pch=16, cex=0.5, cex.axis=0.5,xlab=NA, ylab=NA, ylim=c(0,1.1*max(zMC)), col="blue")
  points(tMC,zmod,pch=16,xlab="Time (mins)")
 
#test   
#t1_MC_vals<-NULL
#lengthz<-seq(1,1600,by=16)
#for (d in lengthz){
#  print (d)
#  t1_MC_vals[d]<-zMC[d]
#}

#hist(t1_MC_vals)
#Parameter "recovery"
x <- log(phiML) # log of starting guess of parameter values. Here, we can use as a guess
# the acutal values used to generate the random counts
exp.kinetics.out<-nlm(LLexpKINETICSlambda, p=as.numeric(x[1,]), y=yMC, lambda=0.30/60, n_delay=8, nu=nuMC,
                      aOFmHillParams=aOFmHillParamsN20, nH=1, t=tMC, print.level=0, fscale=200, iterlim=200)
#names(exp.kinetics.out)
#[1] "minimum"    "estimate"   "gradient"   "code"       "iterations"
exp.kinetics.out$estimate -> xest
phiMLmc <- print( exp(xest) ) #parameter estimation from Monte Carlo data
phiMLmc_df[i,]<-phiMLmc
tp <- seq(from=0, to=max(t), length.out=100)
ExpSyn.out <- phiMLmc[1]*ExpSynWithDelay(tp,alpha=phiMLmc[2]+lambda, t_d=n/phiMLmc[3], n=8)
lines(x=tp, y=ExpSyn.out, lwd=2, col='red', lty="solid")
legend(x="bottomright", legend=c("data", "ML"), lty=c(NA,"solid"), lwd=c(NA,2), pch=c(16,NA),
       col=c("black","red"), bty="o", cex=1.5, pt.cex=c(0.5, NA))

}
#Comparison of actual and "recovered" parameters
phiMLmc_df<-data.frame(phiMLmc_df)
merged_df<-cbind(phiML,phiMLmc_df)
head(merged_df,n=8)
#                [,1]        		[,2]     					[,3]
#phiML      2310.863 		0.006060686 		1.744571
#phiMLmc 2371.503 		0.006049316 		1.779578

#Dan mentioned not to use really high alpha values, is the reason that convergence doesn't take place anymore?
#Do I need to now use the predict function to get all the values for alpha? and subtract those from the mu then square and add them?
