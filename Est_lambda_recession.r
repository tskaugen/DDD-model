# This script estimates the gamma parameters of Lambda and lambda and plots the CDF of Lambda and lambda . Note that lambda written with capital or lowercase l is significant (see Skaugen and Onof, 2014)
#.trPaths <- paste(paste(Sys.getenv('APPDATA'), '\\Tinn-R\\tmp\\', sep=''),c('', 'search.txt', 'objects.txt', 'file.r','selection.r', 'block.r', 'lines.r'), sep='') 

rm(list=ls())
catchment <- "2.11"
q_file<- paste(catchment,"_24hptq_kal.txt", sep="")
path <- paste("f:\\HB\\HB-modellering\\DDDtestbenk\\DDD_pakke\\",sep="")
out_path <- paste(path, sep="")
outfile = paste(out_path,"Lam_24h_kal_",catchment,".png",sep="")

ptq.file <-paste(path,"indata\\",q_file,sep="")#Observations of runoff

ptqinn <- read.table(ptq.file)
navn <- list("yr","mnt","day","hr","p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","q")
names(ptqinn)[1:25] <-navn[1:25]
days <- length(ptqinn$yr)
meanQ <- mean(ptqinn$q)

Lfrac=0

par(mfrow=c(4,2))
  
lambda <- vector("numeric",(days-1))
komplambda<-  vector("numeric",(days-1))
lambda[1:(days-1)] <-0.0      #recession
komplambda[1:(days-1)] <- 1.0 #Complementary values (not recession)

for (i in 1:(days-1))
{
 if(ptqinn$q[i]==0)ptqinn$q[i] <- 0.001
 if(ptqinn$q[i+1]==0)ptqinn$q[i+1] <- 0.001
# if( (log(ptqinn$q[i+1])-log(ptqinn$q[i])) > 0.0) freqv[i] <-log(ptqinn$q[i+1])-log(ptqinn$q[i])
 lambda[i] <-log(ptqinn$q[i])-log(ptqinn$q[i+1])# sjekker om dagen verdi er en resesjonsverdi
 #if we "log" tranform, we look for positive valueswhich are recession events
 if(lambda[i] > 0.0) komplambda[i] <- NA # assigns NA for all days with recession
 if(lambda[i] <= 0.0) lambda[i] <- NA # assigns NA for all days NOT with recession
}

recess <-na.action(na.omit(komplambda))# gives us vector positions in time series for all events with recession
nylambda <-na.omit(lambda)# equal length as "recess"

plot(nylambda,ptqinn$q[recess],col="blue", xlab="Lambda", ylab="Q(t)")#Plots sampled Lambda against runoff at time t
mtext(paste("Lambda plots:",ptq.file),line=0,cex=0.6)

bins <- seq(0,(max(nylambda)+0.2),0.001) 
tull2 <-hist(nylambda, freq=FALSE, TRUE,plot=TRUE, xlim=c(0.0,(max(nylambda)+ 0.2)), breaks =bins) #plots histograms of Lambda

sumlambda <-vector("numeric", length(tull2$density))

for (i in 1:length(tull2$density)) sumlambda[i] <-sum(tull2$density[1:i]*(bins[2]-bins[1]))

plot(tull2$mids,sumlambda, type="l", lwd=2.0,xlab="Lambda",ylab="Fraction", ylim=c(0.0,1.0)) #plots empirical CDF og Lambda 

dy <- sumlambda
dx <- tull2$mids
nlmodexp <-nls(dy~(1-exp(-(dx/p1))),start=list(p1=mean(nylambda)), trace=T)
par1exp <- round(summary(nlmodexp)$coefficients[1],4)
feilexp <- round(summary(nlmodexp)$sigma,4)

#lines(dx, 1-exp(-(dx/par1exp)), col="black")
#lines (dx,(1/par1exp)*exp(-(dx/par1exp)),col=5)
nlmodgamma <- nls(dy~pgamma(dx,p2,p3),start=list(p2=par1exp,p3=par1exp), trace=T) # uses exponential parameter (par1exp) as startingpoint for fitting
par2gamma <- round(summary(nlmodgamma)$coefficients[1],3)#Shape
par3gamma <- round(1/summary(nlmodgamma)$coefficients[2],3)#Scale
feilgamma <- round(summary(nlmodgamma)$sigma,4)
lines(dx, pgamma(dx,par2gamma,1/par3gamma), col= "red", lwd=1.5)
mtext("CDF Lambda", line=2,cex=1.0)
#mtext(paste("Exp_par=",par1exp, "ResSE_exp=",feilexp),line=0,cex=0.8)
mtext(paste("Shape=",par2gamma, "Scale=",par3gamma,"ResSE_gam=",feilgamma ),line=1,cex=0.8)
#legend(0.2,0.2, "-- exponential",bty="n",text.col=2)
legend(0.2,0.4, "--gamma", bty="n", text.col="red")

arg <- seq(0.1,0.9,0.1)
arg[10] <-0.99
lam <-vector("numeric",10)# number of levels
avglam <-vector("numeric",10)# number of levels
glam <- vector("numeric",10)# number of levels
lam[1:10] <-qexp(arg,1/par1exp)# Estimated from runoff data, exponential model, NOT good
glam[1:10] <-qgamma(arg,par2gamma,1/par3gamma) 

avglam <-glam # Estimated from runoff datat, gamma model

trlam <-vector("numeric",10)  # Lowercase lambda for level
trlam[1] <- avglam[1]
Llam <- avglam[1]
Lwgt <- Lfrac

len <-20 # length of UH  for comparison

avgUH<- vector("numeric", len)
lagUH1<- vector("numeric", len)
lagUH2<- vector("numeric", len)
lagUH3<- vector("numeric", len)
lagUH4<- vector("numeric", len)
lagUH5 <- vector("numeric", len)
lagUH6 <- vector("numeric", len)
lagUH7<- vector("numeric", len)
lagUH8 <- vector("numeric", len)
lagUH9<- vector("numeric", len)
lagUH10 <- vector("numeric", len)


for(lag in 1:10)#Loekke for alle lag
{

sumwgt <-sum(avglam[1:lag])
wgt <- vector("numeric", lag)

for (i in 1:lag)
{
wgt[i] <-avglam[i]/sumwgt
}
wgt <-wgt*(1-Lfrac)
wgt[1] <- wgt[1]+Lfrac

avgUH[1:len] <-avglam[lag]*exp(-avglam[lag]*(1:len))# den vi skal sammenlikne mot
dy <-avgUH
dx <-seq(1,len,1)  

#lag1
if(lag==1)
{
lagUH1[1:len] <-avglam[1]*exp(-avglam[1]*(1:len))#avglam[1] og trlam[1]og Llam antas å være like
}

#lag2
if(lag==2)
{
nlmodexp <-nls(dy~((trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (p1*exp(-p1*(dx))*wgt[2])),start=list(p1=0.1), trace=T)
trlam[2] <- round(summary(nlmodexp)$coefficients[1],4)
lagUH2[1:len] <- (trlam[1]*exp(-trlam[1]*(1:len))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(1:len))*wgt[2])
#plot(avgUH, type="l",col="black", ylim=c(0, max(avgUH)))
#lines(lagUH2,col="blue")
#mtext(paste("lambda=",trlam[lag], "Level no=",lag, " Blue estimated LHS eq. 14., Black, observed RHS Eq. 14"),line=0,cex=0.8)
}
#lag3
if(lag==3)
{
nlmodexp <-nls(dy~((trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(dx))*wgt[2])
+ (p1*exp(-p1*(dx))*wgt[3])),start=list(p1=trlam[lag-1]), trace=T)
trlam[3] <- round(summary(nlmodexp)$coefficients[1],4)
lagUH3[1:len] <- (trlam[1]*exp(-trlam[1]*(1:len))*wgt[1])+(trlam[2]*exp(-trlam[2]*(1:len))*wgt[2])+ (trlam[3]*exp(-trlam[3]*(1:len))*wgt[3])
plot(avgUH, type="l",col="black", ylim=c(0, max(avgUH)))
lines(lagUH3,col="blue")
mtext(paste("lambda=",trlam[lag], "Level no=",lag, " Blue estimated LHS eq. 14., Black, observed RHS Eq. 14"),line=0,cex=0.8)
}
#lag4
if(lag==4)
{
nlmodexp <-nls(dy~((trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(dx))*wgt[2])
+ (trlam[3]*exp(-trlam[3]*(dx))*wgt[3])+(p1*exp(-p1*(dx))*wgt[4])),start=list(p1=trlam[lag-1]), trace=T)
trlam[4] <- round(summary(nlmodexp)$coefficients[1],4)
lagUH4[1:len] <- ((trlam[1]*exp(-trlam[1]*(1:len))*wgt[1])+(trlam[2]*exp(-trlam[2]*(1:len))*wgt[2])+ (trlam[3]*exp(-trlam[3]*(1:len))*wgt[3])
+ (trlam[4]*exp(-trlam[4]*(1:len))*wgt[4]))
#plot(avgUH, type="l",col="black", ylim=c(0, max(avgUH)))
#lines(lagUH4,col="blue")# Estimert
#mtext(paste("lambda=",trlam[lag], "Level no=",lag, "Blue estimated LHS eq. 14., Black, observed RHS Eq. 14"),line=0,cex=0.8)
}
#lag5
if(lag==5)
{
nlmodexp <-nls(dy~((trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(dx))*wgt[2])
+ (trlam[3]*exp(-trlam[3]*(dx))*wgt[3])+(trlam[4]*exp(-trlam[4]*(dx))*wgt[4])+(p1*exp(-p1*(dx))*wgt[5])),start=list(p1=trlam[lag-1]), trace=T)
trlam[5] <- round(summary(nlmodexp)$coefficients[1],4)
lagUH5[1:len] <- ((trlam[1]*exp(-trlam[1]*(1:len))*wgt[1])+(trlam[2]*exp(-trlam[2]*(1:len))*wgt[2])+ (trlam[3]*exp(-trlam[3]*(1:len))*wgt[3])
+ (trlam[4]*exp(-trlam[4]*(1:len))*wgt[4])+(trlam[5]*exp(-trlam[5]*(1:len))*wgt[5]))
#plot(avgUH, type="l",col="black", ylim=c(0, max(avgUH)))
#lines(lagUH5,col="blue")
#mtext(paste("lambda=",trlam[lag], "Level no=",lag, " Blue estimated LHS eq. 14., Black, observed RHS Eq. 14"),line=0,cex=0.8)
}
#lag6
if(lag==6)
{
nlmodexp <-nls(dy~((trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(dx))*wgt[2])
+ (trlam[3]*exp(-trlam[3]*(dx))*wgt[3])+(trlam[4]*exp(-trlam[4]*(dx))*wgt[4])+(trlam[5]*exp(-trlam[5]*(dx))*wgt[5])
+ (p1*exp(-p1*(dx))*wgt[6])),start=list(p1=trlam[lag-1]), control=list(maxiter =100, tol=1e-02),trace=T)
trlam[6] <- round(summary(nlmodexp)$coefficients[1],4)
lagUH6[1:len] <- ((trlam[1]*exp(-trlam[1]*(1:len))*wgt[1])+(trlam[2]*exp(-trlam[2]*(1:len))*wgt[2])+ (trlam[3]*exp(-trlam[3]*(1:len))*wgt[3])
+ (trlam[4]*exp(-trlam[4]*(1:len))*wgt[4])+(trlam[5]*exp(-trlam[5]*(1:len))*wgt[5])+(trlam[6]*exp(-trlam[6]*(1:len))*wgt[6]))
plot(avgUH, type="l",col="black", ylim=c(0, max(avgUH)))
lines(lagUH6,col="blue")# Estimert
mtext(paste("lambda=",trlam[lag], "Level no=",lag, "Blue estimated LHS eq. 14., Black, observed RHS Eq. 14"),line=0,cex=0.8)
}
#lag7
if(lag==7)
{
nlmodexp <-nls(dy~((trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(dx))*wgt[2])
+ (trlam[3]*exp(-trlam[3]*(dx))*wgt[3])+(trlam[4]*exp(-trlam[4]*(dx))*wgt[4])+(trlam[5]*exp(-trlam[5]*(dx))*wgt[5])
+ (trlam[6]*exp(-trlam[6]*(dx))*wgt[6])+(p1*exp(-p1*(dx))*wgt[7])),start=list(p1=trlam[lag-1]), control=list(maxiter =100, tol=1e-02),trace=T)
trlam[7] <- round(summary(nlmodexp)$coefficients[1],4)
lagUH7[1:len] <- ((trlam[1]*exp(-trlam[1]*(1:len))*wgt[1])+(trlam[2]*exp(-trlam[2]*(1:len))*wgt[2])+ (trlam[3]*exp(-trlam[3]*(1:len))*wgt[3])
+ (trlam[4]*exp(-trlam[4]*(1:len))*wgt[4])+(trlam[5]*exp(-trlam[5]*(1:len))*wgt[5])+(trlam[6]*exp(-trlam[6]*(1:len))*wgt[6])
+ (trlam[7]*exp(-trlam[7]*(1:len))*wgt[7]))
#plot(avgUH, type="l",col="black", ylim=c(0, max(avgUH)))
#lines(lagUH7,col="blue")
#mtext(paste("Exp_par=",trlam[lag], "Lag no=",lag),line=0,cex=0.8)
}
#lag8
if(lag==8)
{
nlmodexp <-nls(dy~((trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(dx))*wgt[2])
+ (trlam[3]*exp(-trlam[3]*(dx))*wgt[3])+(trlam[4]*exp(-trlam[4]*(dx))*wgt[4])+(trlam[5]*exp(-trlam[5]*(dx))*wgt[5])
+ (trlam[6]*exp(-trlam[6]*(dx))*wgt[6])+(trlam[7]*exp(-trlam[7]*(dx))*wgt[7])+(p1*exp(-p1*(dx))*wgt[8])),start=list(p1=trlam[lag-1]), control=list(maxiter =100, tol=1e-02),trace=T)
trlam[8] <- round(summary(nlmodexp)$coefficients[1],4)
lagUH8[1:len] <- ((trlam[1]*exp(-trlam[1]*(1:len))*wgt[1])+(trlam[2]*exp(-trlam[2]*(1:len))*wgt[2])+ (trlam[3]*exp(-trlam[3]*(1:len))*wgt[3])
+ (trlam[4]*exp(-trlam[4]*(1:len))*wgt[4])+(trlam[5]*exp(-trlam[5]*(1:len))*wgt[5])+(trlam[6]*exp(-trlam[6]*(1:len))*wgt[6])
+ (trlam[7]*exp(-trlam[7]*(1:len))*wgt[7])+(trlam[8]*exp(-trlam[8]*(1:len))*wgt[8]))
plot(avgUH, type="l",col="black", ylim=c(0, max(avgUH)))
lines(lagUH8,col="blue")# Estimert
mtext(paste("lambda=",trlam[lag], "Level no=",lag, "Blue estimated LHS eq. 14., Black, observed RHS Eq. 14"),line=0,cex=0.8)
}
#lag9
if(lag==9)
{
print(paste("trlam[8]=",trlam[8]))
nlmodexp <-nls(dy~((trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(dx))*wgt[2])
+ (trlam[3]*exp(-trlam[3]*(dx))*wgt[3])+(trlam[4]*exp(-trlam[4]*(dx))*wgt[4])+(trlam[5]*exp(-trlam[5]*(dx))*wgt[5])
+ (trlam[6]*exp(-trlam[6]*(dx))*wgt[6])+(trlam[7]*exp(-trlam[7]*(dx))*wgt[7])+(trlam[8]*exp(-trlam[8]*(dx))*wgt[8])
+ (p1*exp(-p1*(dx))*wgt[9])),start=list(p1=trlam[lag-1]),control=list(maxiter =100, tol=1e-02), trace=T)
trlam[9] <- round(summary(nlmodexp)$coefficients[1],4)
lagUH9[1:len] <- ((trlam[1]*exp(-trlam[1]*(1:len))*wgt[1])+(trlam[2]*exp(-trlam[2]*(1:len))*wgt[2])+ (trlam[3]*exp(-trlam[3]*(1:len))*wgt[3])
+ (trlam[4]*exp(-trlam[4]*(1:len))*wgt[4])+(trlam[5]*exp(-trlam[5]*(1:len))*wgt[5])+(trlam[6]*exp(-trlam[6]*(1:len))*wgt[6])
+ (trlam[7]*exp(-trlam[7]*(1:len))*wgt[7])+(trlam[8]*exp(-trlam[8]*(1:len))*wgt[8])+(trlam[9]*exp(-trlam[9]*(1:len))*wgt[9]))
#plot(avgUH, type="l",col="black", ylim=c(0, max(avgUH)))
#lines(lagUH9,col="blue")
#mtext(paste("Exp_par=",trlam[lag], "Lag no=",lag),line=0,cex=0.8)
}
#lag10
if(lag==10)
{
nlmodexp <-nls(dy~((trlam[1]*exp(-trlam[1]*(dx))*wgt[1])+ (trlam[2]*exp(-trlam[2]*(dx))*wgt[2])
+ (trlam[3]*exp(-trlam[3]*(dx))*wgt[3])+(trlam[4]*exp(-trlam[4]*(dx))*wgt[4])+(trlam[5]*exp(-trlam[5]*(dx))*wgt[5])
+ (trlam[6]*exp(-trlam[6]*(dx))*wgt[6])+(trlam[7]*exp(-trlam[7]*(dx))*wgt[7])+(trlam[8]*exp(-trlam[8]*(dx))*wgt[8])
+ (trlam[9]*exp(-trlam[9]*(dx))*wgt[9])+(p1*exp(-p1*(dx))*wgt[10])),start=list(p1=trlam[lag-1]),control=list(maxiter =100, tol=1e-02), trace=T)
trlam[10] <- round(summary(nlmodexp)$coefficients[1],4)
lagUH10[1:len] <- ((trlam[1]*exp(-trlam[1]*(1:len))*wgt[1])+(trlam[2]*exp(-trlam[2]*(1:len))*wgt[2])+ (trlam[3]*exp(-trlam[3]*(1:len))*wgt[3])
+ (trlam[4]*exp(-trlam[4]*(1:len))*wgt[4])+(trlam[5]*exp(-trlam[5]*(1:len))*wgt[5])+(trlam[6]*exp(-trlam[6]*(1:len))*wgt[6])
+ (trlam[7]*exp(-trlam[7]*(1:len))*wgt[7])+(trlam[8]*exp(-trlam[8]*(1:len))*wgt[8])+(trlam[9]*exp(-trlam[9]*(1:len))*wgt[9])
+ (trlam[10]*exp(-trlam[10]*(1:len))*wgt[10]))
plot(avgUH, type="l",col="black", ylim=c(0, max(avgUH)))
lines(lagUH10,col="blue")# Estimert
mtext(paste("lambda=",trlam[lag], "Level no=",lag, "Blue estimated LHS eq. 14., Black, observed RHS Eq. 14"),line=0,cex=0.8)
}

}#loop for levels

#Estimate parameters for trlam (lambda)
dy <- trlam
dx <- arg
#exponentiell
nlmodexp <-nls(dx~(1-exp(-(dy/p2))),start=list(p2=par1exp), trace=T)
par2exp <- round(summary(nlmodexp)$coefficients[1],4)
feilexp <- round(summary(nlmodexp)$sigma,4)
#Gamma
nlmodgamma3 <- nls(dx~pgamma(dy,p2,p3),start=list(p2=par1exp,p3=par1exp), trace=T)
par2gamma2 <- round(summary(nlmodgamma3)$coefficients[1],3)#Shape
par3gamma2 <- round(1/summary(nlmodgamma3)$coefficients[2],3)#Scale
feilgamma2 <- round(summary(nlmodgamma3)$sigma,4)

plot(trlam,dx, xlab="lambda",ylab="Fraction",type="l", lwd=2.0, ylim=c(0.0,1.0))  
mtext(paste("CDF lambda, MAD=", round(mean(ptqinn$q),3)), line=2,cex=1.0)
lines(trlam, pgamma(trlam,par2gamma2,1/par3gamma2), col="red") #gamma model
mtext(paste("Shape=",par2gamma2, "Scale=",par3gamma2,"ResSE_gam=",feilgamma2 ),line=1,cex=0.8)

legend(0.2,0.4, "--gamma", bty="n", text.col="red")

savePlot(filename = outfile, type = c("png"),
            device = dev.cur(), restoreConsole = TRUE)
            
