##This program takes a file "Bottle_20c.txt" of newick local genealogies (with n>2)
#and estimates the effective population size trajectory. The code for n=2 is available upon request 
#contact: julia.pal.r@gmail.com
#Reference: Bayesian Nonparametric Inference of Population Size Changes from Sequential Genealogies
########### by Julia A Palacios, John Wakeley, Sohini Ramachandran
#DOI:http://dx.doi.org/10.1101/019216

library("ape")
path<-getwd()
lib<-paste(path,"/data/",sep="")
lib2<-paste(path,"/source/",sep="")
source(paste(lib2,"coal_lik_BSMC.R",sep=""))
MyTree <- read.tree(paste(lib,"Bottle_20c.txt",sep=""))



##Data preparation
n<-20 #number of samples
sim<-length(MyTree) #Total number of trees to consider. It can be changed to a smaller number
sim<-50 #For this example, the first 50 genealogies
scaling<-10
D<-read_times(MyTree,n,sim,scaling)
#we scale times by 10 to avoid numerical problems

tol<-.00001 #tolerance factor to detect diference between branch lengths

##Grid definition
summary(D[,dim(D)[2]])
window<-max(D)+.0001
grid.size<-100
grid<-seq(0,window,length.out=grid.size)
grid<-c(grid,max(D)+.0002)

  
#Data setup: Reads all sufficient statistics in matrix format
info<-find_info2(MyTree,D,sim,n,tol,scaling)


###Bayesian MCMC Algorithm
set.seed(2014)
alg=1
TrjL=1; Nleap=15; stepsz=TrjL/Nleap
stepsz<-.01
#NSAMP=50000; NBURNIN=1000 #Original parameters used in the paper
NSAMP=500; NBURNIN=50
Ngrid<-length(grid)-1
# hyperparameter in prior of tau
alpha=1e-3; beta=1e-3
#MCMC sampling preparation
SAMP=list(2)
SAMP[[1]]=matrix(NA,NSAMP-NBURNIN,Ngrid) # transformed effective population size
SAMP[[2]]=rep(NA,NSAMP-NBURNIN) # precision parameter in Brownian motion
acpi=0;acpt=0
PRINT<-T
#Initial Values
f_init=rep(0.5,Ngrid)
theta <- c(log(f_init),-1.6)+.0001 
Nleap<-15
stepsz<-.1
alldata<-get.data(grid,D,n,coal_lik_init,info$info_times,info$Fl,info$latent,info$t_new,info$t_del)
U<-function(theta,grad=F)U_split_smc(theta,alldata$lik_init,alldata$invC,alpha,beta,grad)
current.u<-U(theta,F)
current.grad<-U(theta,T)

start_time = Sys.time()
for(Iter in 1:NSAMP){
  
  if(PRINT&&Iter%%50==0){
    cat(Iter, ' iterations have been finished!\n' )
    cat('Online acceptance rate is ',acpi/50,'\n')
    acpi=0
  }
  res=eval(parse(text='splitHMC'))(theta,function(theta,grad=F)U_split_smc(theta,alldata$lik_init,alldata$invC,alpha,beta,grad),alldata$rtEV,alldata$EVC,stepsz,Nleap,current.u,current.grad)
  theta=res$q;
  current.u<-res$current.u
  current.grad<-res$current.grad
  N<-exp(theta[1:(length(theta)-1)])
  acpi=acpi+res$Ind
  if(Iter>NBURNIN){
    SAMP[[1]][Iter-NBURNIN,]<-theta[1:(length(theta)-1)]
    SAMP[[2]][Iter-NBURNIN]<-theta[length(theta)]
    acpt<-acpt+res$Ind
  }
}
stop_time = Sys.time()
time=stop_time-start_time
cat('\nTime consumed : ',time)




##Bayesian Summary of N(t)

#med=apply(exp(SAMP[[1]])[1:(Iter-NBURNIN-1),],2,median); #se=apply(exp(SAMP[[1]]),2,sd)/sqrt(dim(SAMP[[1]])[1])
#low=apply(exp(SAMP[[1]])[1:(Iter-NBURNIN-1),],2,function(x)quantile(x,.025))
#up=apply(exp(SAMP[[1]])[1:(Iter-NBURNIN-1),],2,function(x)quantile(x,.975))
#results<-cbind(grid,c(low[1],low),c(med[1],med),c(up[1],up))

#Bayesian Summary of log(N(t))
ini<-1
med=apply(SAMP[[1]][ini:(Iter-NBURNIN-1),],2,median); #se=apply(exp(SAMP[[1]]),2,sd)/sqrt(dim(SAMP[[1]])[1])
low=apply(SAMP[[1]][ini:(Iter-NBURNIN-1),],2,function(x)quantile(x,.025))
up=apply(SAMP[[1]][ini:(Iter-NBURNIN-1),],2,function(x)quantile(x,.975))

effs1<-rep(0,Ngrid)
for (j in 1:Ngrid){
  print(effs1[j]<-effectiveSize(SAMP[[1]][1:(Iter-NBURNIN-1),j]))
}
results<-cbind(grid/scaling,c(low[1]-log(scaling),low-log(scaling)),c(med[1]-log(scaling),med-log(scaling)),c(up[1],up)-log(scaling),c(effs1[1],effs1))


##Plot results

plot(results[,1],results[,3],type="l",xlim=c(1,0),ylim=c(-3,3),ylab="log N(t)",xlab="No generations",col="white")
plot.res(results)

##True trajectory
x<-sort(c(0.299999,0.3,0.49999,0.5,seq(0,4,length.out=100)))
y<-x
y[x<.3]<-log(.5)
y[x>=.3 & x<.5]<-log(.1/2)
y[x>=.5]<-log(.5)
points(x,y,lty=2,type="l",lwd=1.5)

  