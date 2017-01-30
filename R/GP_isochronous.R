#################
####### Metropolis-Hastings steps for Isochronous sampled data
#################


###Definition of Q matrix and Gaussian Processes

#library("spam")

##Compare time with defining Q first and then declaring it as sparse

##The browninan motion, efficient 
#Q_matrix<-function(input,s.noise,signal){
#	n2<-nrow(input)
#	diff<-(1/(signal*diff(input)))
#	Q<-spam(0,n2,n2)	
#	Q[cbind(seq(1,n2),seq(1,n2))]<-c(diff[1],diff[1:(n2-2)]+diff[2:(n2-1)],diff[n2-1])
##+(1/signal)*rep(s.noise,n2)
#	Q[cbind(seq(1,n2-1),seq(2,n2))]<--diff[1:(n2-1)]
#	Q[cbind(seq(2,n2),seq(1,n2-1))]<--diff[1:(n2-1)]
#	Q[1,1]<-(1/input[1]-Q[1,2])
#	return(Q)}


##The intrinsic one
# Q_matrix<-function(input,s.noise,signal){
# 	n2<-nrow(input)
# 	diff<-(1/(signal*diff(input)))
# 	Q<-spam(0,n2,n2)	
# 	Q[cbind(seq(1,n2),seq(1,n2))]<-c(diff[1],diff[1:(n2-2)]+diff[2:(n2-1)],diff[n2-1])+(1/signal)*rep(s.noise,n2)
# 	Q[cbind(seq(1,n2-1),seq(2,n2))]<--diff[1:(n2-1)]
# 	Q[cbind(seq(2,n2),seq(1,n2-1))]<--diff[1:(n2-1)]
# 	return(Q)}
#   

# Redundant: Use Q_matrix
# Qmatrix<-function(input,s.noise,signal)
# {
#   n2<-nrow(input)
#   diff1<-diff(input)
#   diff1[diff1==0]<-s.noise #correction for dividing over 0
#   diff<-(1/(signal*diff1))
#   Q <- spam::spam(0,n2,n2)	
#   if (n2>2)
#   {
#     Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1],diff[1:(n2-2)]+diff[2:(n2-1)],
#                                        diff[n2-1])+(1/signal)*rep(s.noise,n2)
#   } 
#   else
#   {
#     Q[cbind(seq(1,n2),seq(1,n2))] <- c(diff[1],diff[n2-1])+(1/signal)*rep(s.noise,n2)
#   }
#   
#   
#   Q[cbind(seq(1,n2-1),seq(2,n2))] <- -diff[1:(n2-1)]
#   Q[cbind(seq(2,n2),seq(1,n2-1))] <- -diff[1:(n2-1)]
#   return(Q)
# }
  

GP.prior<-function(signal,s.input,s.noise)
{
  ##We assume s.input is already ordered
	Q<-Q_matrix(as.matrix(s.input),s.noise,signal)
  L<-spam::chol.spam(Q)
	y2<-stats::rnorm(nrow(Q),0,1)
	g<-spam::backsolve(L,y2)
	return(g)
}


GP.posterior<-function(X,y,signal,s.tilde,s.noise)
{
	X2<-as.matrix(c(s.tilde,X))
	X2.2<-sort(X2,index.return=TRUE)
	n2<-length(X2)
	ind1<-rep(0,n2)
	n1<-length(s.tilde)
	ind1[X2.2$ix<=n1]<-1
	ind2<-c(ind1[2:n2],0)
	ind3<-c(0,ind1[1:(n2-1)])
	ind1<-ind1+ind2+ind3
	ind1[ind1>1]<-1
	
  #Q<-Q_matrix(as.matrix(X2.2$x),s.noise,signal)
	Q<-Q_matrix(as.matrix(X2.2$x[ind1==1]),s.noise,signal)
	induse<-X2.2$ix[ind1==1]
	
  #var.f<-Q[X2.2$ix<=n1,X2.2$ix<=n1]
	var.f<-Q[induse<=n1,induse<=n1]
	
	if (n1>1)
  {
		L.1 <- spam::chol.spam(var.f)
		part1 <- spam::forwardsolve.spam(L.1, -Q[induse<=n1, induse>n1] %*% y[induse[induse>n1] - n1],
		                                 transpose = TRUE, upper.tri = TRUE)
		mean.f2 <- spam::backsolve.spam(L.1, part1)
		y2 <- stats::rnorm(length(s.tilde), 0, 1)
		g2 <- mean.f2 + spam::backsolve.spam(L.1, y2)
		one <- induse[induse <= n1]
		return(list(mean_f = mean.f2[order(one)], g = g2[order(one)]))
	}
	if (n1==1)
	{
		part1<--Q[induse<=n1,induse>n1]%*%y[induse[induse>n1]-n1]
		mean.f2<-part1/var.f
		y2<-stats::rnorm(length(s.tilde),0,1)
		g2<-mean.f2+y2/sqrt(var.f)
		one<-induse[induse<=n1]
		return(list(mean_f=mean.f2[order(one)],g=g2[order(one)]))
	}
}


#GP.posterior<-function(X,y,signal,s.tilde,s.noise){
#	X2<-as.matrix(c(s.tilde,X))
#	X2.2<-sort(X2,index.return=TRUE)
#	Q<-Q_matrix(as.matrix(X2.2$x),s.noise,signal)
#	n1<-length(s.tilde)
#	var.f<-Q[X2.2$ix<=n1,X2.2$ix<=n1]
#	if (n1>1){
#		L.1<-chol.spam(var.f)
#		part1<-forwardsolve.spam(L.1,-Q[X2.2$ix<=n1,X2.2$ix>n1]%*%y[X2.2$ix[X2.2$ix>n1]-n1],transpose=T,upper.tri=T)
#		mean.f2<-backsolve.spam(L.1,part1)
#		y2<-stats::rnorm(length(s.tilde),0,1)
#		g2<-mean.f2+backsolve.spam(L.1,y2)
#		one<-X2.2$ix[X2.2$ix<=n1]
#		return(list(mean_f=mean.f2[order(one)],g=g2[order(one)]))
#	}
#	if (n1==1)
#	{
#		part1<--Q[X2.2$ix<=n1,X2.2$ix>n1]%*%y[X2.2$ix[X2.2$ix>n1]-n1]
#		mean.f2<-part1/var.f
#		y2<-stats::rnorm(length(s.tilde),0,1)
#		g2<-mean.f2+y2/sqrt(var.f)
#		one<-X2.2$ix[X2.2$ix<=n1]
#		return(list(mean_f=mean.f2[order(one)],g=g2[order(one)]))
#	}
#}


lambda.prior<-function(x,lambda.mean,alpha.1,alpha.2,alpha.3)
{
  z<-lambda.mean*alpha.2+(alpha.3/alpha.1)*exp(-alpha.1*lambda.mean)
  #indicator<-rep(1,length(x))
  #indicator[x>=lambda.mean]=rep(0,sum(x>=lambda.mean))
  if (x<=lambda.mean)
  {
    indicator=1
  }
  else
  {
    indicator=0
  }
  return((1/z)*((alpha.2*indicator)+alpha.3*exp(-alpha.1*x)*(1-indicator)))
}



GP.posterior2<-function(X,y,signal,s.tilde,s.noise)
{
	X2<-as.matrix(c(s.tilde,X))
	X2.2<-sort(X2,index.return=TRUE)
	Q<-Q_matrix(as.matrix(X2.2$x),s.noise,signal)
	n1<-length(s.tilde)
	n2<-nrow(X2)
	L2<-chol(Q[X2.2$ix>n1,X2.2$ix>n1])
	loglik<--.5*(t(y[X2.2$ix[X2.2$ix>n1]-n1])%*%Q[X2.2$ix>n1,X2.2$ix>n1]%*%y[X2.2$ix[X2.2$ix>n1]-n1])+sum(log(diag(L2)))-.5*(n2-n1)*log(2*pi)
	return(loglik)
}


########Metropolis-Hastings
###

sigmoidal<-function(x){(1+exp(-x))^-1}


number.thinned<-function(s,T,lambda,signal,info,s.noise,n,t,coal.factor,m)
{
	# m<-rep(0,n-1)
	# for (k in n:2){
	# 	m[n-k+1]<-sum(info[,2]==k)-1
	# }
  ind<-2*info[,3]-1
  loglik<-(-lambda)*sum(t*coal.factor)+sum(log(sigmoidal(info[,4]*ind)))+nrow(info)*log(lambda)+sum(m*log(coal.factor))

	add<-stats::rbinom(n-1,1,.5)
	new<-stats::runif(n-1,min=s[1:n-1],max=s[2:n])
	l<-length(s) #Just updated in January 2017. Needs to be checked
  g.new<-rep(0,l-1)
	g.new[add==1]<-GP.posterior(as.matrix(info[,1]),info[,4],signal,as.matrix(new[add==1]),s.noise)$g
	accept<-stats::runif(n-1)<add*lambda*t*sigmoidal(-g.new)*coal.factor/(m+1)
	info<-rbind(info,cbind(new[accept],seq(n,2)[accept],rep(0,sum(accept)),g.new[accept]))
  m<-m+accept
  #	info<-rbind(info,to.add)
  #	m<-rep(0,n-1)
  #	for (k in n:2){
  #		m[n-k+1]<-sum(info[,2]==k)-1
  #	}
  part1<-add==0 & m>0
	if (sum(part1)>0)
  {
    where<-ceiling(stats::runif(n-1)*m)[part1]
		co<-seq(n,2)[part1]

		if (m[part1][1]==1)
    {
#			who<-info[info[,2]==co[1] & info[,3]==0,]
			find<-seq(1,nrow(info))[info[,2]==co[1] & info[,3]==0]
		}
    else
    {
#			who<-info[info[,2]==co[1] & info[,3]==0,][where[1],]
			find<-seq(1,nrow(info))[info[,2]==co[1] & info[,3]==0][where[1]]
		}
		if (length(co)>=2)
    {
			for (j in 2:length(co))
      {
				if (m[part1][j]==1)
        {
#					who<-rbind(who,info[info[,2]==co[j] & info[,3]==0,])
					find<-rbind(find,seq(1:nrow(info))[info[,2]==co[j] & info[,3]==0])
				}
        else
        {
#					who<-rbind(who,info[info[,2]==co[j] & info[,3]==0,][where[j],])
					find<-rbind(find,seq(1:nrow(info))[info[,2]==co[j] & info[,3]==0][where[j]])
				}
			}
		}
		to.delete<-stats::runif(length(co))<m[part1]/(lambda*t[part1]*sigmoidal(-info[find,4])*coal.factor[part1])
		if (sum(to.delete)>0)
    {
			info<-info[-(find[to.delete]),]
      m[part1][to.delete]<-m[part1][to.delete]-1
		}
	}
	return(list(info=info,loglik=loglik,m=m))
}



location.thinned.uniform<-function(s,T,lambda,signal,N,info,s.noise)
{
	n<-sum(info[,3])
	m<-nrow(info)-sum(info[,3])
	coal<-info[,2][info[,3]==0]
	new<-stats::runif(m,min=s[N-coal+1],max=s[N-coal+2])
	g.new<-GP.posterior(as.matrix(info[,1]),info[,4],signal,as.matrix(new),s.noise)$g
	change<-stats::runif(length(coal))<sigmoidal(-g.new)/sigmoidal(-info[,4][info[,3]==0])
	
	info[,c(1,4)][info[,3]==0][change==TRUE]<-cbind(new[change==TRUE],g.new[change==TRUE])
	
#info[,1][info[,3]==0][change==TRUE]<-new[change==TRUE]
#	info[,4][info[,3]==0][change==TRUE]<-g.new[change==TRUE]
	return(info)
}




location.thinned.uniform2<-function(s,T,lambda,signal,N,info,s.noise,t,m)
{
  l<-length(s)
  m2<-m
  m2[m>0]<-1 #only where there are latent points
  n<-sum(info[,3])
  where<-seq(1:length(t))[stats::rmultinom(1,1,t*m2)==1] #sample interval proportional to time length
  who<-info[,1][info[,3]==0 & info[,1]<s[(where+1)] & info[,1]>=s[where]]
  ma<-length(who)
  #m<-nrow(info)-n
	#who<-info[,1][info[,3]==0]
	new<-stats::runif(ma,min=s[where],max=s[where+1])
#  for (j in 1:ma){
#   k<-sum(s<who[j])
#		new[j]<-stats::runif(1,min=s[k],max=s[k+1])
#	}
	g.new<-GP.posterior(as.matrix(info[,1]),info[,4],signal,as.matrix(new),s.noise)$g
	change<-stats::runif(ma)<sigmoidal(-g.new)/sigmoidal(-info[,4][info[,3]==0 & info[,1]<s[(where+1)] & info[,1]>=s[where]])
	info[,1][info[,3]==0 & info[,1]<s[(where+1)] & info[,1]>=s[where]][change==TRUE]<-new[change==TRUE]
	info[,4][info[,3]==0 & info[,1]<s[(where+1)] & info[,1]>=s[where]][change==TRUE]<-g.new[change==TRUE]
	return(info)
}




slice.sampling<-function(data,signal,s.noise)
{
	theta<-stats::runif(1,0,2*pi)
	Y<-as.matrix(data[,1])
	X2.2<-sort(Y,index.return=TRUE)
	Q<-Q_matrix(as.matrix(X2.2$x),s.noise,signal)
	cholQ<-spam::chol.spam(Q)
	v <- spam::backsolve.spam(cholQ,stats::rnorm(nrow(Y),0,1))
	
	v<-v[order(X2.2$ix)]
	v0<-data[,4]*sin(theta)+v*cos(theta)
	v1<-data[,4]*cos(theta)-v*sin(theta)
	theta.min<-0
	theta.max<-2*pi
	u<-stats::runif(1)
	indicator<-1-2*data[,3]
	loglik<-sum(log(1+exp(indicator*data[,4])))*-1
	logh<-log(u)+loglik
	keep<-1
	while (keep==1)
  {
		theta.prime<-stats::runif(1,theta.min,theta.max)
		proposal<-v0*sin(theta.prime)+v1*cos(theta.prime)
		loglik2<--sum(log(1+exp(indicator*proposal)))
		if (loglik2>logh)
    {
      keep<-2
    }
		else
    {
			if (theta.prime<theta)
      {
        theta.min<-theta.prime
      }
      else
      {
        theta.max<-theta.prime
      }
		}
	}
	data[,4]<-proposal
	return(list(data=data,Q=Q,order=X2.2$ix))
}


