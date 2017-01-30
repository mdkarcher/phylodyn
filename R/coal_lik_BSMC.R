#### Sequentially Markov Coalescent likelihood functions (SMC') ####
# Compute log likelihood under the SMC' model and other auxiliary funcitons


find.children.length <- function(tree, tdel, cor = 1) {
  #find the branch length of the branches that coalesced at time tdel
  #look at tdel found at tree2 but search for it in tree1
  edges<-tree$edge
  x<-tree$Nnode+1
  values0<-max(ape::node.depth.edgelength(tree))-ape::node.depth.edgelength(tree)
  values0<-values0*cor
  values1<-values0[(x+1):length(values0)]
  values<-sort(values1,decreasing=T)
  correct.label<-seq(x+1,2*x-1)
  ##do the crossreference
  reference<-matrix(c(values1,seq(x+1,2*x-1),rep(0,length(values1))),ncol=3)
  
  for (j in 1:length(values1)){
    reference[reference[,1]==values[j],3]<-correct.label[j]
  }
  y<-min(reference[reference[,1]>=tdel,1])
  findtdel<-seq(1:nrow(reference))[reference[,1]==y]
  #then
  children<-edges[edges[,1]==reference[findtdel,2],2]
  dist.to.children<-values0[children]
  
  return(dist.to.children)
}




#' Generates a matrix with ranked tree shape information
#' 
#' @param tree a phylo object
#'   
#' @return An F matrix containing ranked tree shape information
#' @export
#' 
#' @examples
#' create_F(ape::rcoal(4))
create_F<-function(tree){
  
  #n is the number of individual samples 
  edges<-tree$edge
  x<-tree$Nnode+1
  n_samples<-x
  values1<-max(ape::node.depth.edgelength(tree))-ape::node.depth.edgelength(tree)
  values1<-values1[(x+1):length(values1)]
  values<-sort(values1,decreasing=T)
  correct.label<-seq(x+1,2*x-1)
  ##do the crossreference
  reference<-matrix(c(values1,seq(x+1,2*x-1),rep(0,length(values1))),ncol=3)
  
  for (j in 1:length(values1)){
    reference[reference[,1]==values[j],3]<-correct.label[j]
  }
  newedges<-matrix(0,nrow=nrow(edges),ncol=4)
  newedges[,1:2]<-edges
  newedges[,3:4]<-edges
  
  for (j in (x+1):(2*x-1)){
    newedges[newedges[,1]==j,3]<-reference[reference[,2]==j,3]
    newedges[newedges[,2]==j,4]<-reference[reference[,2]==j,3]
  }
  edges<-newedges[,3:4]
  #tree.edge is edges correcting the labels
  F<-matrix(0,nrow=n_samples-1,ncol=n_samples-1)
  diag(F)<-n_samples:2
  F[cbind(2:(n_samples-1),1:(n_samples-2))]<-diag(F)[-(n_samples-1)]-2
  x<-nrow(edges)
  y<-max(edges)
  
  
  for (cini in 1:(n_samples-2)){
    ini<-cini+1
    for (j in y:(n_samples+2) ){
      F[ini,cini]<-F[(ini-1),cini]-sum(edges[edges[,1]==j,2]<=n_samples)
      ini<-ini+1
    }
    #remove those two rows
    find<-seq(1,nrow(edges))[edges[,1]==y]
    edges<-edges[-find,]
    edges[edges[,2]==y,2]<-1
    y<-y-1
  }
  return(F)
}

#I don't think this function is used anywhere
# get_F_values<-function(latent,init,Fl){
#   ntrees<-nrow(init$C)
#   latentloc<-seq(1:(ntrees-1))[latent[-ntrees]==1]
#   m.1<-sum(init$B[latentloc[1],])
#   j<-latentloc[1]
#   Fl.out<-matrix(0,nrow=m.1-1,ncol=m.1)
#   for (i in 1:(m.1-1)){
#     myF.in<-C2[(j-1),i]
#     myF.out<-n-C2[(j-1),((i+1):m.1)]+1
#     Fl.out[i,]<-c(rep(0,i),init$Fl[[j]][,(n-myF.in+1)][myF.out]/init$C2[(j-1),(i+1):m.1])
#   }
#   myF.list<-list(Fl.out[,-1])
#   
#   for (j in latentloc[2:length(latentloc)]){
#     m.1<-sum(init$B[j,])
#     Fl.out<-matrix(0,nrow=m.1-1,ncol=m.1)
#     for (i in 1:(m.1-1)){
#       myF.in<-C2[(j-1),i]
#       myF.out<-n-C2[(j-1),((i+1):m.1)]+1
#       Fl.out[i,]<-c(rep(0,i),init$Fl[[j]][,(n-myF.in+1)][myF.out]/init$C2[(j-1),(i+1):m.1])
#     }
#     myF.list<-c(myF.list,list(Fl.out[,-1]))
#   }
#   return(myF.list)
# }




coal_loglik_smc = function(init, f, grad = FALSE)
  #Revised: Jan 27,2017
  #TODO: Think about changing the ocurrances at the exact time and instead at the defined interval
{
  
nsize<-nrow(init$Fl[[1]])+1  
  if (init$ng != length(f))
    stop(paste("Incorrect length for f; should be", init$ng))
  
  #f.ext<-matrix(0,nrow=nrow(init$C),ncol=ncol(init$C))
  ntrees<-nrow(init$C)
  ExpZ<-matrix(0,nrow=ntrees,ncol=ncol(init$C))
  ExpDelta<-matrix(0,nrow=ntrees,ncol=ncol(init$C))
  loglik<-0
  f1<-rep(f,init$Nrep[1,])
  rep_idx = cumsum(init$Nrep[1,])
  rep_idx = cbind(rep_idx-init$Nrep[1,]+1,rep_idx)
  
  llnocoal = init$Delta[1,] * init$C[1,] * exp(-f1)
  lls<--init$Z[1,]*f1-llnocoal
  dll = apply(rep_idx,1,function(idx)sum(-init$Z[1,idx[1]:idx[2]]+llnocoal[idx[1]:idx[2]])) # gradient of log-likelihood wrt f_midpts
  
  loglik<-sum(lls[!is.nan(lls)])
  for (j in 2:(ntrees)){
      
    f.ext<-rep(f,init$Nrep[j,])
    rep_idx = cumsum(init$Nrep[j,])
    rep_idx = cbind(rep_idx-init$Nrep[j,]+1,rep_idx)
    
    if (init$latent[j]==1){
      m.1<-sum(init$C[j,]>1)
      logq<--init$Delta[j,]*init$C[j,]*exp(-f.ext)
      q<-exp(logq)
      Q<-(exp(f.ext)/init$C[j,])*(1-q)
      N1<-init$Delta[j,]*init$I[j,]
      N1[1]<-0
      N2<-(1-(1-logq)*q)*((exp(f.ext)/init$C[j,])^{2})*init$I[j,]
      N3<-((exp(f.ext)/init$C[j,])^{2})*(1/init$C[j,])*(-logq-2+(2-logq)*q)
      Pr.N1<-rep(0,ncol(init$C))
      Pr.N2<-rep(0,ncol(init$C))
      Pr.N4<-rep(0,ncol(init$C))
      
      Pii<-(init$Delta[j,]-Q)/init$C[j,]
      Pii[Pii<0]<-0
      ExpZ[j,1:m.1]<-init$C[j,1:m.1]*Pii[1:m.1]*init$I[j,1:m.1]
      lik1<-sum(ExpZ[j,])
      for (i in 1:(m.1-1)){
        myF.in<-init$C[j,i]
        myF.out<-nsize-init$C[j,((i+1):m.1)]+1
        myF<-init$Fl[[j]][,(nsize-myF.in+1)][myF.out]
        if (i==(m.1-1)){
          lik1<-lik1+Q[i]*(sum((1-q[(i+1):m.1])*myF))
          ExpZ[j,(i+1):m.1]<-ExpZ[j,(i+1):m.1]+Q[i]*(1-q[(i+1):m.1])*myF
          Pr.N2[i]<-sum((1-q[(i+1):m.1])*myF)
          Pr.N4[(i+1):m.1]<-Pr.N4[(i+1):m.1]+Q[i]*myF
        }else{
          lik1<-lik1+Q[i]*(sum((1-q[(i+1):m.1])*myF*c(1,exp(cumsum(logq[(i+1):(m.1-1)])))))
          ExpZ[j,(i+1):m.1]<-ExpZ[j,(i+1):m.1]+Q[i]*(1-q[(i+1):m.1])*myF*exp(cumsum(c(0,logq[(i+1):(m.1-1)])))
          Pr.N2[i]<-sum((1-q[(i+1):m.1])*myF*exp(cumsum(c(0,logq[(i+1):(m.1-1)]))))
          Pr.N4[(i+1):m.1]<-Pr.N4[(i+1):m.1]+Q[i]*myF*exp(cumsum(c(0,logq[(i+1):(m.1-1)])))
          if (i<(m.1-1)){
            Pr.N1[(i+1):(m.1-1)]<-Pr.N1[(i+1):(m.1-1)]+Q[i]*rev(cumsum(rev((1-q[(i+2):m.1])*myF[-1]*exp(cumsum(logq[(i+1):(m.1-1)]))))) #
            
          }
          
        }
      }
      ExpZ[j,]<-ExpZ[j,]/lik1
      ExpDelta[j,]<-(N1*Pr.N1+N2*Pr.N2+N3*init$C[j,]*init$I[j,]+N2*Pr.N4)/lik1
      
      if (!is.na(lik1)){
        loglik<-loglik+log(lik1)}
      
      
      
    } else { #nonlatent
      m.1<-sum(init$I[j,]>0)
      m.0<-min(seq(1,length(init$I[j,]))[init$I[j,]>0])
      indic<-seq(1,ncol(init$I))[init$I[j,]==0 & init$B[j,]==1]
      logq<--init$Delta[j,]*init$C[j,]*exp(-f.ext)
      q<-exp(logq)
      Q<-(exp(f.ext)/init$C[j,])*(1-q)
      N2<-(1-(1-logq)*q)*((exp(f.ext)/init$C[j,])^{2})*init$I[j,]
      if (length(indic)>0){
        if (min(indic)>m.0){ExpDelta[j,indic]<-init$Delta[j,indic]}
      }
      Lik1<-sum(init$I[j,]*Q*c(rev(exp(cumsum(rev(init$B[j,-1]*logq[-1])))),1))
      ExpZ[j,]<-init$Z[j,]
      if (m.1>1){
        ExpDelta[j,m.0:(m.0+m.1-1)]<-N2[m.0:(m.0+m.1-1)]*c(rev(exp(cumsum(rev(init$B[j,-1]*logq[-1])))),1)[m.0:(m.0+m.1-1)]/Lik1+
          init$Delta[j,m.0:(m.0+m.1-1)]*c(0,cumsum((init$I[j,]*Q*c(rev(exp(cumsum(rev(init$B[j,-1]*logq[-1])))),1)))[m.0:(m.0+m.1-2)])/Lik1
      }
      if (m.1==1){
        ExpDelta[j,m.0]<-N2[m.0]/(init$I[j,m.0]*Q[m.0])
      }
      loglik<-loglik+log(Lik1)-sum(init$Z[j,]*f.ext)
      
    }
    llnocoal = ExpDelta[j,] * init$C[j,] * exp(-f.ext)
    dll = dll+apply(rep_idx,1,function(idx)sum(-ExpZ[j,idx[1]:idx[2]]+llnocoal[idx[1]:idx[2]])) # gradient of log-likelihood wrt f_midpts
    
  }
  
  if (grad){
    
    return(dll)
    
  } else{
    return(loglik)}
}




U_split_smc = function(theta, init, invC, alpha, beta, grad=F)
{
  D=length(theta)
  f=theta[-D];tau=theta[D]
  invCf=invC%*%f
  if(!grad){
    loglik = coal_loglik_smc(init, f)
    logpri = ((D-1)/2+alpha)*tau - (t(f)%*%invCf/2+beta)*exp(tau)
    return(list(loglik = -loglik, logpri = -logpri, logpos = -(loglik+logpri)))
  }
  else{
    dU_res = -c(coal_loglik_smc(init, f, grad),((D-1)/2+alpha)-beta*exp(tau))
    return(dU_res)
    
  }
}


#Duplicated. It is coal_lik.R
#precBM = function(times, delta=1e-6)
#{
#  D=length(times)
#  diff1<-diff(times); diff1[diff1==0]<-delta;
#  diff<-1/diff1
#  Q<-spam(0,D,D)
#  if (D>2) Q[cbind(1:D,1:D)]<-c(diff[1]+ifelse(times[1]==0,1/delta,1/times[1]),diff[1:(D-2)]+diff[2:(D-1)],diff[D-1])
#  else Q[cbind(1:D,1:D)]<-c(diff[1]+ifelse(times[1]==0,1/delta,1/times[1]),diff[D-1])

# Q[cbind(1:(D-1),2:D)]=-diff[1:(D-1)]; Q[cbind(2:D,1:(D-1))]=-diff[1:(D-1)]
# return(Q)
#}


#'Extracts coalescent times from a Multi-Phylo object
#' @param MyTree \code{multiPhylo} object containing coalescent trees
#' @param sim integer specifying number of trees to read
#' @param factor numeric value to scale all the times
#' 
#' @return A matrix of sim rows. Entry x_{i,j} has the n-j+1-th coalescent time of the i-th tree
#' @export
#' 
#' @examples 
#' read_times(ape::rmtree(3,5),3,1)
read_times<-function(MyTree,sim,factor){
  ##This function reads a multi-phylo object of local genealogies and returns a matrix with coalescent times
  
  n<-MyTree[[1]]$Nnode+1
  if (n>2){
    D<-matrix(nrow=sim,ncol=n-1)
    D[1,]<-cumsum(ape::coalescent.intervals(MyTree[[1]])$interval.length)*factor
    if (max(ape::node.depth.edgelength(MyTree[[1]]))>ape::coalescent.intervals(MyTree[[1]])$total.depth) {
      D[1,]<-D[1,]+factor*(max(ape::node.depth.edgelength(MyTree[[1]]))-ape::coalescent.intervals(MyTree[[1]])$total.depth)
    }
    
    for (j in 1:sim){
      D[j,]<-factor*cumsum(ape::coalescent.intervals(MyTree[[j]])$interval.length)
      
      if (max(ape::node.depth.edgelength(MyTree[[j]]))>ape::coalescent.intervals(MyTree[[j]])$total.depth) {
        D[j,]<-D[j,]+factor*(max(ape::node.depth.edgelength(MyTree[[j]]))-ape::coalescent.intervals(MyTree[[j]])$total.depth)
      }
    }
    
    return(D)  }
  else{
    D<-rep(0,sim)
    for (j in 1:sim){
      D[j]<-factor*max(ape::node.depth.edgelength(MyTree[[j]]))*factor
    }
    return(D)
    
  }
}


#' Extracts all sufficient statistics for inferring Ne from sequential local genealogies
#' 
#' @param MyTree \code{multiPhylo} object containing coalescent trees
#' @param D a matrix of coalescent times
#' @param sim integer specifying number of trees to read
#' @param tol tolerance for detecting difference between two numeric values
#' @param cor (internal)
#' 
#' @return A list of sufficient statistics
#' @export
#' 
find_info2<-function(MyTree,D,sim,tol,cor=1){
  
  n<-MyTree[[1]]$Nnode+1
  suff<-find_sufficient(D,sim,tol*cor)
  sum(suff$latent)/sim ## % of invisible transitions
  latent<-suff$latent
  t_del<-suff$t_del
  t_new<-suff$t_new
  where_new<-suff$where_hap_new
  where_del<-suff$where_hap_del
  
  
  Fl<-list(create_F(MyTree[[1]]))
  for (j in 2:sim){
    Fl<-c(Fl,list(create_F(MyTree[[j]])))
  }
  
  #Uses F-matrices and coalescent times
  info_times<-matrix(0,nrow=sim,ncol=6)
  #number of branches from-to
  for (j in 2:sim){
    if (latent[j]==0){
      if (where_del[j]<(n-1)){
        children<-Fl[[(j-1)]][where_del[j],1:(min(where_new[j],where_del[j]))]-Fl[[(j-1)]][(where_del[j]+1),1:(min(where_new[j],where_del[j]))]
        
      }else{
        children<-Fl[[(j-1)]][where_del[j],1:(min(where_new[j],where_del[j]))]
      }
      if (where_new[j]<(n-1)){
        children_new<-Fl[[(j)]][where_new[j],1:(min(where_new[j],where_del[j]))]-Fl[[(j)]][(where_new[j]+1),1:(min(where_new[j],where_del[j]))]
        
      }else{
        children_new<-Fl[[(j)]][where_new[j],1:(min(where_new[j],where_del[j]))]
        
      }
      if (min(children)==2){
        #cherry
        info_times[j,1:3]<-c(2,0,min(t_del[j],t_new[j]))
      }else{
        if (min(children)==max(children) & max(children)==1){
          info_times[j,1:3]<-c(1,0,min(t_del[j],t_new[j]))
          
        }else{
          #if (where_del[j]>=where_new[j]){
          val0<-0
          val0_new<-0
          if (min(children)==0){val0<-max(seq(1,length(children),by=1)[children==0])}
          val<-max(seq(1,length(children),by=1)[children==1])
          if (min(children_new)==0){val0_new<-max(seq(1,length(children_new),by=1)[children_new==0])}
          if (min(children_new)<=1){val_new<-max(seq(1,length(children_new),by=1)[children_new==1])}else{val_new=0}
          if (val0>0 & val0_new>0){
            if (  ((D[j-1,val0]==D[j,val0_new]) & (D[j-1,val]==D[j,val_new])) || (sum(children-children_new)==0) & (val_new<where_new[j]) ){
              info_times[j,]<-c(1,D[j,val0],D[j,val],2,D[j,val],min(t_del[j],t_new[j]))
              #children are the same even when topology is not the same (2)
            }else{
              if (val0_new<val0) {info_times[j,1:3]<-c(1,D[j,val_new],min(t_del[j],t_new[j]))}
              if (val0_new==val0){info_times[j,1:3]<-c(1,D[j,val0_new],min(t_del[j],t_new[j]))}
              if (val0_new>val0){info_times[j,1:3]<-c(1,D[j,val],min(t_del[j],t_new[j]))}
            }
          }else{if (val0_new<val0) {info_times[j,1:3]<-c(1,D[j,val_new],min(t_del[j],t_new[j]))}
            if (val0_new==val0 & val0_new>0){info_times[j,1:3]<-c(1,D[j,val0_new],min(t_del[j],t_new[j]))}
            if (val0_new==val0 & val0_new==0 & val_new==val & val<where_new[j]){info_times[j,]<-c(1,0,D[j,val],2,D[j,val],min(t_del[j],t_new[j]))}
            if (val0_new>val0 & val0_new==val){info_times[j,1:3]<-c(1,D[j,val],min(t_del[j],t_new[j]))}
            if (val0_new>val0 & val_new==val){info_times[j,1:3]<-c(1,D[j,val],min(t_del[j],t_new[j]))}
            if (val0_new==val0 & val0_new==0 & val_new!=val){info_times[j,1:3]<-c(1,0,min(t_del[j],t_new[j]))}
            
            
          }
          
        }
      }}}
  return(list(info_times=info_times,Fl=Fl,latent=latent,t_new=t_new,t_del=t_del,D=D,sim=sim,n=n))
}

#This function needs to be improved

find_sufficient <- function(D, sim, tol) {
  n<-ncol(D)+1
  ##This function returns three vectors of statistics needed for likelihood calculations
  if (n>2){
    t_del<-rep(0,sim)
    t_new<-rep(0,sim)
    where_hap_new<-rep(0,sim)
    where_hap_del<-rep(0,sim)
    for (j in 2:sim){
      diff<-D[j,]-D[(j-1),]
      if (sum(diff>tol)==1){
        val1<-min(seq(1,n-1)[diff>tol])
        val2<-max(seq(1,n-1)[diff>tol])
        t_del[j]<-D[j-1,val1]
        t_new[j]<-D[j,val2]
        where_hap_new[j]<-val2
        where_hap_del[j]<-val1
      }else {
        if (sum(-diff>tol)==1){
          val1<-min(seq(1,n-1)[-diff>tol])
          t_del[j]<-D[j-1,val1]
          t_new[j]<-D[j,val1]
          where_hap_new[j]<-val1
          where_hap_del[j]<-val1
        } else{
          
          if (sum(abs(diff)>tol)>1){
            val1<-min(seq(1,n-1)[abs(diff)>tol])
            if (diff[val1]>0) {
              t_del[j]<-D[j-1,val1]
              val2<-max(seq(1,n-1)[diff>tol])
              t_new[j]<-D[j,val2]
              where_hap_new[j]<-val2
              where_hap_del[j]<-val1} else{
                if (diff[val1]<0) { ##OK
                  val2<-max(seq(1,n-1)[abs(diff)>tol])
                  t_new[j]<-D[j,val1]
                  t_del[j]<-D[(j-1),val2]
                  where_hap_del[j]<-val2
                  where_hap_new[j]<-val1
                }
              }
          }
        }
      }
    }
    
    latent<-rep(0,sim)
    latent[t_del==t_new]<-1
    latent[1]<-0
    return(list(latent=latent,t_new=t_new,t_del=t_del,where_hap_new=where_hap_new,where_hap_del=where_hap_del))
  } else{
    t_del<-rep(0,sim)
    t_new<-rep(0,sim)
    diff<-diff(D)
    for (j in 1:(sim-1)){
      if (sum(diff[j]>tol)==1){
        t_del[j]<-D[j]
        t_new[j]<-D[j+1]
      }else {
        if (sum(-diff[j]>tol)==1){
          t_del[j]<-D[j]
          t_new[j]<-D[j+1]
        } 
      }
    }
    
    latent<-rep(0,sim)
    latent[t_del==t_new]<-1
    latent[1]<-0
    return(list(latent=latent,t_new=t_new,t_del=t_del))
    
  }
  
}


#Declared a little bit different in Phyloinfer.R
#splitHMC = function (current.q, U, rtEV, EVC, eps=.1, L=5,current.U,current.grad){
#Function from Lan S, Palacios JA, Karcher, M, Minin VN, Babak Shahbaba. An efficient Bayesian inference framework for coalescent-based nonparametric phylodynamics. 2015
#   browser()

# initialization
#D=length(current.q)
#   q=current.q

# sample momentum
#   p <- rnorm(D)

# calculate current energy

#   current.E <- current.U + sum(p^2)/2

# current.E <- U(q) + sum(p^2)/2


#   randL = ceiling(runif(1)*L)
#   p = p - eps/2*current.grad #corrected Feb 22,2015 p-(eps/2)*current.grad

# p = p - (eps/2)*U(q,T)/sqrt(sum(U(q,T)^2))
#   qT = rtEV*(t(EVC)%*%q[-D]); pT = t(EVC)%*%p[-D]
#   A=t(qT)%*%qT;
# Alternate full steps for position and momentum
#   for (l in 1:randL)
#   {


# Make a half step for the initial half dynamics
#	  A=t(qT)%*%qT;
#	  C1=p[D]^2+A*exp(q[D]); C2=2*atanh(p[D]/sqrt(C1))
# 	  p[D] = sqrt(C1)*tanh((-sqrt(C1)*eps/2+C2)/2)
#	  q[D] = log((C1-p[D]^2)/A)

#       p[D] <- p[D] - eps/2*A/2*exp(q[D])
#       q[D] <- q[D] + eps/2*p[D]

# Make a full step for the middle dynamics
#       Cpx = complex(mod=1,arg=-rtEV*exp(q[D]/2)*eps)*complex(re=qT*exp(q[D]/2),im=pT)
#       qT = Re(Cpx)*exp(-q[D]/2); pT = Im(Cpx)
#       q[-D] = EVC%*%(qT/rtEV)

# Make a half step for the last half dynamics
#       A=t(qT)%*%qT;
#	  C1=p[D]^2+A*exp(q[D]); C2=2*atanh(p[D]/sqrt(C1))
#	  p[D] = sqrt(C1)*tanh((-sqrt(C1)*eps/2+C2)/2)
#	  q[D] = log((C1-p[D]^2)/A)

#       q[D] <- q[D] + eps/2*p[D]
#       p[D] <- p[D] - eps/2*A/2*exp(q[D])

#       g = U(q,T)
#  g<-g/sqrt(sum(g^2))
#       if(l!=randL){
#           pT = pT - eps*(t(EVC)%*%g[-D]); p[D] = p[D] - eps*g[D]
#       }
#   }
#   p[-D] = EVC%*%pT - eps/2*g[-D]; p[D] = p[D] - eps/2*g[D]



# Evaluate potential and kinetic energies at start and end of trajectory
#   new.u<-U(q)
#   proposed.E <- new.u + sum(p^2)/2

# Accept or reject the state at end of trajectory, returning either
# the position at the end of the trajectory or the initial position
#   logAP = -proposed.E + current.E

#   if( is.finite(logAP)&(log(runif(1))<min(0,logAP)) ) return (list(q = q, Ind = 1,current.u=new.u,current.grad=g))
#   else return (list(q = current.q, Ind = 0,current.u=current.u,current.grad=current.grad))

#}


get_data<-function(grid,sim,D,n,info_times,Fl,latent,t_new,t_del,tol){
  
  grid.extended.list<-list(sort(unique(c(grid,as.vector(D[1,])))) )
  check<-rep(0,nrow(D))
  check[1]<-length(grid.extended.list[[1]])
  for (j in 2:(nrow(D))){
    myy<-sort(unique(c(grid,as.vector(D[j,]))))
    #  myy<-myy[myy<=maxt] ##Truncation
    grid.extended.list<-c(grid.extended.list,list(myy))
    check[j]<-length(grid.extended.list[[j]])
  }
  
  
  Z<-matrix(0,nrow=sim,ncol=max(check)-1)
  C<-matrix(1,nrow=sim,ncol=max(check)-1)
  Delta<-matrix(nrow=sim,ncol=max(check)-1)
  Nrep<-matrix(0,nrow=sim,ncol=length(grid)-1)
  I<-matrix(0,nrow=sim,ncol=max(check)-1)
  B<-matrix(0,nrow=sim,ncol=max(check)-1)
  o<-coal_lik_init(grid.extended.list[[1]], c(n,rep(0,length(grid.extended.list[[1]])-1)), D[1,], grid)
  
  Nrep[1,]<-o$gridrep
  Z[1,1:(check[1]-1)]<-o$y
  C[1,1:(check[1]-1)]<-o$C ##the quadratic factor for the first one
  Delta[1,1:(check[1]-1)]<-o$D
  
  
  for (j in 2:sim){
    
    prevo<-coal_lik_init(grid.extended.list[[(j-1)]],c(n,rep(0,length(grid.extended.list[[(j-1)]])-1)),D[(j-1),],grid)
    Nrep[j,]<-prevo$gridrep
    C[j,]<-prevo$l
    Delta[j,1:(check[j]-1)]<-prevo$D
    if (latent[j]==1){
      m.1<-sum(C[j,]>1)
      I[j,]<-c(rep(1,m.1),rep(0,ncol(I)-m.1))
      if (sum(I[j,])==0) {break}
    }else{
      
      find01<-min(seq(1:length(Delta[j,]))[cumsum(Delta[j,])>info_times[j,2]])
      if (Delta[j,1]>info_times[j,3]) {
        find2<-1
        I[j,1]<-info_times[j,1] } else{
          find2<-min(seq(1:length(Delta[j,]))[cumsum(Delta[j,])>=info_times[j,3]])
          I[j,find01:find2]<-info_times[j,1]
        }
      if (info_times[j,4]==2){
        find1<-max(seq(1:length(Delta[j,]))[cumsum(Delta[j,])<=(info_times[j,5]+tol)])+1
        find2<-max(seq(1:length(Delta[j,]))[cumsum(Delta[j,])<=info_times[j,6]])
        if (find1<find2){I[j,find1:find2]<-2}
      }
      find<-min(seq(1:length(Delta[j,]))[cumsum(Delta[j,])>t_new[j]])
      find2<-min(seq(1:length(Delta[j,]))[cumsum(Delta[j,])>=t_del[j]])
      Z[j,find]<-1
      B[j,find01:find]<-1 ##from tp to tnew
      if (find>1){
        Delta[j,find]<-t_new[j]-sum(Delta[j,1:(find-1)]) #correction on Delta for t_new
      } else {
        Delta[j,find]<-t_new[j]
      }
      if (sum(I[j,])==0) {break}
    }
  }
  
  midpts<-grid[-length(grid)]+diff(grid)/2
  
  #Note that this function is a bit different than Q_matrix. Fix this code
  #There is a wraper
  #Q.matrix<-function(input,s.noise,signal){
  #   n2<-nrow(input)
  #   diff1<-diff(input)
  #   diff1[diff1==0]<-s.noise #correction for dividing over 0
  #   diff<-(1/(signal*diff1))
  #   Q<-spam::spam(0,n2,n2)
  #   if (n2>2){
  #       Q[cbind(seq(1,n2),seq(1,n2))]<-c(diff[1],diff[1:(n2-2)]+diff[2:(n2-1)],diff[n2-1])+(1/signal)*rep(s.noise,n2)} else {Q[cbind(seq(1,n2),seq(1,n2))]<-c(diff[1],diff[n2-1])+(1/signal)*rep(s.noise,n2)}
  #   Q[cbind(seq(1,n2-1),seq(2,n2))]<--diff[1:(n2-1)]
  #   Q[cbind(seq(2,n2),seq(1,n2-1))]<--diff[1:(n2-1)]
  #   return(Q)}
  
  
  
  invC<-Q_matrix(as.matrix(midpts),0,1)
  invC[1,1] <- invC[1,1]+.00001 #fudge to be able to compute the cholC
  
  #diag(invC)<-diag(invC)+.000001 #fudge to be able to compute the cholC
  #invC[1,1]<-invC[1,1]+.00001
  eig=spam::eigen.spam(invC,T)
  EV=eig$values; EVC=eig$vectors
  
  Cm=spam::solve.spam(invC)
  cholC<-chol(Cm)
  rtEV=sqrt(EV)
  lik_init =list(Z=Z,Delta=Delta,I=I,Nrep=Nrep,C=C,B=B,ng=(length(grid)-1),Fl=Fl,latent=latent)
  
  return(list(lik_init=lik_init,invC=invC,rtEV=rtEV,EVC=EVC))
  
  
}

#' Plot SMCP Results
#' 
#' @param results the output of smcp_sampling().
#'
#' @export
plot_res <- function(results) {
  temp.time = c(0,rep(results[,1], each = 2))
  temp.time<-temp.time[-(length(temp.time))]
  temp.popsize1 = rep(results[,2], each = 2)
  temp.popsize2 = rep(results[,4], each = 2)
  xx.alpha = c(temp.time, rev(temp.time))
  yy.alpha = c(temp.popsize1, rev(temp.popsize2))
  graphics::polygon(xx.alpha,yy.alpha,col=grDevices::adjustcolor("pink",alpha.f=0.5),border=NA)
  #points(results[,1],results[,3],col="#3182bd",type="l",lwd=2.5)
  graphics::points(results[,1],results[,3],col="red",type="l",lwd=1.75)
  
  graphics::points(results[,1],results[,2],col="red",type="l",lwd=0.5)
  graphics::points(results[,1],results[,4],col="red",type="l",lwd=0.5)
}
