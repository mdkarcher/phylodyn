##Code for Tajima-based inference from a single locus

sufficient.stats<-function(groups,n){
  ##Dan Gusfield Algorithm 1.1 and 1.2 
  l<-length(groups$mut.groups)
  O<-matrix(0,nrow=n,ncol=l)
  frequency<-rep(1,n)
  haplotypes<-paste0(rep(0,l),collapse="")
  Index<-matrix(0,nrow=n,ncol=l)
  Lmat<-matrix(0,nrow=n,ncol=l)
  for (j in 1:l){
    O[,j]<-as.numeric(strsplit(groups$mut.groups[l-j+1],NULL)[[1]])
  }
  Index<- O%*%diag(seq(1,l))
  #(3)For each cell
  for (j in 1:n){
    haplotypes<-c(haplotypes,paste0(O[j,],collapse=""))
    for (i in 2:l){
      if (O[j,i]==1){Lmat[j,i]<-max(Index[j,1:(i-1)])}
    }
  }
  #correct for multiple haplotypes
  sort.hap<-sort(haplotypes[-1],index.return=T)
  base<-sort.hap$x[1]
  base.ind<-sort.hap$ix[1]
  remove<-0
  for (j in 2:n){
    if (base==sort.hap$x[j]){
      remove<-c(remove,sort.hap$ix[j])
      frequency[base.ind]<-frequency[base.ind]+1
      frequency[sort.hap$ix[j]]<-0
      
    }else{
      base<-sort.hap$x[j]
      base.ind<-sort.hap$ix[j]
    }
  }
  leaves<-apply(Index,1,max)
  card<-rev(groups$cardinality)
  carriers<-rev(groups$carriers) #I don't use this
  mylist<-list(list(d=0,x=0,y=0,parent=0))
  L<-apply(Lmat,2,max) #this vector has the nesting information, it has parents nodes
  parents<-sort(unique(L),d=T)
  i<-2
  for (j in parents){
    offspring<-seq(1,l)[L==j]
    offspringsize<-0
    for (no in offspring){
      offspringsize<-c(offspringsize,sum(leaves==no))
    }
    offspringsize<-offspringsize[-1]
    #offspringsize<-frequency[offspring]
    offspringsizelist<-unique(offspringsize)
    for (k in offspringsizelist){
      mylist[[i]]<-list(d=k,x=sum(offspringsize==k),y=card[offspring[offspringsize==k]],parent=j)  
      for (le in offspring[offspringsize==k]){
        leaves[leaves==le]<-j
      }
      i<-i+1
    }
    
  }
  mylist[[1]]<-NULL
  return(mylist)
}

group.data<-function(sort.main,n){
  mut.groups<-unique(sort.main$x)
  new.label<-seq(1,length(mut.groups))
  cardinality<-rep(0,length(mut.groups))
  carriers<-rep(0,length(mut.groups))
  for (j in 1:length(mut.groups)){
    cardinality[j]<-sum(sort.main$x==mut.groups[j])
    carriers[j]<-sum(as.numeric(strsplit(mut.groups[j],NULL)[[1]]))
  }
  if (sum(cardinality)!=max(sort.main$ix)){print("Error"); break}
  if (max(carriers)>n){print("Error");break}
  return(list(carriers=carriers,cardinality=cardinality,mut.groups=mut.groups))
}

bring.branch.lengths<-function(u,F){
  #u is the vector of intercoalescent times
  #given u and F, returns branch lengths ln+1,ln,ln-1,...,l3 and family sizes
  dimf<-nrow(F)
  diffM<-F[2:(dimf),2:dimf]-F[2:dimf,1:(dimf-1)]
  d<-rep(u[1],dimf) #for corresponding ln+1,ln,ln-1,...,l3
  firstzero<-rep(0,dimf-2) #for corresponding Sn,Sn-1,..S2
  familysize<-rep(2,dimf)
  coal_times<-cumsum(u)
  for (j in 1:(dimf-2)){
    firstzero[j]<-min(seq(j,(dimf-1))[diffM[j:(dimf-1),j]==0])
    d[j+1]<-coal_times[firstzero[j]]-coal_times[j]
  }
  
  d[dimf]<-u[dimf]
  firstzerotrans<-dimf+2-firstzero
  for (i in 1:(dimf-2)){
    count<-sum(firstzerotrans==(dimf-i+1))
    if (count==0){
      familysize[i+1]<-2
    }
    if (count==1){
      familysize[i+1]<-familysize[seq(1,dimf-2)[firstzerotrans==(dimf-i+1)]]+1
    }
    if (count==2){
      familysize[i+1]<-familysize[min(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)])]+familysize[max(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)])]
    }
    
  }
  familysize[length(familysize)]<-dimf+1
  return(list(d=d,familysize=familysize))
}
