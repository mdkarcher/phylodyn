##Code for Tajima-based inference from a single locus
group_sufficient<-function(sufficient){
  mylist<-list(list(x=0,y=0))
  isleaf<-rep(0,nrow(sufficient$nodes))
  i<-1
  for (j in sufficient$nodes[,1]){
    isleaf[i]<-sum(sum(sufficient$leaves_o==j)>0)
    i<-i+1
  }
  parents<-unique(sufficient$nodes[,2])
  for (j in parents){
    sizes<-unique(sufficient$nodes[sufficient$nodes[,2]==j,3])
    for (i in sizes){
      who<-seq(1:nrow(sufficient$nodes))[sufficient$nodes[,2]==2 & sufficient$nodes[,3]==i]
      if (length(who)==sum(isleaf[who])) {
        #remove
      }
    }
    
  }
}

sufficient_stats2<-function(data){
  ##This is a modified version that does not group nodes with the same size
  #data is a matrix of m rows (snps) and 1 column of concateted 0s and 1s. The number of characters in each row is n, the number of individuals
  sorted<-sort(data[,1],index.return=T)
  n<-nchar(data[1,])
  groups<-group_data(sorted,n)
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
  leaves_o<-leaves
  card<-rev(groups$cardinality)
  carriers<-rev(groups$carriers) #I don't use this
 # mylist<-list(list(x=0,y=0))
  orderlist<-0
  parentlist<-0
  famsize<-0
  muts<-0
  L<-apply(Lmat,2,max) #this vector has the nesting information, it has parents nodes
  parents <- sort(unique(L), decreasing = TRUE)
  i<-2
  for (j in parents){
    offspring<-seq(1,l)[L==j]
    offspringsize<-0
    for (no in offspring){
      # offspringsize<-c(offspringsize,sum(leaves==no))
      #mylist[[i]]<-list(x=1,y=card[no])
      muts<-c(muts,card[no])
      famsize<-c(famsize,sum(leaves==no))
      parentlist<-c(parentlist,j)
      orderlist<-c(orderlist,no)
      leaves[leaves==no]<-j
      i<-i+1
    }
  }
  mylist[[1]]<-NULL
  return(list(leaves_o=leaves_o,nodes=cbind(orderlist[-1],parentlist[-1],famsize[-1],muts[-1])))
}

sufficient_stats<-function(data,groups,n){
  #data is a matrix of m rows (snps) and 1 column of concateted 0s and 1s. The number of characters in each row is n, the number of individuals
  sorted<-sort(data[,1],index.return=T)
  n<-nchar(data[1,])
  groups<-group_data(sorted,n)
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
  leaf<-leaves
  card<-rev(groups$cardinality)
  carriers<-rev(groups$carriers) #I don't use this
  mylist<-list(list(x=0,y=0))
  orderlist<-0
  parentlist<-0
  famsize<-0
  L<-apply(Lmat,2,max) #this vector has the nesting information, it has parents nodes
  parents <- sort(unique(L), decreasing = TRUE)
  i<-2
  for (j in parents){
    offspring<-seq(1,l)[L==j]
    offspringsize<-0
    leafindicator<-0
    for (no in offspring){
      offspringsize<-c(offspringsize,sum(leaves==no))
      leafindicator<-c(leafindicator,sum(leaf==no))
    }
    offspringsize<-offspringsize[-1]
    leafindicator<-leafindicator[-1]
    #offspringsize<-frequency[offspring]
    offspringsizelist<-unique(offspringsize)
    for (k in offspringsizelist){
      dups<-sum(leafindicator[offspringsize==k]>=1)
      if (dups>0){
        mylist[[i]]<-list(x=sum(offspringsize==k),y=card[offspring[offspringsize==k & leafindicator[offspringsize==k]>=1]]) 
        famsize<-c(famsize,k)
        parentlist<-c(parentlist,j)
        orderlist<-c(orderlist,min(offspring[offspringsize==k & leafindicator[offspringsize==k]>=1]))
        for (le in offspring[offspringsize==k & leafindicator[offspringsize==k]==1]){
          leaves[leaves==le]<-j
        }
        i<-i+1
      }
      left<-sum(leafindicator[offspringsize==k]==0)
      if (left>0){
        for (s in offspring[offspringsize==k & leafindicator[offspringsize==k]==0]){
          mylist[[i]]<-list(x=1,y=card[s])
          famsize<-c(famsize,k)
          parentlist<-c(parentlist,j)
          orderlist=c(orderlist,s)
          leaves[leaves==s]<-j
          i<-i+1
        }}
          
        }
      }
  mylist[[1]]<-NULL
  return(list(mylist=mylist,nodes=cbind(orderlist[-1],parentlist[-1],famsize[-1])))
}

group_data<-function(sort.main,n){
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

bring_branch_lengths<-function(u,F){
  #u is the vector of intercoalescent times
  #given u and F, returns branch lengths ln+1,ln,ln-1,...,l3 and family sizes
  dimf<-nrow(F)
  diffM<-F[2:(dimf),2:dimf]-F[2:dimf,1:(dimf-1)]
  singletons<-seq(1,dimf)[diff(F[,1])<0]
  d<-rep(u[1],dimf) #for corresponding ln+1,ln,ln-1,...,l3
  firstzero<-rep(0,dimf-2) 
  familysize<-rep(2,dimf)#for corresponding Sn,Sn-1,..S2
  clades_offspring<-0
  clades_parents<-0
  coal_times<-cumsum(u)
  for (j in 1:(dimf-2)){
    condition<-diffM[j:(dimf-1),j]==0
    if (sum(condition)>0){ firstzero[j]<-min(seq(j,(dimf-1))[condition])}else{firstzero[j]<-dimf}
    
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
      clades_offspring<-c(clades_offspring,seq(1,dimf-2)[firstzerotrans==(dimf-i+1)])
      clades_parents<-c(clades_parents,i+1)
      
    }
    if (count==2){
      familysize[i+1]<-familysize[min(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)])]+familysize[max(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)])]
      clades_offspring<-c(clades_offspring,min(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)]),max(seq(1,(dimf-2))[firstzerotrans==(dimf-i+1)]))
      clades_parents<-c(clades_parents,i+1,i+1)
    }
    
  }
  familysize<-familysize[-dimf]
  familysize<-c(1,familysize)
  return(list(d=d,familysize=familysize,nodes=cbind(clades_offspring[-1],clades_parents[-1]),singletons=singletons))
}