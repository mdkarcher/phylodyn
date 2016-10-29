#' Returns a phylo object from the argumentes generated with coalsim
#' 
#' @param args is a list containing vectors of coalescent times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}, etc. This list is the output of coalsim
#' @return A list with two elements \code{newikck} contains the tree in phylo format, \code{lables} a vector with tip labels 
#' @export
#' 
#' @examples
#' constant<-function(x){return (rep(1,length(x)))}
#' simulation1<-coalsim(0,10,constant)
#' tree1<-generate_newick(simulation1)
#' plot(tree1$newick)
generate_newick<-function(args)
{
  n<-sum(args$n_sampled)
  labels<-paste(rep("t",n),seq(1,n,1),rep("_",n),rep(args$samp_times[1],args$n_sampled[1]),sep="")
  
  #we could chose labels at random to coalesce but since the process is exchangeable, we don't care. At least not for now
  tb<-args$n_sampled[1] #Total branches (initial)
  s<-0 #time for branch lengths
  temp_labels<-labels[1:tb]
  temp_times<-rep(args$samp_times[1],args$n_sampled[1])
  initial.row<-2
  args2<-gen_INLA_args(args$samp_times,args$n_sampled,args$coal_times)
  for (j in 2:length(args2$event))
  {
    if (args2$event[j]==1)
    {
      s<-args2$s[j]; 
      ra<-sample(tb,1) #choose at random one of them, the other is the one to the right so not really random
      if (ra<tb)
      {
        new_label<-paste("(",temp_labels[ra],":",s-temp_times[ra],",",temp_labels[ra+1],":",s-temp_times[ra+1],")",sep="")
        temp_labels[ra]<-new_label
        temp_labels<-temp_labels[-(ra+1)]
        temp_times[ra]<-s
        temp_times<-temp_times[-(ra+1)]
      }
      else
      {
        new_label<-paste("(",temp_labels[ra],":",s-temp_times[ra],",",temp_labels[1],":",s-temp_times[1],")",sep="")
        temp_labels[1]<-new_label
        temp_labels<-temp_labels[-(ra)]
        temp_times[1]<-s
        temp_times<-temp_times[-(ra)]
      }
      tb<-tb-1
    }
    else
    { #I will be adding samples at 
      s<-args2$s[j]; 
      if (args$n_sample[initial.row]==1)
      {
        temp_labels<-c(temp_labels,labels[cumsum(args$n_sampled)[initial.row]])
        initial.row<-initial.row+1
        tb<-tb+1
        temp_times<-c(temp_times,s)        
      }else{
        end<-cumsum(args$n_sampled)[initial.row]
        ini<-cumsum(args$n_sampled)[initial.row-1]+1
        for (k in ini:end){
          temp_labels<-c(temp_labels,labels[k])
          tb<-tb+1
          temp_times<-c(temp_times,s)      
        }
        initial.row<-initial.row+1
      }
    }
  }  
  
  out.tree<-ape::read.tree(text=paste(temp_labels,";",sep=""))
  return(list(newick=out.tree,labels=labels))
}

generate_newick_old<-function(args,sample)
{
  n<-sum(sample[,1])
  labels<-paste(rep("t",n),seq(1,n,1),rep("_",n),rep(sample[,2],sample[,1]),sep="")
  
  #we could chose labels at random to coalesce but since the process is exchangeable, we don't care. At least not for now
  tb<-sample[1,1] #Total branches (initial)
  s<-0 #time for branch lengths
  temp_labels<-labels[1:tb]
  temp_times<-rep(sample[1,2],sample[1,1])
  initial.row<-2
  for (j in 2:length(args$event))
  {
    if (args$event[j]==1)
    {
      s<-args$s[j]; 
      ra<-sample(tb,1) #choose at random one of them, the other is the one to the right so not really random
      if (ra<tb)
      {
        new_label<-paste("(",temp_labels[ra],":",s-temp_times[ra],",",temp_labels[ra+1],":",s-temp_times[ra+1],")",sep="")
        temp_labels[ra]<-new_label
        temp_labels<-temp_labels[-(ra+1)]
        temp_times[ra]<-s
        temp_times<-temp_times[-(ra+1)]
      }
      else
      {
        new_label<-paste("(",temp_labels[ra],":",s-temp_times[ra],",",temp_labels[1],":",s-temp_times[1],")",sep="")
        temp_labels[1]<-new_label
        temp_labels<-temp_labels[-(ra)]
        temp_times[1]<-s
        temp_times<-temp_times[-(ra)]
      }
      tb<-tb-1
    }
    else
    { #I will be adding samples at 
      s<-args$s[j]; 
      if (sample[initial.row,1]==1)
      {
        temp_labels<-c(temp_labels,labels[cumsum(sample[,1])[initial.row]])
        initial.row<-initial.row+1
        tb<-tb+1
        temp_times<-c(temp_times,s)        
      }
    }
  }
  out.tree<-ape::read.tree(text=paste(temp_labels,";",sep=""))
  return(list(newick=out.tree,labels=labels))
}


##Generates xml file needed to run BEAST. Note that BEAUTI doesn't allow you to generate these xml files
#since we don't have an alignment.
#model can take the following values
#model= 1 is Bayesian Skyline plot with 10 change points and Exponential priors (like in the original paper)
#model= 2 is Bayesian Skytrack
#model =3 is Bayesian Skygrid
#model =4 is SMC

# write.xml<-function(...,file = "",model) 
# {
#   #TODO Adapt for multiple "unlinked" trees
#   obj <- list(...)
#   if (length(obj) == 1)
#   {
#     if (class(obj[[1]]) == "phylo")
#     {
#       ntree <- 1
#     }
#     else
#     {
#       obj <- obj[[1]]
#       ntree <- length(obj)
#     }
#   }
#   else
#   {
#     ntree <- length(obj)
#   }
#   # TODO: Remove stringr dependency
#   library(stringr)
#   dates<-str_split_fixed(obj[[1]]$tip.label,"_",n=2)[,2]
#   #The labels
#   cat("<?xml version='1.0' standalone='yes'?>\n", file = file)
#   cat(paste("<!--Generated by phylodyn ", date(), "-->\n\n", sep = ""), 
#       file = file, append = TRUE)
#   cat(paste("<beast>\n", sep = ""), 
#       file = file, append = TRUE)
#   N <- length(obj[[1]]$tip.label)
#   n<-N
#   cat(paste("<!--ntax=", N, "-->\n\n", sep = ""), 
#       file = file, append = TRUE)
#   cat("\t<taxa id='taxa'>\n", file = file, append = TRUE)
#   for (j in 1:N)
#   {
#     cat(paste("\t\t<taxon id='",obj[[1]]$tip.label[j],"'>\n\t\t\t<date value='",dates[j],"' direction='backwards' units='years'/>\n\t\t</taxon>\n", sep = ""), file = file, 
#         append = TRUE)
#   }
#   cat("\t</taxa>\n\n", file = file, append = TRUE)
#   ###The tree
# 
#   
#   for (i in 1:ntree)
#   {
#     ##TODO: Adjust headers when more than one tree. We will need to update it for GP-BEAST
#     cat("<newick id='startingTree' usingDates='true' units='substitutions'>\n", file = file, append = TRUE)
#     #if (class(obj[[i]]) != "phylo") 
#     #  next
#     #root.tag <- if (is.rooted(obj[[i]])) 
#     #  "= [&R] "
#     #else "= [&U] "
#     #cat("\tTREE *", title[i], root.tag, file = file, append = TRUE)
#     cat(write.tree(obj[[i]], file = ""), "\n", sep = "", 
#         file = file, append = TRUE)
#     cat("</newick>\n\n\n", file = file, append = TRUE)
#   }
#   cat("<treeModel id='treeModel'>\n\t<coalescentTree idref='startingTree'/>\n\t<rootHeight>\n\t\t<parameter id='treeeModel.rootHeight'/>\n\t</rootHeight>\n", file = file, append = TRUE)
#   cat("\t<nodeHeights internalNodes='true'>\n\t\t<parameter id='treeModel.internalNodeHeights'/>\n\t</nodeHeights>\n\t<nodeHeights internalNodes='true' rootNode='true'>\n\t\t<parameter id='treeModel.allInternalNodeHeights'/>\n\t</nodeHeights>\n</treeModel>\n\n", file = file, append = TRUE)
#   if (model==1)
#   {
#     cat("<generalizedSkyLineLikelihood id='skyline' linear='true'>\n\t<populationSizes>\n\t\t<parameter id='skyline.popSize' dimension='11' value='1.2' lower='0.0'/>\n\t</populationSizes>\n\t<groupSizes>\n\t\t<parameter id='skyline.groupSize' dimension='10'/>\n\t</groupSizes>\n\t<populationTree>\n\t\t<treeModel idref='treeModel'/>\n\t</populationTree>\n</generalizedSkyLineLikelihood>\n\n", file = file, append = TRUE)
#     cat("<exponentialMarkovLikelihood id='eml1' jeffreys='true'>\n\t<chainParameter>\n\t\t<parameter idref='skyline.popSize'/>\n\t</chainParameter>\n</exponentialMarkovLikelihood>\n\n", file = file, append = TRUE)
#     #cat("<strictClockBranchRates id='branchRates'>\n\t<rate>\n\t\t<parameter id='clock.rate' value='1.0'/>\n\t </rate>\n</strictClockBranchRates>\n\n", file = file, append = TRUE)
#     #cat("<!--not sure it is needed-->\n\n<HKYModel id='hky'>\n\t<frequencies>\n\t\t<frequencyModel dataType='nucleotide'>\n\t\t\t<frequencies>\n\t\t\t\t<parameter id='frequencies' value='0.25 0.25 0.25 0.25'/>\n\t\t\t</frequencies>\n\t\t</frequencyModel>\n\t </frequencies>\n\t<kappa>\n\t\t<parameter id='kappa' value='2.0' lower='0.0'/>\n\t</kappa>\n</HKYModel>\n\n",file=file,append=TRUE)
#     #cat("<!--not sure it is needed-->\n\n<siteModel id='siteModel'>\n\t<substitutionModel>\n\t\t<hkyModel idref='hky'/>\n\t</substitutionModel>\n</siteModel>\n\n",file=file,append=TRUE)
#     cat("<!-- Define operators. Note we only want to sample popsizes -->\n\n",file=file,append=TRUE)
#     cat("<operators id='operators'>\n\t<scaleOperator scaleFactor='0.75' weight='15'>\n\t\t<parameter idref='skyline.popSize'/>\n\t</scaleOperator>\n\t<deltaExchange delta='1' integer='true' weight='6' autoOptimize='false'>\n\t\t<parameter idref='skyline.groupSize'/>\n\t</deltaExchange>\n</operators>\n\n",file=file,append=TRUE)
#     cat("<mcmc id='mcmc' chainLength='1500000' autoOptimize='true'>\n",file=file,append=TRUE)
#     cat("\t<posterior id='posterior'>\n\t\t<prior id='prior'>\n\t\t\t<uniformPrior lower='0.0' upper='1.0E100'>\n\t\t\t\t<parameter idref='skyline.popSize'/>\n\t\t\t</uniformPrior>\n\t\t\t<exponentialMarkovLikelihood idref='eml1'/>\n\t\t</prior>\n\t\t<likelihood id='likelihood'>\n\t\t\t<generalizedSkyLineLikelihood idref='skyline'/>\n\t\t</likelihood>\n\t</posterior>\n\n",file=file,append=TRUE)
#     cat("\t<operators idref='operators'/>\n\t<log id='screenLog' logEvery='1000'>\n\t\t<column label='Posterior' dp='4' width='12'>\n\t\t\t<posterior idref='posterior'/>\n\t\t</column>\n\t\t <column label='Prior' dp='4' width='12'>\n\t\t\t<prior idref='prior'/>\n\t\t</column>\n\t\t<column label='Likelihood' dp='4' width='12'>\n\t\t\t<likelihood idref='likelihood'/>\n\t\t</column>\n\t</log>\n\n",file=file,append=TRUE)
#     cat(paste("\t<log id='fileLog' logEvery='100' fileName='",file,".log'>\n",sep=""),file=file,append=TRUE)
#     cat("\t\t<posterior idref='posterior'/>\n\t\t<prior idref='prior'/>\n\t\t<likelihood idref='likelihood'/>\n\t\t<parameter idref='skyline.popSize'/>\n\t\t<parameter idref='skyline.groupSize'/>\n\t\t<generalizedSkyLineLikelihood idref='skyline'/>\n\t</log>\n</mcmc>\n\t<report>\n\t\t<property name='timer'>\n\t\t\t<mcmc idref='mcmc'/>\n\t\t</property>\n\t</report>\n</beast>\n",file=file,append=TRUE)
#   }
#   if (model==2)
#   {
#     cat(paste("<gmrfSkyrideLikelihood id='skyride' timeAwareSmoothing='true' randomizeTree='false'>\n\t<populationSizes>\n\t\t<parameter id='skyride.logPopSize' dimension='",n-1,"' value='0.2'/>\n",sep=""),file=file,append=TRUE)
#     cat("\t</populationSizes>\n\t<precisionParameter>\n\t\t<parameter id='skyride.precision' value='1.0' lower='0.0'/>\n\t</precisionParameter>\n\t<populationTree>\n\t\t<treeModel idref='treeModel'/>\n\t</populationTree>\n</gmrfSkyrideLikelihood>",file=file,append=TRUE)
#     cat("\n<!-- Define operators. Note we only want to sample popsizes -->\n\n",file=file,append=TRUE)
#     cat("<operators id='operators'>\n\t<gmrfBlockUpdateOperator scaleFactor='2.0' weight='2'>\n\t\t<gmrfSkyrideLikelihood idref='skyride'/>\n\t</gmrfBlockUpdateOperator>\n</operators>\n\n",file=file,append=TRUE)
#     cat("<mcmc id='mcmc' chainLength='1500000' autoOptimize='true'>\n",file=file,append=TRUE)
#     cat("\t<posterior id='posterior'>\n\t\t<prior id='prior'>\n\t\t\t<gammaPrior shape='0.0010' scale='1000.0' offset='0.0'>\n\t\t\t\t<parameter idref='skyride.precision'/>\n\t\t\t</gammaPrior>\n\t\t</prior>\n\t\t<likelihood id='likelihood'>\n\t\t\t<gmrfSkyrideLikelihood idref='skyride'/>\n\t\t</likelihood>\n\t</posterior>\n\n",file=file,append=TRUE)
#     cat("\t<operators idref='operators'/>\n\t<log id='screenLog' logEvery='1000'>\n\t\t<column label='Posterior' dp='4' width='12'>\n\t\t\t<posterior idref='posterior'/>\n\t\t</column>\n\t\t <column label='Prior' dp='4' width='12'>\n\t\t\t<prior idref='prior'/>\n\t\t</column>\n\t\t<column label='Likelihood' dp='4' width='12'>\n\t\t\t<likelihood idref='likelihood'/>\n\t\t</column>\n\t</log>\n\n",file=file,append=TRUE)
#     cat(paste("\t<log id='fileLog' logEvery='100' fileName='",file,".log'>\n",sep=""),file=file, append=TRUE)
#     cat("\t\t<posterior idref='posterior'/>\n\t\t<prior idref='prior'/>\n\t\t<likelihood idref='likelihood'/>\n\t\t<parameter idref='skyride.precision'/>\n\t\t<parameter idref='skyride.logPopSize'/>\n\t\t<gmrfSkyrideLikelihood idref='skyride'/>\n\t</log>\n</mcmc>\n\t<report>\n\t\t<property name='timer'>\n\t\t\t<mcmc idref='mcmc'/>\n\t\t</property>\n\t</report>\n</beast>\n",file=file,append=TRUE)
#   }
#   if (model==3)
#   {
#     cat(paste("<gmrfSkyGridLikelihood id='skyride'>\n\t<populationSizes>\n\t\t<parameter id='skyride.logPopSize' dimension='",n-1,"' value='0.2'/>\n",sep=""),file=file,append=TRUE)
#     cat(paste("\t</populationSizes>\n\t<numGridPoints>\n\t\t<parameter id='skyride.numGridPoints' value='",n-2,"'/>\n",sep=""),file=file,append=TRUE)
#     cat("\t</numGridPoints>\n\t<cutOff>\n\t\t<parameter id='skyride.cutOff' value='2.9'/>\n\t</cutOff>\n\t<precisionParameter>\n\t\t<parameter id='skyride.precision' value='0.1' lower='0.0'/>\n\t</precisionParameter>\n\t<populationTree>\n\t\t<treeModel idref='treeModel'/>\n\t</populationTree>\n</gmrfSkyGridLikelihood>",file=file,append=TRUE)
#     cat("\n\n<!-- Define operators. Note we only want to sample popsizes -->\n\n",file=file,append=TRUE)
#     cat("<operators id='operators'>\n\t<gmrfGridBlockUpdateOperator scaleFactor='2.0' weight='2'>\n\t\t<gmrfSkyGridLikelihood idref='skyride'/>\n\t</gmrfGridBlockUpdateOperator>\n</operators>\n\n",file=file,append=TRUE)
#     cat("<mcmc id='mcmc' chainLength='1500000' autoOptimize='true'>\n",file=file,append=TRUE)
#     cat("\t<posterior id='posterior'>\n\t\t<prior id='prior'>\n\t\t\t<gammaPrior shape='0.0010' scale='1000.0' offset='0.0'>\n\t\t\t\t<parameter idref='skyride.precision'/>\n\t\t\t</gammaPrior>\n\t\t</prior>\n\t\t<likelihood id='likelihood'>\n\t\t\t<gmrfSkyGridLikelihood idref='skyride'/>\n\t\t</likelihood>\n\t</posterior>\n\n",file=file,append=TRUE)
#     cat("\t<operators idref='operators'/>\n\t<log id='screenLog' logEvery='1000'>\n\t\t<column label='Posterior' dp='4' width='12'>\n\t\t\t<posterior idref='posterior'/>\n\t\t</column>\n\t\t <column label='Prior' dp='4' width='12'>\n\t\t\t<prior idref='prior'/>\n\t\t</column>\n\t\t<column label='Likelihood' dp='4' width='12'>\n\t\t\t<likelihood idref='likelihood'/>\n\t\t</column>\n\t</log>\n\n",file=file,append=TRUE)
#     cat(paste("\t<log id='fileLog' logEvery='100' fileName='",file,".log'>\n",sep=""),file=file, append=TRUE)
#     cat("\t\t<posterior idref='posterior'/>\n\t\t<prior idref='prior'/>\n\t\t<likelihood idref='likelihood'/>\n\t\t<parameter idref='skyride.precision'/>\n\t\t<parameter idref='skyride.logPopSize'/>\n\t\t<gmrfSkyGridLikelihood idref='skyride'/>\n\t</log>\n</mcmc>\n\t<report>\n\t\t<property name='timer'>\n\t\t\t<mcmc idref='mcmc'/>\n\t\t</property>\n\t</report>\n</beast>\n",file=file,append=TRUE)
#   }
# }
