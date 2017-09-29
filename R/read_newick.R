#' Summarize a phylogeny.
#' 
#' @param phy a \code{phylo} object containing a phylogeny.
#'   
#' @return A list containing vectors of sampling times \code{samp_times}, number 
#'   sampled per sampling time \code{n_sampled}, and coalescent times
#'   \code{coal_times}.
#' @export
#' 
#' @examples
#' data("NY_flu")
#' summarize_phylo(NY_flu)
summarize_phylo <- function(phy)
{
  hgpstat <- heterochronous_gp_stat(phy)
  return(list(samp_times = hgpstat$samp_times,
              n_sampled  = hgpstat$n_sampled,
              coal_times = hgpstat$coal_times))
}

summarize_phylo2 <- function(phy, backwards = TRUE)
{
  if (class(phy) != "phylo")
    stop("object \"phy\" is not of class \"phylo\"")
  
  n_nodes = phy$Nnode
  n_tips = length(phy$tip.label)
  
  root_node <- phy$edge[1,1]
  raw_times <- dist.nodes(phy)[root_node, ]
  raw_samp_times <- head(raw_times, n_tips)
  raw_coal_times <- tail(raw_times, n_nodes)
  
  if (backwards)
  {
    raw_coal_times <- max(raw_samp_times) - raw_coal_times
    raw_samp_times <- max(raw_samp_times) - raw_samp_times
  }
  
  samp_tab <- table(raw_samp_times)
  samp_times <- as.numeric(names(samp_tab))
  n_sampled <- as.numeric(samp_tab)
  
  coal_times <- sort(as.numeric(raw_coal_times))
  
  return(list(samp_times = samp_times,
              n_sampled  = n_sampled,
              coal_times = coal_times))
}

branching_sampling_times <- function(phy)
{
  phy = ape::new2old.phylo(phy)

  if (class(phy) != "phylo")
    stop("object \"phy\" is not of class \"phylo\"")

  tmp <- as.numeric(phy$edge)
  nb.tip <- max(tmp)
  nb.node <- -min(tmp)
  xx <- as.numeric(rep(NA, nb.tip + nb.node))
  names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
  xx["-1"] <- 0

  for (i in 2:length(xx))
  {
    nod <- names(xx[i])
    ind <- which(phy$edge[, 2] == nod)
    base <- phy$edge[ind, 1]
    xx[i] <- xx[base] + phy$edge.length[ind]
  }

  depth <- max(xx)
  branching_sampling_times <- depth - xx
  
  return(branching_sampling_times)
}

heterochronous_gp_stat <- function(phy, tol=0.0)
{
  #Update Aug 2015 by Julia. Adhoc for simulation with a tolerance parameters
  b.s.times = branching_sampling_times(phy)
  int.ind = which(as.numeric(names(b.s.times)) < 0)
  tip.ind = which(as.numeric(names(b.s.times)) > 0)
  num.tips = length(tip.ind)
  num.coal.events = length(int.ind)
  sampl.suf.stat = rep(NA, num.coal.events)
  coal.interval = rep(NA, num.coal.events)
  coal.lineages = rep(NA, num.coal.events)
  sorted.coal.times = sort(b.s.times[int.ind])
  names(sorted.coal.times) = NULL
  sampling.times = sort((b.s.times[tip.ind]))
  
  for (i in 2:length(sampling.times))
  {
    if ((sampling.times[i] - sampling.times[i - 1]) < tol)
    {
      sampling.times[i] <- sampling.times[i - 1]
    }
  }
  
  unique.sampling.times <- unique(sampling.times)
  sampled.lineages = NULL
  
  for (sample.time in unique.sampling.times)
  {
    sampled.lineages = c(sampled.lineages, sum(sampling.times == sample.time))
  }
  
  return(list(coal_times = sorted.coal.times,
              samp_times = unique.sampling.times,
              n_sampled = sampled.lineages))
}

heterochronous_gp_stat_old <- function(phy)
{
  b.s.times = branching_sampling_times(phy)
  int.ind = which(as.numeric(names(b.s.times)) < 0)
  tip.ind = which(as.numeric(names(b.s.times)) > 0)
  num.tips = length(tip.ind)
  num.coal.events = length(int.ind)
  sampl.suf.stat = rep(NA, num.coal.events)
  coal.interval = rep(NA, num.coal.events)
  coal.lineages = rep(NA, num.coal.events)
  sorted.coal.times = sort(b.s.times[int.ind])
  names(sorted.coal.times) = NULL
  #unique.sampling.times = sort(unique(b.s.times[tip.ind]))
  sampling.times = sort((b.s.times[tip.ind]))

  for (i in 2:length(sampling.times))
  {
    if ((sampling.times[i]-sampling.times[i-1])<0.1)
    {
      sampling.times[i]<-sampling.times[i-1]
    }
  }
  unique.sampling.times<-unique(sampling.times)
  sampled.lineages = NULL
  for (sample.time in unique.sampling.times)
  {
    sampled.lineages = c(sampled.lineages, sum(sampling.times == sample.time))  
  }
  
  return(list(coal_times=sorted.coal.times, samp_times = unique.sampling.times, n_sampled=sampled.lineages))
}
