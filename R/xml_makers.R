get_tips_and_times = function(tree, backwards=FALSE) {
  root_node = tree$edge[1,1]
  tip_labels = tree$tip.label
  tip_times = dist.nodes(tree)[root_node, 1:(root_node - 1)]
  
  if (backwards)
    tip_times = max(tip_times) - tip_times
  
  return(data.frame(label = tip_labels, time = tip_times))
}

#' Make taxa BEAST XML block
#'
#' @param tree a phylo tree to be processed.
#' @param units character units of time for BEAST.
#'
#' @return a character string to be pasted into the taxa block of the BEAST XML.
#' @export
#'
#' @examples
make_taxa = function(tree, units = "years") {
  tip_times = get_tips_and_times(tree)
  n = dim(tip_times)[1]
  items = list()
  for (i in 1:n) {
    items[[i]] = sprintf('\t\t<taxon id="%s">\n\t\t\t<date value="%f" direction="forwards" units="%s"/>\n\t\t</taxon>',
                         tip_times[i, 1], tip_times[i, 2], units)
  }
  return(paste('\t<taxa id="taxa">', paste0(items, collapse = '\n'), '\t</taxa>', sep = '\n'))
}

#' Make alignment BEAST XML block
#' 
#' @param nexus_data character filename of nexus file containing the sequence 
#'   data to be processed.
#'   
#' @return a character string to be pasted into the alignment block of the BEAST
#'   XML.
#' @export
#' 
#' @examples
make_alignment = function(nexus_data) {
  alig = read.nexus.data(nexus_data)
  n = length(alig)
  ids = names(alig)
  
  items = list()
  for (i in 1:n) {
    id = ids[i]
    seq = paste0(alig[[i]], collapse = "")
    seq = str_to_upper(seq)
    items[[i]] = sprintf('\t\t<sequence>\n\t\t\t<taxon idref="%s"/>\n\t\t\t%s\n\t\t</sequence>',
                         id, seq)
  }
  return(paste('\t<alignment id="alignment" dataType="nucleotide">', paste0(items, collapse = '\n'), '\t</alignment>', sep = '\n'))
}
