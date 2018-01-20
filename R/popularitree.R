#' Get a smaller tree with popular taxa spanning the diversity
#'
#' This assumes you will have a data.frame of popularity measures
#' and then sum to combine the columns into a single measure.
#' @param phy A phylo object
#' @param ntax The final number of taxa you want
#' @param measures Data.frame of the measures of popularity
#' @param pd_weighting How much to weight the pd vs the popularity
#' @param nrep How many replicates to use
#' @return A phylo object pruned to the number of taxa
#' @examples
#' big.phy <- rotl::get_study("ot_308")[[1]]
#' big.phy <- geiger::drop.random(big.phy, ape::Ntip(phy) - 20)
#' measures <- tree_compute_measures(big.phy)
#' small.phy <- popularitree(phy=big.phy, ntax=5, measures, nrep=100)
#' par(mfcol=c(1,3))
#' plot(big.phy)
#' plot(small.phy)
#' plot(tree_convert_to_common(small.phy))
#' @export
popularitree <- function(phy, ntax, measures, pd_weighting=.5, nrep=10000) {
  tree.scaling <- (ntax/ape::Ntip(phy))/sum(phy$edge.length) # want the measures and the tree brlen to be on the same scale after sampling
  phy$edge.length <- phy$edge.length*tree.scaling
  taxon_scores <- apply(measures, 1, sum)
  taxon_scores <- taxon_scores*(ntax/ape::Ntip(phy))/sum(taxon_scores)
  names(taxon_scores) <- rownames(measures)
  population <- unique(replicate(nrep, generate_sample(taxon_scores, ntax, phy, pd_weighting)))
  best_index <- which.max(unlist(population['final_score',]))
  best_phy <- ape::drop.tip(phy, phy$tip.label[!phy$tip.label %in% population['selected',best_index]$selected])
  best_phy$edge.length <- best_phy$edge.length/tree.scaling #restore brlen
  return(best_phy)
}

#' Generate a sample of species, weighted by their score
#' @param taxon_scores The vector of scores for taxa
#' @param ntax The number of taxa you want returned
#' @param phy Phylogeny to prune
#' @param pd_weighting How much to weight the pd vs the popularity
#' @param min.prob The minimum probability of sampling for a taxon
#' @param max.prob The maximum probability of sampling for a taxon
#' @return A list: selected=vector of names of sampled taxa to include, and score=sum of the taxon metrics for these taxa
#' @export
generate_sample <- function(taxon_scores, ntax, phy, pd_weighting, min.prob=0.01, max.prob=0.99) {
  taxon_scores[is.na(taxon_scores)] <- 0
  prob <- taxon_scores/sum(taxon_scores)
  prob[prob<min.prob] <- min.prob
  prob[prob>max.prob] <- max.prob
  names.to.choose <- sample(names(taxon_scores), size=ntax, replace=FALSE, prob=prob)
  pd_sample <- data.frame(matrix(0, nrow=1, ncol=length(taxon_scores)))
  colnames(pd_sample) <- names(taxon_scores)
  pd_sample[1,names.to.choose] <- 1
  pd <- picante::pd(pd_sample, phy, include.root=FALSE)[1,1]
  score <- sum(taxon_scores[names.to.choose])
  final_score <- pd_weighting*pd + (1-pd_weighting)*score
  return(list(selected=names.to.choose, score=score, pd=pd, final_score=final_score))
}

#' Get a sum of hits of a species relative to a standard
#' @param species The focal species to look at search interest over
#' @param standard Since Google Trends always compares popularity, give a standard species to use
#' @param verbose If TRUE, chat about what species is being examined
#' @return The sum of relative hits over the time period
#' @export
species_gtrends_popularity <- function(species, standard="Myrmeocystus mexicanus", verbose=TRUE) {
  if(verbose) {
    print(paste0("Retrieving trends data for taxon '", species, "'"))
  }
  result <- gtrendsR::gtrends(gsub("_", " ",c(species, standard)), time="2017-01-01 2017-12-31")$interest_over_time
  result <- subset(result, result$keyword==species)
  return(sum(result$hits))
}

#' Return 1 if a species has a common name in NCBI, 0 otherwise
#' @param species The focal species to examine
#' @return 1 or 0
#' @examples
#' species_ncbi_common_exists('Puma concolor')
#' species_ncbi_common_exists('Pomatomus saltatrix')
#'
#' @export
species_ncbi_common_exists <- function(species) {
  returned.name <- NULL
  has_common <- 0
  try(returned.name <- taxize::sci2comm(scinames=species, db='ncbi'))
  if(!is.null(returned.name)) {
    has_common <- length(nchar(returned.name[[1]]))
  }
  return(has_common)
}

#' Return common name if known, otherwise return scientific name
#' @param species Species name
#' @return Common name (ideally) or species name
#' @export
species_get_common <- function(species) {
  returned.name <- NULL
  has_common <- 0
  try(returned.name <- taxize::sci2comm(scinames=species, db='ncbi'))
  if(!is.null(returned.name)) {
    has_common <- length(nchar(returned.name[[1]]))
  }
  final_name <- species
  if(has_common==1) {
    final_name <- unname(returned.name[[1]])
  }
  return(final_name)
}

#' Compute popularity measures
#' @param phy An input phylogeny
#' @return data.frame of measures
#' @export
tree_compute_measures <- function(phy) {
  taxa <- phy$tip.label
  common <- sapply(taxa, species_ncbi_common_exists)
  common <- common/max(1,sum(common)) #to normalize the different measures
  trends <- sapply(taxa, species_gtrends_popularity)
  trends <- trends/max(1,sum(trends))  #to normalize the different measures
  result <- data.frame(common=common, trends=trends)
  rownames(result) <- taxa
  return(result)
}

#' Convert tips to common names if possible
#' @param phy A phylogeny
#' @return Phylogeny with names changed
#' @export
tree_convert_to_common <- function(phy) {
  phy$tip.label <- gsub(" ", "_",sapply(phy$tip.label,species_get_common))
  return(phy)

}
