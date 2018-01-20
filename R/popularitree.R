
popularitree <- function(phy, ntax, measures, formula, weighting=1) {
  taxon.metrics <- with(measures, formula)
  names(taxon.metrics) <- rownames(measures)
  #an optimization occurs: generate a set of taxa and get the score for the included set, sample the tree to get the phylogenetic diversity, combine the scores using the weighting, repeat, repeat, repeat.
  return(phy.best)
}

#' Generate a sample of species, weighted by their score
#' @param taxon.metrics The vector of scores for taxa
#' @param ntax The number of taxa you want returend
#' @param min.prob The minimum probability of sampling for a taxon
#' @param max.prob The maximum probability of sampling for a taxon
#' @return A list: selected=vector of names of sampled taxa to include, and score=sum of the taxon metrics for these taxa
#' @export
generate.sample <- function(taxon.metrics, ntax, min.prob=0.01, max.prob=0.99) {
  taxon.metrics[is.na(taxon.metrics)] <- 0
  prob <- taxon.metrics/sum(taxon.metrics)
  prob[prob<min.prob] <- min.prob
  prob[prob>max.prob] <- max.prob
  names.to.choose <- sample(names(taxon.metrics, size=ntax, replace=FALSE, prob=prob))
  return(list(selected=names.to.choose, score=sum(taxon.metrics[names.to.choose])))
}
