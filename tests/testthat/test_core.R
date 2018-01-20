  test_that("doing a single analysis", {
    phy <- ape::rcoal(12)
    measures <- data.frame(a=runif(12), b=runif(12))
    rownames(measures) <- phy$tip.label
    taxon_scores <- apply(measures, 1, sum)
    result <- generate_sample(taxon_scores, 7, phy, pd_weighting=0.5)
    expect_equal("list", class(result))
    expect_gte(result$score, 1)
    expect_gte(result$pd, 0.1)
    expect_gte(result$final_score, 0.1)
  })

test_that("can optimize with popularity", {
  phy <- ape::rcoal(50)
  ntax_desired <- 20
  measures <- data.frame(a=runif(50), b=runif(50))
  rownames(measures) <- phy$tip.label
  result <- popularitree(phy, ntax_desired, measures, pd_weighting=.5, nrep=1000)
  expect_equal("phylo", class(result))
  expect_equal(ape::Ntip(result), ntax_desired)
})
