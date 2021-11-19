#' @importFrom stats runif
simulateSNPtree <- function(n_obs, n_snps, maf, prevalence, effect_terms, genotypes, odds_ratios) {
  data <- data.frame(y = rep(0, n_obs))
  if(length(maf) == 1)
    maf <- rep(maf, n_snps)
  # if(length(genotype) == 1)
  #   genotypes <- rep(genotype, length(effect_snps))

  for(j in 1:n_snps) {
    prob_table <- c((1-maf[j])^2, 2*maf[j]*(1-maf[j]), maf[j]^2)
    data[, paste("SNP", j, sep="")] <- sample(1:3, n_obs, replace = TRUE, prob = prob_table)
  }

  x <- odds_ratios * prevalence/(1-prevalence)
  probs <- x/(1+x)

  effect_frame <- data.frame(matrix(nrow = n_obs, ncol = 0))

  for(j in 1:length(effect_terms)) {
    sub_data <- data[, effect_terms[[j]] + 1, drop=FALSE]
    for(k in 1:length(effect_terms[[j]])) {
      if(genotypes[[j]][k] < 0) {
        sub_data[, k] <- sub_data[, k] != -genotypes[[j]][k]
      } else {
        sub_data[, k] <- sub_data[, k] == genotypes[[j]][k]
      }
    }
    effect_frame[, j] <- as.numeric(Reduce("&", sub_data))
  }

  for(i in 1:n_obs) {
    prob <- probs[t(as.numeric(effect_frame[i,]) + 1)]
    if(prob > runif(1))
      data$y[i] <- 1
  }

  return(data)
}

# data <- simulateSNPtree(2000, 20, 0.35, 0.4, list(c(4,5), c(1,6,7)), list(c(-1,-1), c(-1,-1,-1)),
#                         aperm(array(c(1.0, 2.0, 2.0, 3.0), rep(2, 2))))
