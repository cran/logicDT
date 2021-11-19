#' @importFrom stats runif rnorm
simulateGRS <- function(n.obs, n.snps, pathway.interactions, pathway.genotypes,
                        maf, outer.beta0, pathway.beta0, pathway.beta,
                        exposure.dist, exposure.dist.params, exposure.dist.probs,
                        outer.interactions, outer.int.directions, outer.beta,
                        outcome = "binary", variance = 1) {
  snp.data <- data.frame(matrix(ncol = 0, nrow = n.obs))
  grs.data <- data.frame(matrix(ncol = 0, nrow = n.obs))
  exp.data <- data.frame(matrix(ncol = 0, nrow = n.obs))
  y <- rep(0, n.obs)
  if(length(maf) == 1)
    maf <- rep(maf, n.snps)

  for(j in 1:n.snps) {
    prob.table <- c((1-maf[j])^2, 2*maf[j]*(1-maf[j]), maf[j]^2)
    snp.data[, paste("SNP", j, sep="")] <- sample(1:3, n.obs, replace = TRUE, prob = prob.table)
  }

  n.pathways <- length(pathway.interactions)
  n.outer.terms <- length(outer.interactions)

  for(i in 1:n.pathways) {
    effect.frame <- matrix(1, nrow = n.obs, ncol = length(pathway.interactions[[i]]) + 1)
    for(j in 1:length(pathway.interactions[[i]])) {
      sub_data <- snp.data[, pathway.interactions[[i]][[j]], drop=FALSE]
      for(k in 1:length(pathway.interactions[[i]][[j]])) {
        if(pathway.genotypes[[i]][[j]][k] < 0) {
          sub_data[, k] <- sub_data[, k] != -pathway.genotypes[[i]][[j]][k]
        } else {
          sub_data[, k] <- sub_data[, k] == pathway.genotypes[[i]][[j]][k]
        }
      }
      effect.frame[, j+1] <- as.numeric(Reduce("&", sub_data))
    }
    grs.data[, paste("GRS", i, sep="")] <- as.vector(effect.frame %*% c(pathway.beta0[i], pathway.beta[[i]]))
  }

  n.exposure <- length(exposure.dist)
  if(n.exposure > 0) {
    for(j in 1:n.exposure) {
      n.inner.dists <- length(exposure.dist[[j]])
      exp.chosen.dists <- sample(1:n.inner.dists, n.obs, replace = TRUE, prob = exposure.dist.probs[[j]])
      for(k in 1:n.inner.dists) {
        if(exposure.dist[[j]][k] == "uniform") {
          exp.data[exp.chosen.dists == k, paste("E", j, sep="")] <- runif(sum(exp.chosen.dists == k), exposure.dist.params[[j]][[k]][1], exposure.dist.params[[j]][[k]][2])
        } else if(exposure.dist[[j]][k] == "normal") {
          exp.data[exp.chosen.dists == k, paste("E", j, sep="")] <- rnorm(sum(exp.chosen.dists == k), exposure.dist.params[[j]][[k]][1], exposure.dist.params[[j]][[k]][2])
        }
      }
    }
  }

  lin.pred <- outer.beta0

  for(i in 1:n.outer.terms) {
    current.term <- outer.interactions[[i]]
    current.grs <- current.term[current.term > 0]
    current.exp <- -current.term[current.term < 0]
    current.int.direction <- outer.int.directions[[i]]
    if(is.null(current.int.direction))
      current.int.frame <- cbind(grs.data[, current.grs, drop=FALSE], exp.data[, current.exp, drop=FALSE])
    else {
      grs.directions <- current.int.direction[match(current.grs, current.term)]
      exp.directions <- current.int.direction[match(-current.exp, current.term)]
      current.grs.data <- grs.data[, current.grs, drop=FALSE]
      current.exp.data <- exp.data[, current.exp, drop=FALSE]

      current.grs.data[current.grs.data < 0, match("+", grs.directions)[!is.na(match("+", grs.directions))]] <- 0
      current.grs.data[current.grs.data > 0, match("-", grs.directions)[!is.na(match("-", grs.directions))]] <- 0
      current.exp.data[current.exp.data < 0, match("+", exp.directions)[!is.na(match("+", exp.directions))]] <- 0
      current.exp.data[current.exp.data > 0, match("-", exp.directions)[!is.na(match("-", exp.directions))]] <- 0

      current.int.frame <- cbind(current.grs.data, current.exp.data)
    }
    lin.pred <- lin.pred + outer.beta[i] * Reduce("*", current.int.frame)
  }

  if(outcome == "binary") {
    probs <- 1/(1+exp(-lin.pred))
    for(i in 1:n.obs) {
      y[i] <- sample(0:1, 1, prob = c(1-probs[i], probs[i]))
    }
  } else if (outcome == "continuous") {
    for(i in 1:n.obs) {
      y[i] <- rnorm(1, lin.pred[i], sd = sqrt(variance))
    }
  }

  ret_obj <- list(y = y, snp.data = snp.data, grs.data = grs.data, exp.data = exp.data, lin.pred = lin.pred)
  class(ret_obj) <- "sim.grs.data"
  return(ret_obj)
}


# data <- simulateGRS(2000, 50,
#                     list(list(1, 2, 3, c(4, 5)), list(26, 27, 28, c(26, 29))),
#                     list(list(-1, -1, -1, c(-1, -1)), list(-1, -1, -1, c(-1, -1))),
#                     0.35, log(0.3/(1-0.3)), log(c(0.2/(1-0.2), 0.2/(1-0.2))),
#                     list(log(c(1.2, 1.2, 1.2, 1.8)), log(c(1.2, 1.2, 1.2, 1.8))),
#                     list(c("normal", "normal")), list(list(c(30, 5), c(20, 5))), list(c(0.5, 0.5)),
#                     list(1, c(1,2), c(1,-1)), log(c(1.2, 1.5, 1.8)),
#                     "binary")


