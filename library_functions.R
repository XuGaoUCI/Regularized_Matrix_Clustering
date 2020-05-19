###This script is to implement the mixture matrix normal with regularization
require(expm)
require(jpeg)
require(matlib)
require(MASS)
# library functions
generateCrossing <- function(sampledim, scale) {
  dim.row <- sampledim[2]
  dim.col <- sampledim[3]
  mu.crossing <- matrix(0, sampledim[2], sampledim[3])
  mu.crossing[floor(0.45 * sampledim[2]) : floor(0.55 * sampledim[2]), 
              floor(0.25 * sampledim[3]) : floor(0.75 * sampledim[3])] <- scale
  mu.crossing[floor(0.25 * sampledim[2]) : floor(0.75 * sampledim[2]), 
              floor(0.45 * sampledim[3]) : floor(0.55 * sampledim[3])] <- scale
  return(mu.crossing)  
}
generateSquare <- function(sampledim, scale, constant = FALSE) {
  dim.row <- sampledim[2]
  dim.col <- sampledim[3]
  mu.square <- matrix(0, sampledim[2], sampledim[3])
  row.chosen <- floor(0.15 * sampledim[2]) : floor(0.85 * sampledim[2])
  if (constant == TRUE) {
    for (row.num in row.chosen) {
      mu.square[row.num, 
                floor(0.15 * sampledim[3]) : floor(0.85 * sampledim[3])] <- scale 
    }
  } else  {
    for (row.num in row.chosen) {
      mu.square[row.num, 
                floor(0.15 * sampledim[3]) : floor(0.85 * sampledim[3])] <- scale + row.num
    }
  }
  return(mu.square)  
}
generateSeparatedSquare <- function(sampledim, scale) {
  dim.row <- sampledim[2]
  dim.col <- sampledim[3]
  mu.square <- matrix(0, sampledim[2], sampledim[3])
  mu.square[floor(0.15 * sampledim[2]) : floor(0.85 * sampledim[2]), 
            floor(0.15 * sampledim[3]) : floor(0.45 * sampledim[3])] <- scale
  mu.square[floor(0.15 * sampledim[2]) : floor(0.85 * sampledim[2]), 
            floor(0.55 * sampledim[3]) : floor(0.85 * sampledim[3])] <- scale
  return(mu.square)  
}
generateTriangle <- function(sampledim, scale) {
  dim.row <- sampledim[2]
  dim.col <- sampledim[3]
  mu.triangle <- matrix(0, sampledim[2], sampledim[3])
  j <- floor(sampledim[3] * 0.35)
  for (i in 0:floor(sampledim[2]/2)) {
    mu.triangle[j, (floor(sampledim[2]/2) - i):(floor(sampledim[2]/2)+i)] <- scale + i
    j <- j + 1
  }
  return(mu.triangle)  
}
generateNormalparam <- function(cluster.nums, sampledim) {
  mus <- array(0, dim = c(cluster.nums, sampledim[2:3]))
  rowcovs <- array(0, dim = c(cluster.nums, sampledim[2], sampledim[2]))
  colcovs <- array(0, dim = c(cluster.nums, sampledim[3], sampledim[3]))
  for (i in 1:cluster.nums) {
    mus[i, , ] <- matrix(runif(sampledim[2] * sampledim[3]) + 10 * i, ncol = sampledim[3]) 
    rowcovs[i, , ] <- matrix(runif(sampledim[2]^2), ncol=sampledim[2]) 
    rowcovs[i, , ] <- t(rowcovs[i, , ]) %*% rowcovs[i, , ]
    colcovs[i, , ] <- matrix(runif(sampledim[3]^2), ncol=sampledim[3]) 
    colcovs[i, , ] <- t(colcovs[i, , ]) %*% colcovs[i, , ]
  }
  return(list(mus = mus, rowcovs = rowcovs, colcovs = colcovs))
}
generateMixMatrixnormal <- function(dim) {
  signals <- array(0, dim = dim$sampledim)
  result.indicate <- c()
  for (i in 1:dim$sampledim[1]) {
    indicate <- which(rmultinom(1, 1, dim$pi) == 1)
    sqrt.rowcovs <- sqrtm(normal.params$rowcovs[indicate, , ])
    sqrt.colcovs <- sqrtm(normal.params$colcovs[indicate, , ])
    signals[i, , ] <- dim$normal.params$mus[indicate, , ] + sqrt.rowcovs %*% matrix(rnorm(dim$sampledim[2] * dim$sampledim[3]), ncol = dim$sampledim[3]) %*%
      sqrt.colcovs
    result.indicate <- c(result.indicate, indicate)
  }
  return(list(signals = signals, indicate = result.indicate))
}
generateParamsSim <- function(cluster.nums, sampledim, scale, diag = FALSE) {
  # generate mus
  mus <- array(0, c(cluster.nums, sampledim[2], sampledim[3]))
  mus[1, , ] <- generateSquare(sampledim, scale, constant = TRUE)
  mus[2, , ] <- generateCrossing(sampledim, scale)
  rowcovs <- array(0, c(cluster.nums, sampledim[2], sampledim[2]))
  colcovs <- array(0, c(cluster.nums, sampledim[3], sampledim[2]))
  if (diag == TRUE) {
    cov.matrix <- diag(1, sampledim[2])
  } else {
    cov.matrix <- outer(1:sampledim[2], 1:sampledim[2], function(x, y) return(0.99^abs(x-y)))
  }
  for (i in 1:cluster.nums) {
    rowcovs[i, , ] <- colcovs[i, , ] <- cov.matrix
  }
  normal.params <- list(mus = mus, rowcovs = rowcovs, colcovs = colcovs)
  return(normal.params)
}
likelihood <- function(obs, mu, rowcov, colcov, log = TRUE, tol = 1e-4) {
  while (kappa(rowcov) > 1e14) {
    rowcov <- rowcov + diag(tol, nrow(rowcov))
    tol <- tol * 2
  }
  while (kappa(colcov) > 1e14) {
    colcov <- colcov + diag(tol, nrow(colcov))
    tol <- tol * 2
  }
  log.like <- - 0.5 * nrow(colcov) * determinant(rowcov, logarithm = TRUE)$modulus - 0.5 * nrow(rowcov) * determinant(colcov, logarithm = TRUE)$modulus -
    0.5 * sum(diag(solve(rowcov) %*% (obs - mu) %*% solve(colcov) %*% t(obs - mu)))
  if (log == TRUE) {
    return (log.like)
  } else {
    return(exp(log.like))
  }
}
expectationStep <- function(obses, mus, rowcovs, colcovs, pis, scale = TRUE) {
  dim.row <- length(pis)
  dim.col <- dim(obses)[1]
  pitimesfij <- matrix(0, dim.row, dim.col)
  for (i in 1:dim.row) {
    pitimesfij[i, ] <- log(pis[i]) + apply(obses, 1, likelihood, mu = mus[i,,], rowcov = rowcovs[i,,], colcov = colcovs[i,,])
  }
  if (scale == FALSE) {
    return(pitimesfij)
  }
  # rescale for each column
  max.vec <- apply(pitimesfij, 2, max)
  pitimesfij <- pitimesfij - 
    matrix(max.vec, nrow = dim.row, ncol = dim.col, byrow = TRUE)
  pitimesfij <- exp(pitimesfij)
  colsum.vec <- colSums(pitimesfij)
  gamma.matrix <- pitimesfij / matrix(colsum.vec, nrow = dim.row, ncol = dim.col, byrow = TRUE)
  return(gamma.matrix)
}
maximizationStep <- function(obses, gamma.matrix, mus, rowcovs, colcovs, pis, lambda, penalty) {
  # update mus, rowcovs and colcovs
  if (penalty == 'lasso') {
    new.mus <- updateMusLassoexact(obses, gamma.matrix, mus, rowcovs, colcovs, pis, lambda)
    #    new.mus <- updateMusexact(obses, gamma.matrix, mus, rowcovs, colcovs, pis, lambda)
    #    new.mus <- updateMus(obses, gamma.matrix, mus, rowcovs, colcovs, pis, step.size = 0.4, lambda)
  } else if (penalty == 'l2') {
    new.mus <- updateMusl2exact(obses, gamma.matrix, mus, rowcovs, colcovs, pis, lambda)
    #    new.mus <- updateMusl2(obses, gamma.matrix, mus, rowcovs, colcovs, pis, step.size = 0.4, lambda)
  } else {
    new.mus <- updateMusnuclearexact(obses, gamma.matrix, mus, rowcovs, colcovs, pis, lambda)
    # new.mus <- updateMusnuclear(obses, gamma.matrix, mus, rowcovs, colcovs, pis, step.size = 0.4, lambda)
  }
  #    updated.rowandcol <- updaterowAndcolcovs(obses, gamma.matrix, new.mus)
  updated.rowandcol <- updaterowAndcolcovsTraceConstrain(obses, gamma.matrix, new.mus)
  new.rowcovs <- updated.rowandcol$rowcovs
  new.colcovs <- updated.rowandcol$colcovs
  # new.rowcovs <- rowcovs
  # new.colcovs <- colcovs
  #new.rowcovs <- updaterowcovs(obses, gamma.matrix, mus, colcovs)
  #new.colcovs <- updatecolcovs(obses, gamma.matrix, mus, rowcovs)
  new.pis <- rowSums(gamma.matrix) / sum(rowSums(gamma.matrix))
  return(list(mus = new.mus, rowcovs = new.rowcovs, colcovs = new.colcovs, pis = new.pis))
}
updaterowAndcolcovs <- function(obses, gamma.matrix, mus, tol = 1e-4, max.iter.loop = 10) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  dim.rowcov <- dim(mus)[2]
  dim.colcov <- dim(mus)[3]
  rowcovs <- array(0, c(dim.cluster, dim.rowcov, dim.rowcov))
  colcovs <- array(0, c(dim.cluster, dim.colcov, dim.colcov))
  for (i in 1:dim.cluster) {
    mu <- mus[i, , ]
    old.colcov <- diag(1, dim.colcov) 
    # old.rowcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(old.colcov + diag(tol, nrow(old.colcov))))
    # new.colcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(old.rowcov + diag(tol, nrow(old.rowcov))), FALSE)
    # new.rowcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(new.colcov + diag(tol, nrow(new.colcov))))
    old.rowcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(old.colcov))
    while (kappa(old.rowcov) > 1e14) {
      old.rowcov <- old.rowcov + diag(tol, nrow(old.rowcov))
      tol <- tol * 2
    }
    tol <- 1e-4
    new.colcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(old.rowcov), FALSE)
    while (kappa(new.colcov) > 1e14) {
      new.colcov <- new.colcov + diag(tol, nrow(new.colcov))
      tol <- tol * 2
    }
    tol <- 1e-4
    new.rowcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(new.colcov))
    while (kappa(new.rowcov) > 1e14) {
      new.rowcov <- new.rowcov + diag(tol, nrow(new.rowcov))
      tol <- tol * 2
    }
    tol <- 1e-4
    for (j in 1:max.iter.loop) {
      old.colcov <- new.colcov
      old.rowcov <- new.rowcov
      while (kappa(old.rowcov) > 1e14) {
        old.rowcov <- old.rowcov + diag(tol, nrow(old.rowcov))
        tol <- tol * 2
      }
      tol <- 1e-4
      new.colcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(old.rowcov), FALSE)
      while (kappa(new.colcov) > 1e14) {
        new.colcov <- new.colcov + diag(tol, nrow(new.colcov))
        tol <- tol * 2
      }
      tol <- 1e-4
      new.rowcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(new.colcov))
      # new.colcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(old.rowcov + diag(tol, nrow(old.rowcov))), FALSE)
      # new.rowcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(new.colcov+ diag(tol, nrow(new.colcov))))
      #     print (sum((old.colcov - new.colcov)^2))
      #     print (sum((old.rowcov - new.rowcov)^2))
      # if (sum((old.colcov - new.colcov)^2) > tol |
      #     sum((old.rowcov - new.rowcov)^2) > tol) {
      #   break
      # }
      if (sum((kronecker(old.colcov, old.rowcov) - kronecker(new.colcov,
                                                             new.rowcov))^2) < tol)
        break
    }
    rowcovs[i, , ] <- new.rowcov
    colcovs[i, , ] <- new.colcov
  }
  return(list(rowcovs = rowcovs, colcovs = colcovs))
}
addDiag <- function(data.matrix, tol) {
  if (kappa(data.matrix) <= 1e14) {
    return(data.matrix)
  } else {
    while (kappa(data.matrix) > 1e14) {
      data.matrix <- data.matrix + diag(tol, nrow(data.matrix))
      tol <- tol * 2
    }
    return(data.matrix)
  }
}
updaterowAndcolcovsTraceConstrain <- function(obses, gamma.matrix, mus, tol = 1e-4, max.iter.loop = 10) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  dim.rowcov <- dim(mus)[2]
  dim.colcov <- dim(mus)[3]
  rowcovs <- array(0, c(dim.cluster, dim.rowcov, dim.rowcov))
  colcovs <- array(0, c(dim.cluster, dim.colcov, dim.colcov))
  for (i in 1:dim.cluster) {
    mu <- mus[i, , ]
    old.colcov <- diag(1, dim.colcov) 
    old.rowcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(old.colcov))
    old.rowcov <- addDiag(old.rowcov, 1e-4)
    old.rowcov <- old.rowcov / tr(old.rowcov)
    new.colcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(old.rowcov), FALSE)
    new.colcov <- addDiag(new.colcov, 1e-4)
    new.rowcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(new.colcov))
    new.rowcov <- addDiag(new.rowcov, 1e-4)
    new.rowcov <- new.rowcov / tr(new.rowcov)
    for (j in 1:max.iter.loop) {
      old.colcov <- new.colcov
      old.rowcov <- new.rowcov
      old.rowcov <- addDiag(old.rowcov, 1e-4)
      old.rowcov <- old.rowcov / tr(old.rowcov)
      new.colcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(old.rowcov), FALSE)
      new.colcov <- addDiag(new.colcov, 1e-4)
      new.rowcov <- initializeUpdatedEM(obses, gamma.matrix[i, ], mu, solve(new.colcov))
      new.rowcov <- new.rowcov / tr(new.rowcov)
      if (sum((kronecker(old.colcov, old.rowcov) - kronecker(new.colcov,
                                                             new.rowcov))^2) < tol)
        break
    }
    rowcovs[i, , ] <- new.rowcov
    colcovs[i, , ] <- new.colcov
  }
  return(list(rowcovs = rowcovs, colcovs = colcovs))
}
updaterowcovs <- function(obses, gamma.matrix, mus, colcovs) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  dim.rowcov <- dim(mus)[2]
  dim.colcov <- dim(mus)[3]
  updated.rowcovs <- array(0, c(dim.cluster, dim.rowcov, dim.rowcov))
  for (i in 1:dim.cluster) {
    inv.colcov <- solve(colcovs[i, , ])
    for (j in 1:dim.sample) {
      updated.rowcovs[i, , ] <- updated.rowcovs[i, , ] + 
        gamma.matrix[i, j] * (obses[j, , ] - mus[i, , ]) %*% inv.colcov %*% t(obses[j, , ] - mus[i, , ])
    }
    updated.rowcovs[i, , ] <- updated.rowcovs[i, , ] / (sum(gamma.matrix[i, ]) * dim.colcov)
  }
  return(updated.rowcovs)
}
updatecolcovs <- function(obses, gamma.matrix, mus, rowcovs) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  dim.rowcov <- dim(mus)[2]
  dim.colcov <- dim(mus)[3]
  updated.colcovs <- array(0, c(dim.cluster, dim.colcov, dim.colcov))
  for (i in 1:dim.cluster) {
    inv.rowcov <- solve(rowcovs[i, , ])
    for (j in 1:dim.sample) {
      updated.colcovs[i, , ] <- updated.colcovs[i, , ] + 
        gamma.matrix[i, j] * t(obses[j, , ] - mus[i, , ]) %*% inv.rowcov %*% (obses[j, , ] - mus[i, , ])
    }
    updated.colcovs[i, , ] <- updated.colcovs[i, , ] / (sum(gamma.matrix[i, ]) * dim.rowcov)
  }
  return(updated.colcovs)
}
updateMus <- function(obses, gamma.matrix, mus, rowcovs, colcovs, pis, step.size = 0.4, lambda) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  updated.mus <- array(0, dim(mus))
  penalized.Q.old <- penalizedQ(obses, mus, rowcovs, colcovs, pis, gamma.matrix, lambda)
  penalized.Q.new <- penalized.Q.old - 9999
  while (penalized.Q.old > penalized.Q.new) {
    for (i in 1:dim.cluster) {
      inv.rowcov <- solve(rowcovs[i, , ])
      inv.colcov <- solve(colcovs[i, , ])
      weighted.sum <- matrix(0, nrow = dim(obses)[2], ncol = dim(obses)[3])
      for (j in 1:dim(obses)[1]) {
        weighted.sum <- weighted.sum + gamma.matrix[i, j] * obses[j, ,]
      }
      updated.mus[i, , ] <- mus[i, , ] + step.size * (inv.rowcov %*% (weighted.sum - 
                                                                        sum(gamma.matrix[i, ]) * mus[i, , ]) %*% 
                                                        inv.colcov - lambda * ((mus[i, , ] > 0) * 2 - 1))
    }
    penalized.Q.new <- penalizedQ(obses, updated.mus, rowcovs, colcovs, pis, gamma.matrix, lambda)
    step.size <- step.size / 2
  }
  return(updated.mus)
}
updateMusexact <- function(obses, gamma.matrix, mus, rowcovs, colcovs, pis, lambda) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  dim.row <- nrow(rowcovs[1, , ])
  dim.col <- nrow(colcovs[1, , ])
  updated.mus <- array(0, dim(mus))
  for (i in 1:dim.cluster) {
    mus.prime <- matrix(0, dim.row, dim.col)
    for (j in 1:dim.sample) {
      mus.prime <- mus.prime + gamma.matrix[i, j] * obses[j, , ]
    }
    mus.prime <- mus.prime / sum(gamma.matrix[i, ])
    temp <- (mus.prime) - (lambda/sum(gamma.matrix[i, ])) * rowcovs[i, , ] %*%
      ((mus[i, , ] > 0) * 2 - 1) %*% colcovs[i, , ]
    updated.mus[i, , ] <- temp
  }
  return(updated.mus)
}
updateMusLassoexact <- function(obses, gamma.matrix, mus, rowcovs, colcovs, pis, lambda) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  dim.row <- nrow(rowcovs[1, , ])
  dim.col <- nrow(colcovs[1, , ])
  updated.mus <- array(0, dim(mus))
  for (i in 1:dim.cluster) {
    mus.prime <- matrix(0, dim.row, dim.col)
    for (j in 1:dim.sample) {
      mus.prime <- mus.prime + gamma.matrix[i, j] * obses[j, , ]
    }
    mus.prime <- mus.prime / sum(gamma.matrix[i, ])
    lambda.second <-  (lambda/sum(gamma.matrix[i, ])) * rowcovs[i, , ] %*%
      matrix(1, dim.row, dim.col) %*% colcovs[i, , ]
    temp <- ((mus.prime > 0) * 2 - 1) * ((abs(mus.prime) - lambda.second) > 0) * (abs(mus.prime) - lambda.second)
    updated.mus[i, , ] <- temp
  }
  return(updated.mus)
}
updateMusl2exact <- function(obses, gamma.matrix, mus, rowcovs, colcovs, pis, lambda) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  dim.row <- nrow(rowcovs[1, , ])
  dim.col <- nrow(colcovs[1, , ])
  updated.mus <- array(0, dim(mus))
  for (i in 1:dim.cluster) {
    mus.prime <- matrix(0, dim.row, dim.col)
    for (j in 1:dim.sample) {
      mus.prime <- mus.prime + gamma.matrix[i, j] * obses[j, , ]
    }
    mus.prime <- mus.prime / sum(gamma.matrix[i, ])
    temp <- (mus.prime) - (2*lambda/sum(gamma.matrix[i, ])) * rowcovs[i, , ] %*%
      (mus[i, , ]) %*% colcovs[i, , ]
    updated.mus[i, , ] <- temp
  }
  return(updated.mus)
}
updateMusl2 <- function(obses, gamma.matrix, mus, rowcovs, colcovs, pis, step.size = 0.4, lambda) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  updated.mus <- array(0, dim(mus))
  penalized.Q.old <- penalizedQl2(obses, mus, rowcovs, colcovs, pis, gamma.matrix, lambda)
  penalized.Q.new <- penalized.Q.old - 9999
  while (penalized.Q.old > penalized.Q.new) {
    for (i in 1:dim.cluster) {
      inv.rowcov <- solve(rowcovs[i, , ])
      inv.colcov <- solve(colcovs[i, , ])
      weighted.sum <- matrix(0, nrow = dim(obses)[2], ncol = dim(obses)[3])
      for (j in 1:dim(obses)[1]) {
        weighted.sum <- weighted.sum + gamma.matrix[i, j] * obses[j, ,]
      }
      updated.mus[i, , ] <- mus[i, , ] + step.size * (inv.rowcov %*% (weighted.sum - 
                                                                        sum(gamma.matrix[i, ]) * mus[i, , ]) %*% 
                                                        inv.colcov - lambda * mus[i, , ])
    }
    penalized.Q.new <- penalizedQl2(obses, updated.mus, rowcovs, colcovs, pis, gamma.matrix, lambda)
    step.size <- step.size / 2
  }
  return(updated.mus)
}
updateMusnuclear <- function(obses, gamma.matrix, mus, rowcovs, colcovs, pis, step.size = 0.4, lambda) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  updated.mus <- array(0, dim(mus))
  penalized.Q.old <- penalizedQnuclear(obses, mus, rowcovs, colcovs, pis, gamma.matrix, lambda)
  penalized.Q.new <- penalized.Q.old - 9999
  while (penalized.Q.old > penalized.Q.new) {
    for (i in 1:dim.cluster) {
      svd.mus <- svd(mus[i, , ])
      inv.rowcov <- solve(rowcovs[i, , ])
      inv.colcov <- solve(colcovs[i, , ])
      weighted.sum <- matrix(0, nrow = dim(obses)[2], ncol = dim(obses)[3])
      for (j in 1:dim(obses)[1]) {
        weighted.sum <- weighted.sum + gamma.matrix[i, j] * obses[j, ,]
      }
      updated.mus[i, , ] <- mus[i, , ] + step.size * (inv.rowcov %*% (weighted.sum - 
                                                                        sum(gamma.matrix[i, ]) * mus[i, , ]) %*% 
                                                        inv.colcov - lambda * svd.mus$u %*% t(svd.mus$v))
    }
    penalized.Q.new <- penalizedQnuclear(obses, updated.mus, rowcovs, colcovs, pis, gamma.matrix, lambda)
    step.size <- step.size / 2
  }
  return(updated.mus)
}
updateMusnuclearexact <- function(obses, gamma.matrix, mus, rowcovs, colcovs, pis, lambda) {
  dim.cluster <- dim(mus)[1]
  dim.sample <- dim(obses)[1]
  dim.row <- nrow(rowcovs[1, , ])
  dim.col <- nrow(colcovs[1, , ])
  updated.mus <- array(0, dim(mus))
  for (i in 1:dim.cluster) {
    svd.mus <- svd(mus[i, , ])
    mus.prime <- matrix(0, dim.row, dim.col)
    for (j in 1:dim.sample) {
      mus.prime <- mus.prime + gamma.matrix[i, j] * obses[j, , ]
    }
    mus.prime <- mus.prime / sum(gamma.matrix[i, ])
    temp <- (mus.prime) - (lambda/sum(gamma.matrix[i, ])) * rowcovs[i, , ] %*%
      (svd.mus$u %*% t(svd.mus$v)) %*% colcovs[i, , ]
    updated.mus[i, , ] <- temp
  }
  return(updated.mus)
}
initializeEM <- function(obses, cluster.nums, shuffle = FALSE, diag = TRUE) {
  obses.vec <- t(apply(obses, 1, as.vector))
  kmean.result <- kmeans(scale(obses.vec), cluster.nums)
  if (shuffle == TRUE) {
    kmean.result$cluster = sample(1:cluster.nums, dim(obses)[1], replace = TRUE)
  }
  mus <- array(0, c(cluster.nums, dim(obses)[2], dim(obses)[3]))
  rowcovs <- array(0, c(cluster.nums, dim(obses)[2], dim(obses)[2]))
  colcovs <- array(0, c(cluster.nums, dim(obses)[3], dim(obses)[3]))
  for (i in 1:cluster.nums) {
    if (sum(kmean.result$cluster == i) == 1) {
      mus[i, , ] <- obses[kmean.result$cluster == i,,]
    } else {
      mus[i, , ] <- apply(obses[kmean.result$cluster == i,,], c(2, 3), mean)
    }
    obs <- obses[kmean.result$cluster == i, , ]
    if (diag == TRUE) {
      rowcovs[i, , ] <- diag(1, dim(obses)[2])
      colcovs[i, , ] <- diag(1, dim(obses)[3])
    } else {
      cov.result <- initializeMLEcovariance(obs, mus[i, , ])
      rowcovs[i, , ] <- cov.result$rowcov 
      colcovs[i, , ] <- cov.result$colcov 
    }
  }
  pis <- table(kmean.result$cluster) / dim(obses)[1]
  return (list(pis = pis, mus = mus, rowcovs = rowcovs, colcovs = colcovs, result = kmean.result$cluster))
}
initializeMLEcovariance <- function(obs, mu, tol = 1e-4) {
  old.colcov <- diag(1, dim(obs)[3]) 
  old.rowcov <- initializeUpdated(obs, mu, solve(old.colcov))
  new.colcov <- initializeUpdated(obs, mu, solve(old.rowcov), FALSE)
  new.rowcov <- initializeUpdated(obs, mu, solve(new.colcov))
  while (sum((old.colcov - new.colcov)^2) > tol | 
         sum((old.rowcov - new.rowcov)^2) > tol) {
    old.colcov <- new.colcov
    old.rowcov <- new.rowcov
    new.colcov <- initializeUpdated(obs, mu, solve(old.rowcov), FALSE)
    new.rowcov <- initializeUpdated(obs, mu, solve(new.colcov))
  }
  return (list(colcov = new.colcov, rowcov = new.rowcov))
}
initializeUpdated <- function(obs, mu, inv.matrix, rowcov = TRUE) {
  if (rowcov == TRUE) {
    updated.rowcovs <- matrix(0, dim(obs)[2], dim(obs)[2])
    for (i in 1:dim(obs)[1]) {
      updated.rowcovs <- updated.rowcovs + (obs[i, , ] - mu) %*% inv.matrix %*%
        t(obs[i, , ] - mu)
    }
    return(updated.rowcovs / (dim(obs)[1] * nrow(inv.matrix)))
  } else {
    updated.colcovs <- matrix(0, dim(obs)[3], dim(obs)[3])
    for (i in 1:dim(obs)[1]){
      updated.colcovs <- updated.colcovs + t(obs[i, , ] - mu) %*% inv.matrix %*%
        (obs[i, , ] - mu)
    }
    return(updated.colcovs / (dim(obs)[1] * nrow(inv.matrix)))
  }
}
initializeUpdatedEM <- function(obs, gamma.vector, mu, inv.matrix, rowcov = TRUE) {
  if (rowcov == TRUE) {
    updated.rowcovs <- matrix(0, dim(obs)[2], dim(obs)[2])
    for (i in 1:dim(obs)[1]) {
      updated.rowcovs <- updated.rowcovs + gamma.vector[i] * (obs[i, , ] - mu) %*% inv.matrix %*%
        t(obs[i, , ] - mu)
    }
    return(updated.rowcovs / (sum(gamma.vector) * nrow(inv.matrix)))
  } else {
    updated.colcovs <- matrix(0, dim(obs)[3], dim(obs)[3])
    for (i in 1:dim(obs)[1]){
      updated.colcovs <- updated.colcovs +  gamma.vector[i] * t(obs[i, , ] - mu) %*% inv.matrix %*%
        (obs[i, , ] - mu)
    }
    return(updated.colcovs / (sum(gamma.vector) * nrow(inv.matrix)))
  }
}
penalizedQ <- function(obses, mus, rowcovs, colcovs, pis, gamma.matrix, lambda) {
  updated.loggamma <- expectationStep(obses, mus, rowcovs, colcovs, pis, scale = FALSE)
  first.Q <- sum(gamma.matrix * updated.loggamma)
  second.Q <- lambda * sum(abs(as.vector(mus)))
  penalized.Q <- first.Q - second.Q
  return(penalized.Q)
}
penalizedQl2 <- function(obses, mus, rowcovs, colcovs, pis, gamma.matrix, lambda) {
  updated.loggamma <- expectationStep(obses, mus, rowcovs, colcovs, pis, scale = FALSE)
  first.Q <- sum(gamma.matrix * updated.loggamma)
  second.Q <- lambda * sum((as.vector(mus))^2)
  penalized.Q <- first.Q - second.Q
  return(penalized.Q)
}
penalizedQnuclear <- function(obses, mus, rowcovs, colcovs, pis, gamma.matrix, lambda) {
  updated.loggamma <- expectationStep(obses, mus, rowcovs, colcovs, pis, scale = FALSE)
  first.Q <- sum(gamma.matrix * updated.loggamma)
  second.Q <- lambda * sum(apply(mus, 1,function(x) {return(sum(sqrt(diag(t(x) %*% x))))}))
  penalized.Q <- first.Q - second.Q
  return(penalized.Q)
}
penalizedLikelihoodobs <- function(obses, mus, rowcovs, colcovs, pis, lambda, penalty, single.cluster = FALSE) {
  if (single.cluster == TRUE) {
    first.L <- 0 
    for (i in 1:dim(obses)[1]) {
      first.L <- first.L + likelihood(obses[i,,], mus[1,,], rowcovs[1,,], colcovs[1,,])
    }
    lambda <- lambda * dim(obses)[1]
  } else {
    dim.row <- dim(mus)[1]
    dim.col <- dim(obses)[1]
    updated.loggamma <- expectationStep(obses, mus, rowcovs, colcovs, pis, scale = FALSE)
    # rescale for each column
    max.vec <- apply(updated.loggamma, 2, max)
    updated.loggamma <- updated.loggamma - 
      matrix(max.vec, nrow = dim.row, ncol = dim.col, byrow = TRUE)
    updated.vec <- log(colSums(exp(updated.loggamma)))
    first.L <- sum(updated.vec + max.vec)
  }
  if (penalty == 'l2') {
    second.L <- lambda * sum((as.vector(mus))^2)
  } else if (penalty == 'lasso') {
    second.L <- lambda * sum(abs(as.vector(mus)))
  } else {
    second.L <- lambda * sum(apply(mus, 1,function(x) {return(sum(sqrt(diag(t(x) %*% x))))}))
  }
  penalized.L <- first.L - second.L
  return(penalized.L)
}
mainEM <- function(obses, lambda, cluster.nums, max.iter, penalty = 'lasso', tol = 1e-4, shuffle.ini = TRUE) {
  loss <- c()
  loss <- c(loss, -999)
  initial.resultkmean <- initializeEM(obses, cluster.nums)
  initial.result <- initializeEM(obses, cluster.nums, shuffle = TRUE)
  if (shuffle.ini == FALSE) {
    initial.result <- initial.resultkmean
  }
  old.pis <- initial.result$pis
  old.mus <- initial.result$mus 
  old.rowcovs <- initial.result$rowcovs 
  old.colcovs <- initial.result$colcovs 
  gamma.matrix.old <- expectationStep(obses, old.mus, old.rowcovs, old.colcovs, old.pis, scale = TRUE)
  for (i in 1:max.iter) {
    # M step
    updated.result <- maximizationStep(obses, gamma.matrix.old, old.mus, old.rowcovs, old.colcovs, old.pis, lambda, penalty)
    new.pis <- updated.result$pis
    new.mus <- updated.result$mus
    new.rowcovs <- updated.result$rowcovs
    new.colcovs <- updated.result$colcovs
    obsLikelihood <- (penalizedLikelihoodobs(obses, new.mus, new.rowcovs, new.colcovs, new.pis, lambda, penalty))
    loss <- c(loss, obsLikelihood)
    # E step
    gamma.matrix.new <- expectationStep(obses, new.mus, new.rowcovs, new.colcovs, new.pis, scale = TRUE)
    if (penalty == 'nuclear') {
      if (
        (sum((new.mus - old.mus)^2) < tol) | (sum(apply(gamma.matrix.new, 2, max)) > (nrow(gamma.matrix.new)-1e-6)) 
      ) {   
        break
      }
    } else {
      if (
        sum((new.mus - old.mus)^2) < tol
      ) {   
        break
      }}
    old.pis <- new.pis
    old.mus <- new.mus
    old.rowcovs <- new.rowcovs
    old.colcovs <- new.colcovs
    gamma.matrix.old <- gamma.matrix.new
  }
  return(list(mixedmodel = 
                list(loss = loss, mus = new.mus, 
                     pis = new.pis, rowcovs = new.rowcovs, 
                     colcovs = new.colcovs, results =  apply(gamma.matrix.new, 2, which.max),
                     initial.result = initial.result), 
              kmeans = initial.resultkmean))
}
CalculateCVL <- function(obses, lambda, cluster.nums, max.iter, 
                         penalty = 'lasso', tol = 1e-4, num.fold = 10) {
  CV.likelihood <- c()
  folds <- cut(seq(1,dim(obses)[1]), breaks=num.fold,labels=FALSE)
  for (j in unique(folds)) {
    single.index <- which(folds == j)
    obses.train <- obses[-single.index, , ]
    obses.test <- obses[single.index, , ]
    if (cluster.nums == 1) {
      result.param <- initializeEM(obses.train, 1, diag = FALSE)
      single.cluster <- TRUE
    } else {
      lambda <- lambda / (cluster.nums * dim(obses.train)[2])
      result <- mainEM(obses.train, lambda, cluster.nums, max.iter, penalty, tol, TRUE)
      result.param <- result$mixedmodel
      single.cluster <- FALSE}
    CV.likelihood <- c(CV.likelihood, 
                       penalizedLikelihoodobs(obses.test, result.param$mus, 
                                              result.param$rowcovs, result.param$colcovs, 
                                              result.param$pis, lambda, penalty, single.cluster = single.cluster))
  }
  return(mean(CV.likelihood[CV.likelihood > 0]))
}



cluster.nums <- 2
sampledim <- c(20, 10, 10)
normal.params <- generateParamsSim(cluster.nums, sampledim, scale = 1, diag = FALSE)
dim <- list(pi = rep(1, cluster.nums) / sum(rep(1, cluster.nums)), sampledim = sampledim, normal.params = normal.params)
lambda <- 1.5/(cluster.nums * sampledim[2])
penalty = 'lasso'

generate.signals <- generateMixMatrixnormal(dim)
signals <- generate.signals$signals
result <- mainEM(signals, max.iter = 10, cluster.nums = cluster.nums, lambda = lambda, penalty = penalty)
CV.likehood <- CalculateCVL(signals, lambda, cluster.nums, max.iter = 10, 
                            penalty = penalty, tol = 1e-4, num.fold = 3)






