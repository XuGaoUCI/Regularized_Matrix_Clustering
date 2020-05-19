# Regularized_Matrix_Clustering
This repository is to implement the regularized mixture matrix clustering method.

## An example

1. Source library functions

```source(library_functions.R)```

2. Generate sample signals

```cluster.nums <- 2
sampledim <- c(20, 10, 10)
normal.params <- generateParamsSim(cluster.nums, sampledim, scale = 1, diag = FALSE)
dim <- list(pi = rep(1, cluster.nums) / sum(rep(1, cluster.nums)), sampledim = sampledim, normal.params = normal.params)
generate.signals <- generateMixMatrixnormal(dim)
```

3. Calculate the clustering result 

```
lambda <- 1/(cluster.nums * sampledim[2])
penalty = 'lasso'
signals <- generate.signals$signals
result <- mainEM(signals, max.iter = 10, cluster.nums = cluster.nums, lambda = lambda, penalty = penalty)
```

4. Calcualte the CVPL values

```
CV.likehood <- CalculateCVL(signals, lambda, cluster.nums, max.iter = 10, penalty = penalty, tol = 1e-4, num.fold = 3)
```
