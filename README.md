# Regularized_Matrix_Clustering
This repository is to implement the regularized mixture matrix clustering method.

## An example

1. Source library functions

```source(library_functions.R)```

2. Generate sample signals

```cluster.nums <- 2
sampledim <- c(20, 10, 10)
cluster.nums <- 2
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
## Real Data Access
1. For Rat Stroke Experiment, source data is available from [here](https://drive.google.com/drive/u/1/folders/1coYfDCnKL04ah2e89MkibK8bOR74iv1s).
2. For Odor Memory Experiment, source data is available from [here](https://drive.google.com/drive/u/1/folders/1xEznXAZ9-eGUNX4XpKrZx0DrKWC-oohi).
