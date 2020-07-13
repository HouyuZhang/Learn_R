library(parallel)
library(MASS)
RNGkind("L'Ecuyer-CMRG")
mc.cores <- detectCores()
results <- mclapply(rep(25, 4),
                    function(nstart) kmeans(Boston, 4, nstart=nstart),
                    mc.cores=mc.cores)
i <- sapply(results, function(result) result$tot.withinss)
result <- results[[which.min(i)]]























