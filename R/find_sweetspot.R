# library(doParallel)
# library(foreach)
# library(parallel)
# library(RcppRoll)

# This function computes all differences in means, for a fixed window size k.
# It returns the value of the largest statistic, and the window to which it corresponds
one.statistic <- function(k, original.sequence, benefits.cumsum, n, overall.mean, stat="mean"){
  result  <- ((benefits.cumsum[((2+k):n)-1] - benefits.cumsum[(2:(n-k))-1]) - k*overall.mean)
  best.ix <- which.max(result)
  
  return(
    c(best.ix, best.ix+k-1, max(result, na.rm=T), k)
  )
}


compute.all.statistics <- function(n, original.sequence, benefits.cumsum, start, end, maxonly=F){
  result <- do.call(rbind,  lapply(start:(end-1), function(k) one.statistic(k, original.sequence, benefits.cumsum, n, benefits.cumsum[n]/n)))
  if(maxonly) return( max(result[, 3]) )
  
  best   <- result[which.max(result[, 3]), ]

  return(
    list(start.index  = best[1],
         end.index    = best[2],
         mean.inside  = sum(original.sequence[best[1]:best[2]])/best[4], #(benefits.cumsum[best[2]] - benefits.cumsum[best[1]-1])/best[4],
         mean.outside = sum(original.sequence[-(best[1]:best[2])])/(n-best[4]), #(benefits.cumsum[n] - benefits.cumsum[best[2]] + benefits.cumsum[best[1]-1])/(n - best[4]),
         statistic    = best[3]
    )
  )
}

find.sweetspot <- function(observed.benefit, 
                           ntrials.significance=1000, ntrials.bias=1000,
                           min.size.fraction = 1/20, max.size.fraction = 1/2, ncores=NULL){
  
  # https://stackoverflow.com/questions/50571325/r-cran-check-fail-when-using-parallel-functions
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cores <- 2L
  } else {
    cores <- if(is.null(ncores)){ detectCores()-1 } else { ncores }
  }
  registerDoParallel(cores = cores)
  
  n                <- length(observed.benefit)
  search.start     <- max(4, floor(min.size.fraction*n))
  search.end       <- floor(max.size.fraction*n)

  # Compute the statistic for all possible windows on our observed data, and choose the largest:
  best.window      <- compute.all.statistics(n, observed.benefit, cumsum(observed.benefit), search.start, search.end, maxonly=F)
  
  # Debias:
  values.inside  <- observed.benefit[best.window$start.index:best.window$end.index]
  values.outside <- observed.benefit[-(best.window$start.index:best.window$end.index)]
  
  pvalue <- foreach(i=1:ntrials.significance, .combine = sum, .inorder=F, .multicombine=T, .maxcombine=ntrials.significance) %dopar% {
    sampled <- sample(observed.benefit)
    compute.all.statistics(n, sampled, cumsum(sampled), search.start, search.end, maxonly=T) >= best.window$statistic
  }
  pvalue <- pvalue/ntrials.significance
 
  debias.boot <- foreach(i=1:ntrials.bias, .combine = rbind, .inorder=F, .multicombine=T, .maxcombine=ntrials.bias) %dopar% {
    sampled <- c(sample(values.outside, best.window$start.index - 1, replace=T),
                        sample(values.inside,  best.window$end.index - best.window$start.index + 1, replace=T),
                        sample(values.outside, max(0, n - best.window$end.index), replace=T))
    compute.all.statistics(n, sampled, cumsum(sampled),  search.start, search.end, maxonly=F)
  }


  # Return:
  list(start.index  = best.window$start.index,
       end.index    = best.window$end.index,
       statistic    = best.window$statistic,
       mean.inside  = best.window$mean.inside,
       mean.outside = best.window$mean.outside,

       mean.inside.debiased  = 2 * best.window$mean.inside  - 1/ntrials.bias * sum(unlist(debias.boot[, "mean.inside"])),
       mean.outside.debiased = 2 * best.window$mean.outside - 1/ntrials.bias * sum(unlist(debias.boot[, "mean.outside"])),

       means         = unname(unlist(debias.boot[, "mean.inside"])),
       means.outside = unname(unlist(debias.boot[, "mean.outside"])),
       start.indices = unname(unlist(debias.boot[, "start.index"])),
       end.indices   = unname(unlist(debias.boot[, "end.index"])),
       bias.inside   = best.window$mean.inside - 1/ntrials.bias * sum(unlist(debias.boot[, "mean.inside"])),
       bias.outside  = best.window$mean.outside- 1/ntrials.bias * sum(unlist(debias.boot[, "mean.outside"])),

       p.value      = pvalue
  )
}
