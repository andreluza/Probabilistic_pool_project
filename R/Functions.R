
## logic function
`%notin%` <- function(x,y) !(x %in% y)

# the null model to be applied across communities
sampling_prevalence_NM <- function (pool, local,nsamples){## 
  sampling <- lapply (as.list(seq(1, nrow (pool))), function (k) ## over all sites
    lapply (as.list(rep(1,nsamples)), function (i) ## replicate 1 x niter
      replicate(i,sample (x=pool[k,],size=local[k],prob=pool[k,], replace=F)))) ## do the prevalence sampling
  # obtain only the sampled probabilities
  result_list <- vector ("list")
  result_list$probabilities <- lapply (sampling, function (i) do.call (cbind, i))
  # obtain the random list of sampled species
  result_list$species <- lapply (sampling, function (k) do.call (cbind, lapply (k, function (i) dimnames(i)[[2]])))
  # return results
  return (result_list)
}



# get legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
