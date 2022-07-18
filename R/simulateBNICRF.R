simulate1bnicrf <- function(nsim, 
                            n,
                            ancestral.order.withoutfirst,
                            prediction,
                            p){
  for (nodeindex in 1:length(ancestral.order.withoutfirst)){
    if (nodeindex == 1){
      common.idx <- sum(2^seq(0, to = length(ancestral.order.withoutfirst[[1]])-1, by = 1))
      prmat <- prediction[, 1:common.idx]
      class(prmat) <- "RandomForestDist.prediction.simulable"
      sim <- randomForestSimulate(prmat, distr = "conditionalBernoulli", n = 1)[,,1]
      colnames(sim) <- ancestral.order.withoutfirst[[1]]
    } else {
      current.nodes <- ancestral.order.withoutfirst[[nodeindex]][1:(length(ancestral.order.withoutfirst[[nodeindex]])-1)]
      ids.cond <- match(current.nodes, colnames(sim))
      aux.sim <- sim[, ids.cond]
      
      prmat.idx <- seq((common.idx + 1), by = 1, length.out = 2^(length(ancestral.order.withoutfirst[[nodeindex]])-1))
      common.idx <- max(prmat.idx)
      
      aux.prmat <- prediction[, prmat.idx]
      
      cn <- colnames(sim)
      sim <- cbind(sim,
                   sapply(X = 1:nrow(aux.sim), FUN = function(ipr){
                     ss <- aux.sim[ipr,]
                     p <- aux.prmat[ipr, sum( sapply(1:length(ss), function(ik) ss[ik]*2^(ik-1) ) )+1]
                     sample(c(0,1), size = 1, prob = c(1-p, p))
                   })
      )
      
      colnames(sim) <- c(cn, ancestral.order.withoutfirst[[nodeindex]][length(ancestral.order.withoutfirst[[nodeindex]])])
    }
  }
  
  p(message = sprintf("Simulation %g/%g", nsim, n))
  
  return( sim[, match(x = attr(prediction, "y_names"), table = colnames(sim))] )
}




simulateBNICRF <- function(prediction,
                           n = 1,
                           parallel.plan,
                           workers = parallel::detectCores(),
                           progress.bar = TRUE){
  
  if (requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)
  ) {
    lapply.opt <- "future_lapply"
    
    if (missing(parallel.plan) || is.null(parallel.plan)){
      parallel.plan <- future::plan()
    }
    else if (is.character(parallel.plan) && parallel.plan == "auto") {
      o.plan <- future::plan()
      parallel.plan <- future::plan(future::multisession, workers = parallel::detectCores())
      on.exit(future::plan(o.plan), add = TRUE)
    }
    else {
      if (!is.list(parallel.plan) &&
          !is.function(parallel.plan) &&
          is.na(parallel.plan)) lapply.opt <- "lapply"
      else{
        o.plan <- future::plan()
        if (missing(workers)) future::plan(parallel.plan)
        else future::plan(parallel.plan,
                          workers = if (is.null(workers) || workers == 0) (availableCores()) else workers
        )
        on.exit(future::plan(o.plan), add = TRUE)
      }
    }
  } else { # package future or future.apply not available
    print("Future or future.apply packages not available. Parallelization disabled.")
    lapply.opt <- "lapply"
  }
  
  
  ancestral.order <- attr(prediction, "ancestral.order")
  
  
  ancestral.order.withoutfirst <- ancestral.order[2:length(ancestral.order)]
  
  progressorids <- 1:n
  with_progress(enable = progress.bar, expr = {
    p <- progressor(along = progressorids)    
    if (lapply.opt == "future_lapply") {
      simS <-
        future_lapply( future.seed = T,
                       X = progressorids, FUN = simulate1bnicrf,
                       n = n,
                       ancestral.order.withoutfirst = ancestral.order.withoutfirst,
                       prediction = prediction,
                       p = p
                       )
    } else {
      simS <-
        lapply(X = progressorids, FUN = simulate1bnicrf,
               n = n,
               ancestral.order.withoutfirst = ancestral.order.withoutfirst,
               prediction = prediction,
               p = p
        )
    }
  })
  
  return(simS)

}