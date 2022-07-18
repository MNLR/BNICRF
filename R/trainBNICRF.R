# LIBRARIES AND SOURCES
if ( !require(BNWeatherGen, quietly = T) ){ library(bnlearn) }
library(RandomForestDist)
library(future.apply)
source("R/sortStationsAncestral.R")
# LIBRARIES AND SOURCES


# INTERNAL
buildCRF <- function(stid, ancestral.order, x, y, ntree, minbucket, mtry, p){
  if (stid == names(ancestral.order)[1]){
    full.data <- x[[paste("st", ancestral.order[[stid]], sep = "_")]]
    
    tbr <- randomForestTrain(x = full.data,
                             y = as.numeric(y$Data[ , ancestral.order[[stid]]]), 
                             method = "binaryCrossEntropy",
                             ntree = ntree, 
                             minbucket = minbucket, 
                             progress.bar = F,
                             parallel.plan = NA)
    
  } else {
    full.data <- lapply(paste("st", ancestral.order[[stid]], sep = "_"),
                        function(iiist) return(x[[iiist]]))
    
    x.train <- do.call(cbind, full.data)
    
    y.train <- matrix(as.numeric(y$Data[, ancestral.order[[stid]]]),
                      nrow = nrow(y$Data[, ancestral.order[[stid]]]),
                      ncol = ncol(y$Data[, ancestral.order[[stid]]]))
    #remove duplicates in x.train
    
    unique.columns <- 1
    if (!identical(x.train[,2], x.train[,1])) unique.columns <- c(unique.columns, 2)
    for (ic in 3:ncol(x.train)){
      identicS <-
        apply(x.train[, 1:(ic-1) ], MARGIN = 2, FUN = function(cc){
          return( identical(cc, x.train[, ic]) )
        })
      
      if (sum(identicS) ==  0) unique.columns <- c(unique.columns, ic)
    }
    
    x.train <- x.train[, unique.columns]
    
    if (is.null(mtry)){ 
      mtry <- floor(sqrt(ncol(x.train)))
    } else { mtry <- floor(mtry(ncol(x.train))) }
    
    tbr <- randomForestTrain(x = x.train,
                             y = y.train,
                             minbucket = minbucket,
                             minsplit = 2*minbucket,
                             mtry = mtry,
                             oob.prunning = F,
                             ntree = ntree,
                             method = "binaryMultiEntropyCond",
                             parallel.plan = NA,
                             progress.bar = F,
                             resample = T,
                             replace = F,
                             sampsize = floor(0.632*nrow(x.train)), 
                             undersample.binary = F,
    )
  }
  
  attr(tbr, "parents&child") <- ancestral.order[[stid]]
  
  p(message = sprintf("Location %s (of %g locations)", stid, length(ancestral.order)))
  
  return(tbr)
}









trainBNICRF <- function(y, x, 
                        max.parents = 3,
                        tabu.list.size = 10^4,
                        ntree = 100,
                        minbucket = 15,
                        mtry = NULL,
                        progress.bar = TRUE,
                        quiet = FALSE,
                        parallel.plan,
                        workers = parallel::detectCores(),
                        return.informativeBN = FALSE,
                        informativeBN = NULL
                        ){
  # y can be either a list with $Data, $xyCoords and $Metadata$station_id
  # $Data column names will be used
  # $
  
  
  
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
  
  
      ## x is a list of length ncol(y), containing the predictors for each marginal variable.
  
  
  if (!is.list(y)){ 
    if (!is.matrix(y)){stop("y must be a list or a matrix.")} 
    else{
      aux <- y 
      if (is.null(colnames(y))) colnames(y) <- paste0("loc_", 1:ncol(y))
      y <- list(Data = y, 
                xyCoords = data.frame(x = 1:ncol(y), y = 1:ncol(y)),
                Metadata = list(station_id = colnames(y) )
                )
    }
  }
  
  aux <- apply(y$Data, MARGIN = 2, function(cl) as.factor(as.character(cl)))
  colnames(aux) <- colnames(y$Data)
  y$Data <- aux
  aux <- NULL
  
  
  if (is.null(informativeBN)){
    if (require(BNWeatherGen, quietly = TRUE)){
    informativeBN <- buildDescriptive(y = y,
                               structure.learning.algorithm = "tabu",
                               structure.learning.args.list = list(maxp = max.parents,
                                                                   tabu = tabu.list.size)
                               )
    } else { # Use BNLearn instead
      if (!quiet) print("Training Bayesian Network. This may take a while...")
      
      aux <- as.data.frame(lapply(1:ncol(y$Data), FUN = function(icol) as.factor(y$Data[, icol])))
      names(aux) <- paste0("D.", colnames(y$Data))
        
      informativeBN <- tabu(x = aux, tabu = tabu.list.size, maxp = max.parents)
      aux <- NULL
    }
  }

  ancestral.order <- sortStationsAncestral(informativeBN,
                                           merge.first.block = F, 
                                           first.alone = T)
  
  if (!return.informativeBN) informativeBN <- NULL
  
  if (!quiet) print("Training Conditional Random Forests...")
  
  progressorids <- names(ancestral.order)
  with_progress(enable = progress.bar, expr = {
    p <- progressor(along = progressorids)
    if (lapply.opt == "future_lapply") {
      BNICRF <- 
        future_lapply(future.seed = TRUE, X = progressorids, FUN = buildCRF,
                      ancestral.order = ancestral.order, 
                      x = x,
                      y = y, 
                      ntree = ntree,
                      minbucket = minbucket,
                      mtry = mtry,
                      p = p)
    } else {
      BNICRF <- 
        lapply( X = progressorids, FUN = buildCRF,
                ancestral.order = ancestral.order, 
                x = x,
                y = y, 
                ntree = ntree,
                minbucket = minbucket,
                mtry = mtry,
                p = p)
    }
  })
  
  names(BNICRF) <- paste0("loc",names(ancestral.order))
  attr(BNICRF, "y_names") <- colnames(y$Data)
  if (return.informativeBN) attr(BNICRF, "informativeBN") <- informativeBN
  if (!quiet) print("Done.")

  return(BNICRF)
}