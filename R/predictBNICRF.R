predictBNICRF <- function(BNICRF, x, 
                          non.informative.threshold = 3,
                          non.informative.p = "marginals"){
  # x must be a named list. Same names as training.
  
  nif <- c()
  ancestral.order <- list()
  for(iCRF in names(BNICRF)){
    pc <- attr(BNICRF[[iCRF]], "parents&child")
    ancestral.order[[iCRF]] <- pc
    if (iCRF == names(BNICRF)[1]){
      pr1 <- randomForestPredict(model = BNICRF[[iCRF]], 
                                 newdata = x[[paste("st", pc, sep = "_")]],
                                 bagging.function = mean)
    } else {
      xx <- do.call(cbind,
                    lapply(paste("st", pc, sep = "_"),
                           function(iiist) return(x[[iiist]])
                           )
                    )
      
      unique.columns <- 1
      if (!identical(xx[,2], xx[,1])) unique.columns <- c(unique.columns, 2)
      if (ncol(xx) >=3){
        for (ic in 3:ncol(xx)){
          identicS <-
            apply(xx[, 1:(ic-1) ], MARGIN = 2, FUN = function(cc){
              return( identical(cc, xx[, ic]) )
            })
          
          if (sum(identicS) ==  0) unique.columns <- c(unique.columns, ic)
        }
      }
      
      xx <- xx[, unique.columns]
      
      if (iCRF == names(BNICRF)[2]){
        pred.ap <-
          randomForestPredict(model = BNICRF[[iCRF]], 
                              newdata = xx
                              , method = "aposteriori",
                              non.informative.threshold = non.informative.threshold,
                              non.informative.p = non.informative.p)
        
        # Rename the columns to Location ID:
        aux.colnames <- colnames(pred.ap)
        for ( ivar in 1:length(pc) ){
          aux.colnames <- gsub(x = aux.colnames, 
                               pattern = paste0("Y", ivar), 
                               replacement = pc[ivar])
          
        }
        colnames(pred.ap) <- aux.colnames
        # E Rename the columns to Location ID:
        nif <- c(nif,
                 attr(pred.ap, "not.informed.values")) # appends not informed values.
        
        
        pred.ap[,1] <- pr1
        
        pr1 <- NULL
      } else {
        
        aux <-
          randomForestPredict(model = BNICRF[[iCRF]],
                              newdata = xx
                              , method = "aposteriori",
                              non.informative.threshold = non.informative.threshold,
                              non.informative.p = non.informative.p)
        
        # Rename the columns to Location ID:
        aux.colnames <- colnames(aux)
        for ( ivar in 1:length(pc) ){
          aux.colnames <- gsub(x = aux.colnames, 
                               pattern = paste0("Y", ivar), 
                               replacement = pc[ivar])
          
        }
        colnames(aux) <- aux.colnames
        # E Rename the columns to Location ID:
        
        
        
        aux.idx <- ncol(aux) - seq(from = 2^(length(pc)-1)-1,
                                   to = 0,
                                   by = -1)
        nif <- c(nif, attr(aux, "not.informed.values")[aux.idx]) # appends not informed values.
        
        aux <- aux[ ,
                    aux.idx]
        pred.ap <- cbind(pred.ap, aux)
      }
    }
  }
  
  attr(pred.ap, "not.informed.values") <- nif
  attr(pred.ap, "ancestral.order") <- ancestral.order
  attr(pred.ap, "y_names") <- attr(BNICRF, "y_names")
  
  return(pred.ap)
}
  