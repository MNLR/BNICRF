
sortStationsAncestral <- function(cbn, 
                                  first.alone = FALSE, 
                                  merge.first.block = FALSE){
  
  if (class(cbn) == "cbn"){ cbn <- cbn$BN } # Only nodes are required.
  
  ancestral.stations.id <-
    sapply(cbn$nodes[
      order(sapply(cbn$nodes, function(nn) length(nn$parents)))
      ], 
      function(nn) sapply(strsplit(nn$parents, split = ".", fixed = T),
                          function(ss) ss[[2]])
    )
  
  first.node <- list(
     strsplit(names(ancestral.stations.id)[1], split = "D.")[[1]][2]
    )
  names(first.node) <- first.node
  
  ancestral.stations.id <-
    mapply(function(c, p) c(p,c),  sapply(strsplit(names(ancestral.stations.id),
                                                   split = ".", fixed = T), 
                                          function(ss) ss[[2]]),
           ancestral.stations.id)
  
  ancestral.stations.id <- ancestral.stations.id[2:length(ancestral.stations.id)]
  
  
  ### Sorting:
  
  for (i in 2:length(ancestral.stations.id)){
    parents_ <- ancestral.stations.id[[i]][1:(length(ancestral.stations.id[[i]])-1)]
    explored <- unlist(ancestral.stations.id[1:(i-1)])
    
    j <- 1
    while(sum(parents_ %in%  explored) != length(parents_)){
      
      parents_ <- ancestral.stations.id[[i + j]][1:(length(ancestral.stations.id[[i + j]])-1)]
      explored <- unlist(ancestral.stations.id[1:(i-1)])
      
      j <- j+1  
    }
    if (j != 1) {
      aux <- ancestral.stations.id[[i]]
      ancestral.stations.id[[i]] <- ancestral.stations.id[[i+j-1]] # suma al final
      ancestral.stations.id[[i+j-1]] <- aux
    }
    
  }
  
  names(ancestral.stations.id) <- sapply(ancestral.stations.id, function(ss) ss[length(ss)])
  
  if (merge.first.block){
    remove <- c()
    for (i in 1:(length(ancestral.stations.id)-1)){
      if (sum(ancestral.stations.id[[i]] %in% ancestral.stations.id[[i+1]]) == length(ancestral.stations.id[[i]])) remove <- c(remove, i)
    }
    
    ancestral.stations.id <- ancestral.stations.id[-remove]
    
    aux <- ancestral.stations.id[[1]]
    
    inode <- 2
    while (inode < length(aux)){
      node <- aux[inode]
      parents <- sapply( strsplit(cbn$nodes[[ paste0("D.", node) ]]$parents, 
                                  split =  "D."), 
                         function(xx) xx[2] )
      
      if (sum(parents %in% aux[1:inode]) != length(parents)){
        aux <- c(aux, aux[inode])
        aux <- aux[-inode]
      } else { inode <- inode + 1 }
    }
    
    ancestral.stations.id[[1]] <- aux
  }
  
  if (first.alone) ancestral.stations.id <- c(first.node, ancestral.stations.id)
  
  return(ancestral.stations.id)
}

