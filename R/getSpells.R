source("R/getConsecutiveDatesIndex.R")

getSpells <- function(Data, dates){
  idxs <- getConsecutiveDatesIndex(dates)
  cuts <- which( (idxs[2:nrow(idxs),"id_T-1"] != idxs[1:(nrow(idxs)-1),"id_T"]))
  cuts <- c(0, cuts, nrow(idxs))
  
  consecutive.indices <-
    lapply(1:(length(cuts)-1),
           function(icut){
             return(
               idxs[cuts[icut]+1,"id_T-1"]:idxs[cuts[icut + 1],"id_T"] 
             )
           })
  
  return(
    apply(Data, MARGIN = 2,
          function(obs){
            spellS <- matrix(nrow = 3, ncol = 0)
            rownames(spellS) <- c("Occurrence", "Days", "Year")
            
            
            for (cut in consecutive.indices){
              consecutive.years <- sapply(strsplit(dates[cut], split = "-"), 
                                          function(dd) dd[1])
              
              consecutive <- obs[cut]  
              previous.day <- consecutive[1]
              
              spells <- matrix(nrow = 3, data = c(previous.day, 1, as.numeric(consecutive.years[1])) )
              rownames(spells) <- c("Occurrence", "Days", "Year")
              
              for (iday in 2:length(consecutive)){
                day <- consecutive[iday]
                year <- consecutive.years[iday]
                if (day == previous.day){ spells["Days", ncol(spells)] <- 1 + spells["Days", ncol(spells)] }
                else {
                  spells <- cbind(spells, 
                                  matrix(data = c(day, 1, as.numeric(year)), nrow = 3, ncol = 1)
                  ) 
                }
                previous.day <- day
              }
              
              spellS <- cbind(spellS, spells)
            }
            
            return(t(spellS))
          }, simplify = FALSE)
  )
  
}
