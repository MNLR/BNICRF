
getConsecutiveDatesIndex <- function(dates){
  previous.present <- c(NA, diff(as.Date(dates))==1)
  id.t <- which(previous.present)
  id.t_1 <- which(previous.present) - 1
  dates_consecutive <- cbind.data.frame(as.Date(dates[id.t_1]),
                                        as.Date(dates[id.t]),
                                        id.t_1,
                                        id.t)
  rownames(dates_consecutive) <- id.t
  colnames(dates_consecutive) <- c("date_T-1", "date_T", "id_T-1", "id_T")
  return(dates_consecutive)
}