#' Calculates the Exponential Weighted Moving Average volatilty
#'
#'@param prices A price series on which you want to compute the EWMA volatility measure 
#'@param lamda The exponential decay factor to be used as weights
#'@param UseLog A boolean to determine if we should take the log of the price series. The series must all be positive to set UseLog = TRUE otherwise set to FALSE to use the differences method. 
#'
#'@author Helena Ristov
#'
#'@export


ewmaVol <- function(prices, dPeriod, lambda, UseLog = TRUE, .combine = FALSE, skipDiff = FALSE){
  
  if(!is.xts(prices)){ stop("ewmaVol requires an xts prices file.") }
  
  incr        <- abs(as.numeric(difftime(index(prices[1]), index(prices[2]), units = "secs")))
  IncrPerDay  <- max(sapply(unique(as.Date(index(prices))), function(x){ length(which(as.Date(index(prices)) == x)) }))
  
  ## Set up stationary price series
  if(!skipDiff){
    if(length(which(c(-1, 0, 1) %in% sign(prices))) >= 2 || !(UseLog)){
      if(UseLog){ warning("Can not use log prices when there are differing signs in the price series.  Using differences instead.") }
      log.prices  <- diff(prices)[-1,]
    }else{
      log.prices  <- log(as.matrix(prices[-1,])/as.matrix(prices[-length(prices),]))  
    }
    UseLog <- FALSE
  }else{
    log.prices <- prices
  }
  
  ## Determine whether to combine all obs or create a rolling window
  if(.combine){
    WeekdayNum <- length(GetWeekDays(as.Date(first(index(prices))), as.Date(last(index(prices)))))
    f.index <- numeric()
    for(i in WeekdayNum:1){ f.index <- append(f.index, rep(WeekdayNum - i, IncrPerDay)) }
    f.index <- rev(f.index[1:nrow(log.prices)])
    
    lambda.decay  <- ((1 - lambda)*lambda^f.index) / sum((1 - lambda)*lambda^f.index)
    
    if(UseLog){
      ewma.vol <-  sqrt(apply(lambda.decay*(log.prices)^2, 2, sum))  
    }else{
      mean.log.prices <- apply(log.prices, 2, mean)
      ewma.vol <-  sqrt(apply(lambda.decay*(log.prices - mean.log.prices)^2, 2, sum))  
    }
    
    
  }else{
    f.index <- numeric()
    for(i in dPeriod:1){ f.index <- append(f.index, rep(dPeriod - i, IncrPerDay)) }
    f.index <- rev(f.index)    
    
    lambda.decay       <- ((1 - lambda)*lambda^f.index) / sum((1 - lambda)*lambda^f.index)
    ewma.vol           <- rollapply(log.prices, width = IncrPerDay*dPeriod, function(x){ sum((x-mean(x))^2 * lambda.decay) })
    ewma.vol           <- sqrt(ewma.vol)
  }
  
  colnames(ewma.vol) <- paste0(colnames(prices), '.ewma.vol')
  return(ewma.vol)
  
}


##prices <- c(29.66, 30.14, 29.85, 30.1, 30.95, 31.72, 31.3, 30.59, 29.06, 29.42, 28.76, 27.77, 29, 26.62, 25.91, 26.05, 26.51, 26.93, 26.26, 27.36, 27.41)

##prices <- as.xts(prices, order.by = Sys.Date()-1:21 )
##colnames(prices) <- 'price'
##dPeriod <- 3
