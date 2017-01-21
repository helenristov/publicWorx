#' Calculates the population stability index for two sample variables
#'
#'@param data The dataset that contains the two variables for comparison
#'@param var  The development variable or score used to set the distribution 
#'@param valvar The validation variable to compare against the development variable
#'@param breaks The number of bins to break the development variable for stability testing. Default is decile distribution.
#'@author Helena Ristov
#'@examples \dontrun{
#' x1      <- runif(1000, 1, 10)
#' x1_val  <- runif(1000, 1, 10)
#' 
#' random <- sample(c(0,1),1)
#' for(i in 1:999){
#' randomnext <- sample(c(0,1),1)
#' random     <- c(random, randomnext)
#'  
#' }
#'
#' test <- cbind(x1, x1_val, random)
#' calPSI(test, 'x1', 'x1_val')
#' }

#'@export

calcPSI  <- function(data, var, valvar, breaks = 10){
  
## get the approprate variable from their respective names in the table 
  if(class(data) == 'matrix'){
    var     <- data[,var]
    valvar  <- data[,valvar]
    ##perf    <- data[,perf]
  
  }else if('data.table' %in% class(data)){
    var    <- data[,var, with=FALSE]
    valvar <- data[,valvar, with=FALSE]
    ##perf   <- data[,perf, with=FALSE]
  
  }else if(class(data) == 'data.frame'){
    var    <- eval(parse(text=paste0( "data$",var)))
    valvar <- eval(parse(text=paste0( "data$",valvar)))
    ##perf   <- eval(parse(text=paste0( "data$",perf)))
  }

  data <- as.data.frame(cbind(var,valvar))
  
  ##compute the breaks
  Buckets <- c(0:breaks)/breaks
  quantiles_var <- quantile(x=data$var, probs=Buckets, na.rm =TRUE)
  quantile_data <- data.frame(id=names(quantiles_var), values=quantiles_var)
  
  table <- data.frame(Score.Range = c(1:breaks), Expected = NA, Actual = NA, Index=NA)
  
  for(i in 1:(length(quantile_data$values)-1)){
   
    if(i==1){
      
      range       <- paste0('Inf -- ', round(quantile_data$values[i+1],3))
      var_dist    <- which(data$var < quantile_data$value[i+1])
      valvar_dist <- which(data$valvar < quantile_data$value[i+1])
    
    }else if(i==length(quantile_data$values)-1){
      range       <- paste0(round(quantile_data$values[i],3), " -- Max")
      var_dist    <- which(data$var >= quantile_data$value[i])
      valvar_dist <- which(data$valvar >= quantile_data$value[i])
    }else{
      range <- paste0(round(quantile_data$values[i],3), "--", round(quantile_data$values[i+1],3))
      var_dist    <- which(data$var >= quantile_data$value[i] & data$var < quantile_data$value[i+1]  )
      valvar_dist <- which(data$valvar >= quantile_data$value[i] & data$valvar < quantile_data$value[i+1]  )
    }
    
    table$Score.Range[i]   <- range
    table$Expected[i] <- length(var_dist)
    table$Actual[i]   <- length(valvar_dist)
  }   

  table$Index <-  (table$Actual/sum(table$Actual) - table$Expected/sum(table$Expected)) * log((table$Actual/sum(table$Actual))/(table$Expected/sum(table$Expected)))
  
  PSI <- sum(table$Index)
    
  if(PSI < .1){
    print("No significant population differences")
  }else if(PSI < .25){
    print("Moderate population differences")
  }else{
    print("Substantial population differences")
  }
  
  return(list(table=table, PSI=PSI))
}    