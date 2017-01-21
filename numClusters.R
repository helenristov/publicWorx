#' Calculates the optimal number of clusters using a cluster R^2 equivalent metric, namely, the between sum of squares over the total sum of squares metric.
#'
#'@param Data  The data that you wish to examine containing the cluster dimensions
#'@param clusterIter  The number of clusters centers for examination. A setting of 10 will run 2 through 10 cluster functions and return the best choice. 
#'@param threshold The difference in the cluster R^2 must meet in order to include the next cluster iteration. 
#'@author Helena Ristov
#'
#'@export

numClusters <- function(Data, clusterIter = 10, threshold = .01){
  
  for(i in 2:clusterIter){
    cData        <- kmeans(Data, centers = i, nstart= 100)  
    
    if(i == 2){ 
      
      explainedVar <- cData$betweenss/cData$totss
    }else{
      nextExpVar   <- cData$betweenss/cData$totss
      explainedVar <- rbind(explainedVar,nextExpVar)
    }
  }
  
  numClusters <- first(which(diff(explainedVar) < threshold))+1
  xaxis       <- seq(2,clusterIter, by=1)
  
  p <-  qplot(x = xaxis, y = explainedVar, geom = c("point", "line"), xlab = 'Clusters', main = 'Explained Variance by Number of Cluster Centers')
  p <- p + geom_line(size=0.8, colour='blue') + geom_point(colour="black", size = 4.5) 
  plot(p)
  
  return(numClusters)
} 