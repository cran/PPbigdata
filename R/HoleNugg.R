##Function to calculate Hole index for nugget projection

HoleNugg<- function(nuggproj,weight){
  nuggproj = as.matrix(nuggproj)
  n <- nrow(nuggproj)
  d <- ncol(nuggproj)

  if (length(weight) != n){
    stop("The length of weights should be the same as the number of nuggets. ")
  }

  if(sum(weight <=0) > 0){
    stop("The nugget weights should be positive values.")
  }

  sum = exp(-0.5*diag(nuggproj%*%t(nuggproj)))%*%matrix(weight,n,1)
  num <- 1 - sum[1,1]/sum(weight)
  den <- 1 - exp(-d / 2)
  num / den
}
