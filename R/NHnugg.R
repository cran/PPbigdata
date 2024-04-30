##Function to calculate PP index for 1d/2d nugget projection

NHnugg<- function(nuggproj, weight, scale, bandwidth = NULL, gridn = 300,lims = NULL, gridnAd = TRUE) {

  nuggproj = as.matrix(nuggproj)
  dim <- ncol(nuggproj)
  if(!dim %in% 1:2){
    stop("The dimension of projected values should be 1 or 2.")
  }

  n = nrow(nuggproj)
  weight = as.vector(weight)
  scale = as.vector(scale)

  if (length(weight) != n){
    stop("The length of weights should be the same as the number of nuggets. ")
  }

  if(sum(weight <=0) > 0){
    stop("The nugget weights should be positive values.")
  }

  if (length(scale) != n){
    stop("The length of scales should be the same as the number of nuggets. ")
  }

  if(sum(scale <0) > 0){
    stop("The nugget scales should be non-negative values.")
  }

  if(sum(scale == 0) >0){
    scale[scale == 0] = 1e-8
  }

  #set defalt limits for density estimation
  if(is.null(lims)){
    lims = if(dim == 1)range(nuggproj[,1]) else c(range(nuggproj[,1]), range(nuggproj[,2]))
  }

  #estimate density based on nuggets
  estw = nuggKDE(nuggproj,weight, scale, h = bandwidth, gridn,lims,gridnAd)

  gridnx = length(estw$x)
  if(dim == 2)gridny = length(estw$y)

  #cell size for numerical integral
  cell_size = if(dim == 1)diff(lims[1:2]) / gridnx else (diff(lims[1:2]) / gridnx) * (diff(lims[3:4]) / gridny)

  if(dim == 1){
    phiw = dnorm(estw$x)
  }else{
    x0w = rep(estw$x,gridny)
    y0w = rep(estw$y,rep(gridnx,gridny))
    phiw = dmvnorm(cbind(x0w,y0w),rep(0,2),diag(2))
    dim(phiw) = c(gridnx,gridny)
  }

  #calculate index value
  funw = (estw$z-phiw)^2*phiw
  ## get numerical integral by summation:
  sum(funw) * cell_size
}
