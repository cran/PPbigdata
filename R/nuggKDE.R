##Function to estimate density function of projected data based on the data nuggets

nuggKDE <- function(nuggproj, weight, scale, h = NULL, gridn = 300,lims = NULL, gridnAd = TRUE){

  #check arguments
  nuggproj = as.matrix(nuggproj)
  dim <- ncol(nuggproj)
  weight = as.vector(weight)
  scale = as.vector(scale)
  gridn = as.vector(gridn)


  if(!dim %in% 1:2){
    stop("The dimension of projected values should be 1 or 2.")
  }

  if (length(weight) != nrow(nuggproj)){
    stop("The length of weights should be the same as the number of nuggets. ")
  }

  if(sum(weight <=0) > 0){
    stop("The nugget weights should be positive values.")
  }

  if (length(scale) != nrow(nuggproj)){
    stop("The length of weights should be the same as the number of nuggets. ")
  }

  if(!is.numeric(gridn) | sum(gridn <= 0) >= 1){
    stop("The number of grid points in each direction should be positive.")
  }

  if(!length(gridn) %in% 1:2){
    stop("The number of grid points should be a scalar or a length-2 integer vector.")
  }


  x = if(dim == 1)as.vector(nuggproj) else nuggproj[,1]
  if(dim ==2)y = nuggproj[,2]


  nx = if(dim == 1)length(nuggproj) else nrow(nuggproj)

  if(is.null(lims)){
    lims = if(dim == 1)range(x) else c(range(x), range(y))
  }else if((dim == 1 & length(lims) != 2) | (dim == 2 & length(lims) != 4)){
    stop("The limits should have upper and lower limits for each projection direction.")
  }

  if (!is.null(lims) & any(!is.finite(lims)))
    stop("only finite values are allowed in 'lims'")

  w = weight

  h <- if (is.null(h) & dim == 1)bandwidth.nrd(rep(x,w))
  else if(is.null(h) & dim == 2)c(bandwidth.nrd(rep(x,w)), bandwidth.nrd(rep(y,w)))
  else if(dim == 2 & length(h) == 1)rep(h, length.out = 2L)
  h <- h/4

  hx <- scale + h[1]
  if(dim ==2)hy <- scale + h[2]

  gridnx = if(length(gridn) == 1)round(gridn) else round(gridn[1])
  if(dim == 2){
    gridny = if(length(gridn) == 1 & gridnAd)round(gridn*diff(lims[3:4])/diff(lims[1:2])) else
      if(length(gridn) == 1)round(gridn) else round(gridn[2])
  }


  gx <- seq(lims[1], lims[2], length = gridnx) # gridpoints x

  #The outer product of the arrays X and Y is the array A with dimension c(dim(X), dim(Y))
  #where element A[c(arrayindex.x, arrayindex.y)] = FUN(X[arrayindex.x], Y[arrayindex.y], ...).
  # distance of each point to each grid point in x-direction
  ax <- matrix(rep(hx,gridnx), nrow=gridnx, ncol=nx, byrow=TRUE)^(-1)*outer(gx, x, "-")

  if(dim == 1){
    # z is the density
    z <- apply((matrix(rep(hx,gridnx), nrow=gridnx, ncol=nx, byrow=TRUE)^(-1)
                *matrix(rep(w/sum(w),gridnx), nrow=gridnx, ncol=nx, byrow=TRUE)
                *matrix(dnorm(ax), gridnx, nx)),1,sum)

    return(list(x = gx, z = z))

  }else if(dim == 2){
    gy <- seq(lims[3], lims[4], length = gridny) # gridpoints y
    # distance of each point to each grid point in y-direction
    ay <- matrix(rep(hy,gridny), nrow=gridny, ncol=nx, byrow=TRUE)^(-1)*outer(gy, y, "-")
    # z is the density
    z <- (matrix(rep(hx,gridnx), nrow=gridnx, ncol=nx, byrow=TRUE)^(-1)
          *matrix(rep(w,gridnx), nrow=gridnx, ncol=nx, byrow=TRUE)
          *matrix(dnorm(ax), gridnx, nx)) %*%
      t(matrix(rep(hy,gridny), nrow=gridny, ncol=nx, byrow=TRUE)^(-1)*matrix(dnorm(ay), gridny, nx))/(sum(w)) # z is the density

    return(list(x = gx, y = gy, z = z))

  }

}

