##Optimize PP index for 1d/2d projection for data nuggets----
#GTSA optimization: grand tour simulated annealing optimization method

PPnuggOptim<- function(FUN, nugg_wsph, dimproj, tempInit = 1, cooling = 0.9, eps = 1e-3, tempMin = 0.01,
                            maxiter = 1000, half = 10, tol = 1e-5, maxc = 15,
                            seed = 3, initP = NULL, ...) {

  FUN <- match.fun(FUN)

  if (!is.data.frame(nugg_wsph) && !is.matrix(nugg_wsph))
    stop("The input 'nugg_wsph' should be of type data frame or matrix.")


  data = as.matrix(nugg_wsph)

  if(!is.numeric(dimproj) | dimproj < 1){
    stop("The dimension of projection should be a positive integar.")
  }

  if(!is.numeric(tempInit) | tempInit <= 0 | tempInit >1){
     stop("The initial temperature 'tempInit' should be a numeric value between (0,1). ")
  }

  if(!is.numeric(cooling) | cooling <= 0 | cooling >= 1){
    stop("The cooling rate 'cooling' should be a numeric value between (0,1). ")
  }

  if (!is.numeric(eps) || eps <= 0 || eps >= 1)
    stop("The approximation accuracy for cooling, 'eps', should be a numeric value between (0,1).")

  if (!is.numeric(tempMin) || tempMin <= 0 || tempMin >= 1)
    stop("The minimal temperature, 'tempMin', should be a numeric value between (0,1).")

  if (!is.numeric(maxiter) || maxiter < 1 || (floor(maxiter)-maxiter) != 0)
    stop("The maximum number of iterations 'maxiter' should be a numeric integer greater than zero.")

  if (!is.numeric(half) || half < 1 || (floor(half)-half) != 0)
    stop("The input 'half' should be a numeric integer greater than zero.")

  if(!is.numeric(tol) | tol <= 0){
    stop("The input 'tol' should be a positive numeric value.")
  }

  if (!is.numeric(maxc) || maxc < 1 || (floor(maxc)-maxc) != 0)
    stop("The input 'maxc' should be a numeric integer greater than zero.")



  #set seeds
  set.seed(seed)

  Dat <- as.matrix(data)
  NumCol <- ncol(Dat)

  if(!is.null(initP)){
    if(nrow(initP) != NumCol | ncol(initP) != dimproj){
      stop("The initial projection matrix 'initP' should have the dimension of ncol(nugg_wsph)*dimproj.")
    }
    Aa <- initP
  }else{
    Aa <- diag(1, nrow = NumCol, ncol = dimproj) # Matrix of orthogonal initialization
  }

  Proj <- Dat %*% Aa # initial projection
  Proj <- as.matrix(Proj)

  proj.data <- Proj

  indexMax <- FUN(Proj, ...)#calculate index value for projection
  index <- as.matrix(indexMax)

  mi <- 1
  h  <- 0	# number of iterations without increase in index
  e <- 0 # number of different temp with no increase or increase smaller than tol in index
  temp <- tempInit

  while (mi <= maxiter && temp > tempMin && e < maxc) {

    Ai <- Base(NumCol, dimproj) # initial base

    Az <- Aa + temp * Ai # target projection

    Az <-  Interpolation(Aa, Az, eps) # Projection matrix through interpolation

    Proj <- Dat %*% Az

    indexC <- FUN(Proj, ...)

    message(paste ("Iteration <-", round(mi,1),"   index <-", round(indexMax,10), "   temp <-", round(temp,9)))

    if (indexC > indexMax) {
      diff      <- indexC - indexMax
      if(diff > tol) e <- 0
      Aa        <- Az
      indexMax  <- indexC
      proj.data <- Proj
      index     <- rbind(index, indexC)
      h <- 0
    } else{
      h <- h + 1
    }

    mi <- mi + 1

    if (h == half) {
      temp <- temp * cooling
      e <- e + 1
      h <- 0
    }



  }


  rownames(index) <- NULL

  rownames(proj.data)  <- rownames(data)

  rownames(Aa) <- colnames(data)

  colnames(proj.data) <- c(paste("Projection", 1:(ncol(proj.data))))

  colnames(Aa) <- c(paste("Axis", 1:(ncol(Aa))))


  Lista <- list(proj.nugg = proj.data, proj.opt = Aa, index = index)

  return(Lista)


}



Base <- function(NumLin, d) {
  #This function helps to find an orthonormal base using the principal components
  dataBase <- matrix(rnorm(NumLin * d), ncol = d)

  return(dataBase)
}


Interpolation <- function(Aa, Az, eps) {
  # This function performs matrix Aa interpolation in Az

  # Input:
  # Aa - Initial projection
  # Az - Projection target

  # Return:
  # A - Matrix of interpolated projection of Aa in Az

  if (!is_orthonormal(Aa)) Aa <- Orthonormalise(Aa)

  if (!is_orthonormal(Az)) Az <- Orthonormalise(Az)

  sv <- svd(t(Aa) %*% Az)

  NumCol <- ncol(Aa)
  lambda <- sv$d[NumCol:1]
  Va     <- sv$u[, NumCol:1]
  Vz     <- sv$v[, NumCol:1]

  Ba <- Aa %*% Va
  Bz <- Az %*% Vz

  Ba <- Orthonormalise(Ba)
  Bz <- Orthonormalise(Bz)
  Bz <- Orthonormalise_by(Bz, Ba)

  Tau <- suppressWarnings(acos(lambda))
  Tau.NaN <- is.nan(Tau) | Tau < eps
  Tau[Tau.NaN] <- 0

  Bz[, Tau.NaN] <- Ba[, Tau.NaN]

  k <- 1
  # for (k in 1:length(Tau))
  while (k <= length(Tau)) {
    Bz[,k] <- Ba[,k] * cos(Tau[k]) + Bz[,k] * sin(Tau[k])
    k <- k + 1
  }

  A = Bz %*% Va

  return(A)
}


Normalise <- function(Base) {
  # This function normalizes a Base
  sweep(Base, 2, sqrt(colSums(Base^2,na.rm = TRUE)), FUN = "/")
}


Orthonormalise_by <- function(MatX, MatY) {
  # This function Orthonormalize one matrix for another.
  # ensures that each MatX column is orthogonal to the column
  # correspondent in by MatY, using the Gram-Shimidt process

  stopifnot(ncol(MatX) == ncol(MatY))
  stopifnot(nrow(MatX) == nrow(MatY))

  MatX <- Normalise(MatX)
  j <- 1
  while(j <= ncol(MatX)) {
    MatX[, j] <- MatX[, j] - c(crossprod(MatX[, j], MatY[, j]) / sum(MatY[, j]^2)) * MatY[, j]
    j <- j + 1
  }

  Normalise(MatX)
}


Orthonormalise <- function(Base) {
  # This function finds an orthogonal or orthonormal basis
  # for Base vectors using the Gram-Shimidt process

  Base <- Normalise(Base) # to be conservative

  if (ncol(Base) > 1) {
    j <- 1
    while(j <= ncol(Base)) {
      i <- 1
      while(i <= (j - 1)) {
        Base[, j] <- Base[, j] - c(crossprod(Base[, j], Base[, i]) / sum(Base[, i]^2)) * Base[, i]
        i <- i + 1
      }
      j <- j + 1
    }
  }

  Normalise(Base)
}


is_orthonormal <- function(data) {

  stopifnot(is.matrix(data))

  Limit <- 0.001

  j <- 1
  while(j <= ncol(data)) {
    if (sqrt(sum(data[, j] ^ 2)) < 1 - Limit) return(FALSE)
    j <- j + 1
  }

  if (ncol(data) > 1) {
    j <- 2
    while(j <= ncol(data)) {
      i <- 1
      while(i <= (ncol(data) - 1)) {
        if (abs(sum(data[, j] * data[, i])) > Limit) return(FALSE)
        i <- i + 1
      }
      j <- j + 1
    }
  }

  TRUE
}


