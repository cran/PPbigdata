##Function to spherize data with weights
wsph <- function(data,weight = rep(1, nrow(data))){

  data = as.matrix(data)
  weight = as.vector(weight)

   # make sure weights is of class numeric or table or integer
  if (!is.numeric(weight)&
      !is.table(weight) &
      !is.integer(weight)){

    stop('weight must be of class "numeric" or "table"')

  }

  # make sure weight is length nrow(data)
  if (length(weight) != nrow(data)){
    stop(stop("The length of weights should be the same as the number of observations. "))
  }

  n = nrow(data)
  wmean = 1/sum(weight)*t(as.matrix(weight))%*%data
  data_wcen = data - as.matrix(rep(1,n))%*%wmean
  wcov = 1/sum(weight)*t(as.matrix(data_wcen))%*%diag(weight)%*%as.matrix(data_wcen)
  ev = eigen(wcov)
  wsph_proj = ev$vectors%*%diag(ev$values^(-0.5))
  data_wsph = as.matrix(data_wcen)%*%wsph_proj
  if(!is.null(colnames(data)))colnames(data_wsph) = colnames(data)

  return(list(data_wsph = data_wsph, data_wcen = data_wcen, wmean = wmean, wcov = wcov, wsph_proj = wsph_proj))
}
