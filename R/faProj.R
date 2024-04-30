##Function to make a factor rotation for nugget projection dim>=2----

faProj<-function(nugg, weight, wsph_proj = NULL, proj, method = c("varimax","promax")){

  nugg = as.matrix(nugg)
  weight = as.vector(weight)
  n = nrow(nugg)
  method = match.arg(method)

  if (length(weight) != n){
    stop("The length of weights should be the same as the number of nuggets. ")
  }

  if(sum(weight <=0) > 0){
    stop("The nugget weights should be positive values.")
  }


  wmean = 1/sum(weight)*t(as.matrix(weight))%*%nugg
  nugg_wcen = nugg - as.matrix(rep(1,n))%*%wmean

  if(is.null(wsph_proj)){
    wcov = 1/sum(weight)*t(as.matrix(nugg_wcen))%*%diag(weight)%*%as.matrix(nugg_wcen)
    ev = eigen(wcov)
    wsph_proj = ev$vectors%*%diag(ev$values^(-0.5))
  }else{
    wsph_proj = as.matrix(wsph_proj)
    if(nrow(wsph_proj) != ncol(nugg) | ncol(wsph_proj) != ncol(nugg)){
      stop("The spherization matrix should have the dimension of ncol(nugg)*ncol(nugg). ")
    }
  }

  proj = as.matrix(proj)
  proj_dim = ncol(proj)
  if(proj_dim < 2 | nrow(proj) != ncol(nugg)){
    stop("The projection matrix should have the dimension of
    ncol(nugg)*dim, where projection dimension dim should be larger than 1.")
  }

  proj_all = wsph_proj%*%proj

  if(method == "varimax"){
    rotat_proj = matrix(varimax(as.matrix(proj_all))$loadings,ncol = proj_dim)
  }else{
    rotat_proj = matrix(promax(as.matrix(proj_all))$loadings,ncol = proj_dim)
  }


  nuggproj_rotat = as.matrix(nugg_wcen)%*%as.matrix(rotat_proj)
  loadings = rotat_proj

  return(list(nuggproj_rotat = nuggproj_rotat, loadings = loadings, nugg_wcen = nugg_wcen))
}
