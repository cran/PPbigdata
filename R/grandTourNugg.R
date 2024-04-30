
grandTourNugg <- function(nugg, weight, dim, qt = 0.8, ...){

  #check arguments
  nugg = as.matrix(nugg)
  n = nrow(nugg)
  weight = as.vector(weight)

  if(!dim %in% 1:2){
    stop("The dimension of projected values should be 1 or 2.")
  }

  if (length(weight) != n){
    stop("The length of weights should be the same as the number of nuggets. ")
  }

  if(sum(weight <=0) > 0){
    stop("The nugget weights should be positive values.")
  }

  if(qt <= 0 | qt > 1 | !is.numeric(qt)){
    stop("The quantile probability for weights transformation qt
           should be a number between 0 and 1.")
  }

  #spherize data nugget centers considering weights
  wsph.res = wsph(nugg,weight)
  nugg_wsph = wsph.res$data_wsph
  colnames(nugg_wsph) = colnames(nugg)


  #get data nugget colors based on weights for 2-d projection plots
  weight_trans = weight/sum(weight)
  weight_trans = weight_trans/quantile(weight_trans,qt)
  colors =sapply(1:n, function(t){rgb(0,min(weight_trans[t],1),0)})

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  #grand tour
  if(dim == 1){
    argnames <- names(list(...))
    whargs <- intersect(argnames, names(c(as.list(args(weights::wtd.hist)),as.list(args(display_dist_nugg)))))
    anargs <- intersect(argnames, names(as.list(args(tourr::animate))))
    do.call(tourr::animate, c(list(nugg_wsph, tourr::grand_tour(d = 1),  display = do.call(display_dist_nugg,c(list(weight),list(...)[whargs])),sphere = FALSE),
                              list(...)[anargs]))
  }else if(dim == 2){
    argnames <- names(list(...))
    anargs <- intersect(argnames, names(c(as.list(args(tourr::animate_xy)),as.list(args(tourr::display_xy)))))
    do.call(tourr::animate_xy, c(list(nugg_wsph, tourr::grand_tour(d = 2),sphere = FALSE,col = colors),
                                 list(...)[anargs]))
  }
}
