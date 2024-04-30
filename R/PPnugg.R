## Function to perform projection persuit for big data

PPnugg <- function(data, index = c("NH","Hole","CM"), dim,
                   h = NULL, fa = TRUE, den = TRUE,
                   R = 5000, DN.num1 = 10^4, DN.num2 = 2000,max.splits = 5, seed_nugg = 5,
                   cooling = 0.9, tempMin = 1e-3, maxiter = 2000,
                   tol = 1e-6, seed_opt = 3, initP = NULL,
                   qt = 0.8,label = colnames(data), ...){

  #check arguments
  index = match.arg(index)
  label = as.vector(label)
  argnames <- names(list(...))

  if(!dim %in% 1:2){
    stop("The dimension of projected values should be 1 or 2.")
  }

  if(!is.null(label) & length(label) != ncol(data)){
    stop("The length of labels should be the same as the number of variables.")
  }

  if(qt <= 0 | qt > 1 | !is.numeric(qt)){
    stop("The quantile probability for weights transformation qt
           should be a number between 0 and 1.")
  }

  if(R > nrow(data) | DN.num1 > nrow(data) | DN.num2 > nrow(data)){
    stop("The inputs 'R', 'DN.num1' and 'DN.num2' for data nuggets should be
         less than the number of observations in data.")
  }

  if(!is.numeric(cooling) | cooling <= 0 | cooling >= 1){
    stop("The cooling rate 'cooling' should be a numeric value between (0,1). ")
  }

    if (!is.numeric(tempMin) || tempMin <= 0 || tempMin >= 1)
    stop("The minimal temperature, 'tempMin', should be a numeric value between (0,1).")

  if (!is.numeric(maxiter) || maxiter < 1 || (floor(maxiter)-maxiter) != 0)
    stop("The maximum number of iterations 'maxiter' should be a numeric integer greater than zero.")

  if(!is.numeric(tol) | tol <= 0){
    stop("The input 'tol' should be a positive numeric value.")
  }

  if(!is.null(initP)){
    initP = as.matrix(initP)
    if(nrow(initP) != ncol(data) | ncol(initP) != dim){
      stop("The initial projection matrix 'initP' should have the dimension of ncol(data)*dim.")
    }
  }


  #standardize raw data
  data = as.data.frame(scale(data))

  message("Create Data nuggets: \n")
  #create data nuggets
  my.DN = do.call(datanugget::create.DN, c(list(x = data,
                                        R = R,
                                        DN.num1 = DN.num1,
                                        DN.num2 = DN.num2, seed = seed_nugg),
                                   list(...)[intersect(argnames, names(as.list(args(datanugget::create.DN))))]))



  #refine data nuggets
  message("Refine Data nuggets: \n")
  my.DN2 = do.call(datanugget::refine.DN, c(list(x = data,
                                        DN = my.DN,
                                        max.splits = max.splits, seed = seed_nugg),
                                   list(...)[intersect(argnames, names(as.list(args(datanugget::refine.DN))))]))



  #get nugget centers, weights, and scales
  nugg = my.DN2$`Data Nuggets`[,2:(ncol(data)+1)]
  weight = my.DN2$`Data Nuggets`$Weight
  scale = my.DN2$`Data Nuggets`$Scale


  #spherize nugget centers considering weights
  wsph.res = wsph(nugg,weight)
  nugg_wsph = wsph.res$data_wsph
  nugg_wcen = wsph.res$data_wcen
  wsph_proj = wsph.res$wsph_proj

  #conduct the same projection on the standardized raw data
  data_cen = as.matrix(data)- as.matrix(rep(1,nrow(data)))%*%wsph.res$wmean


  message("\nOptimize PP index: \n")
  if(index == "NH"){
    res = do.call(PPnuggOptim, c(list(NHnugg, nugg_wsph, dimproj = dim,
                                cooling = cooling, tempMin = tempMin,
                                maxiter = maxiter, tol = tol,
                                seed = seed_opt , initP = initP, weight = weight, scale = scale, bandwidth = h),
                                list(...)[intersect(argnames, names(as.list(args(PPnuggOptim))))],
                                list(...)[intersect(argnames, names(as.list(args(NHnugg))))]))

  }else{
    func = if(index == "Hole")HoleNugg else CMNugg

    res = do.call(PPnuggOptim, c(list(func, nugg_wsph, dimproj = dim, weight = weight,
                                cooling = cooling, tempMin = tempMin,
                                maxiter = maxiter, tol = tol,
                                seed = seed_opt , initP = initP),
                           list(...)[intersect(argnames, names(as.list(args(PPnuggOptim))))]))

  }

  if(fa & dim == 2){
    fa_rotat = do.call(faProj, c(list(nugg,weight,proj = res$proj.opt),list(...)[intersect(argnames, names(as.list(args(faProj))))]))

    nuggproj = fa_rotat$nuggproj_rotat
    loadings = fa_rotat$loadings
    #nuggproj = nugg_wcen %*% loadings
  }else{
    nuggproj = res$proj.nugg
    loadings = wsph_proj%*%res$proj.opt
  }


  dataproj = data_cen%*%loadings

  do.call(plotNugg, c(list(nuggproj,weight,qt = qt),list(...)[intersect(argnames, names(as.list(args(plotNugg))))]))
  do.call(plotLoadings, c(list(loadings,label = label),list(...)[intersect(argnames, names(as.list(args(plotLoadings))))]))


  if(den)estw = do.call(nuggKDE, c(list(nuggproj,weight,scale),list(...)[intersect(argnames, names(as.list(args(nuggKDE))))]))

  output = list(DN = my.DN2, nuggproj = nuggproj, loadings = loadings,
                nugg_wcen = nugg_wcen, dataproj = dataproj, data_cen = data_cen,
                index.opt =  res$index[nrow(res$index),])

  if(den) output$density = estw

  return(output)

}
