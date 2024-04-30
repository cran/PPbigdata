#Function to perform guided tour for data nuggets
guidedTourNugg <- function(nugg, weight, scale, dim, index = c("NH","Hole","CM"), qt = 0.8, ...){

  #check arguments
  index = match.arg(index)
  nugg = as.matrix(nugg)
  n = nrow(nugg)
  weight = as.vector(weight)
  scale = as.vector(scale)

  if(!dim %in% 1:2){
    stop("The dimension of projected values should be 1 or 2.")
  }

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

  #PP index function for guided tour
  if(index == "NH"){
    argnames <- names(list(...))
    NHargs <- intersect(argnames, names(as.list(args(NHnugg))))
    func = function(){function(mat){do.call(NHnugg,c(list(nuggproj = mat,weight,scale),list(...)[NHargs]))}}
  }else if(index == "Hole"){
    func = function(){function(mat){HoleNugg(nuggproj = mat,weight)}}
  }else if(index == "CM"){
    function(){function(mat){CMNugg(nuggproj = mat,weight)}}
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  #guided tour
  if(dim == 1){
    argnames <- names(list(...))
    gtargs <- intersect(argnames, names(as.list(args(tourr::guided_tour))))
    whargs <- intersect(argnames, names(c(as.list(args(weights::wtd.hist)),as.list(args(display_dist_nugg)))))
    anargs <- intersect(argnames, names(as.list(args(tourr::animate))))
    do.call(tourr::animate, c(list(nugg_wsph, do.call(tourr::guided_tour, c(list(func(),d = 1),list(...)[gtargs])),  display = do.call(display_dist_nugg,c(list(weight),list(...)[whargs])),sphere = FALSE),
                              list(...)[anargs]))
  }else if(dim == 2){
    argnames <- names(list(...))
    gtargs <- intersect(argnames, names(as.list(args(tourr::guided_tour))))
    anargs <- intersect(argnames, names(c(as.list(args(tourr::animate_xy)),as.list(args(tourr::display_xy)))))
    do.call(tourr::animate_xy, c(list(nugg_wsph, do.call(tourr::guided_tour, c(list(func(),d = 2),list(...)[gtargs])),sphere = FALSE,col = colors),
                              list(...)[anargs]))
  }
}


#subfunctions used

max_dist <- function(data, center = FALSE) {
  max(sqrt(rowSums(data^2)))
}

center <- function(x) {
  scale(x, center = TRUE, scale = FALSE)
}

compute_half_range <- function(half_range, data, center) {
  if (!is.null(half_range)) {
    return(half_range)
  }

  if (center) {
    data <- center(data)
  }
  half_range <- max_dist(data, center)
  half_range
}

display_dist_nugg <- function(weights, center = TRUE, half_range = NULL,
                              col = "black", density_max = 3.8, breaks = seq(-1, 1, 0.05), ...) {

  labels <- NULL
  init <- function(data) {
    half_range <- compute_half_range(half_range, data, center)
    labels <- abbreviate(colnames(data), 2)
    message("Using half_range ", format(half_range, digits = 2))
  }
  render_frame <- function() {
    par(pty = "m", mar = c(4, 4, 0, 1), mfrow = c(2, 1))
  }
  render_transition <- function() {
    rect(-1, -1.1, 1.2, density_max, col = "#FFFFFFE6", border = NA)
  }
  render_data <- function(data, proj, geodesic) {

    x <- data %*% proj
    if (center) x <- center(x)
    x <- x / compute_half_range(half_range, data, center)

    # Render projection data
    par(mar = c(0, 4, 1, 1))
    plot(
      x = NA, y = NA, xlim = c(-1, 1.2), ylim = c(0, density_max),
      xlab="", ylab = "Density",
      xaxt = "n", yaxt = "n"
    )
    axis(2, seq(0, density_max, by = 1))
    abline(h = seq(0, density_max, by = 0.5), col = "grey80")

    #weighted histogram
    bins <- weights::wtd.hist(x, weight = weights, plot = FALSE,breaks = breaks, ...)
    with(bins, rect(mids - unique(diff(bins$breaks))[1]/2, 0, mids + unique(diff(bins$breaks))[1]/2, density,
                    col = col
    ))

    par(mar = c(4, 4, 0, 1))
    plot(
      x = NA, y = NA, xlim = c(-1, 1.2), ylim = c(-1.1, 0),
      xaxs = "i", yaxs = "i",
      xlab = "", ylab = "Projection",
      yaxt = "n"
    )
    lines(c(0, 0), c(-1, 0), col = "grey80")
    lines(c(-1, -1), c(-1.1, 0), col = "grey80")
    lines(c(1, 1), c(-1.1, 0), col = "grey80")
    abline(h = 0)
    box(col = "grey70")

    # Render tour axes
    ax <- seq_along(proj) / length(proj)
    segments(0, -ax, proj, -ax, col = "black", lwd = 3)
    text(1.0, -ax, abbreviate(colnames(data), 2), pos = 4)
  }

  list(
    init = init,
    render_frame = render_frame,
    render_transition = render_transition,
    render_data = render_data,
    render_target = function(){}
  )
}


