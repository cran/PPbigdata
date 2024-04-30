##Function to plot loadings of 1d or 2d projection

plotLoadings <- function(loadings, label = NULL,pch = 16,main = "Loadings",
                         xlab = NULL,ylab = NULL, xlim = NULL, ylim = NULL,
                         textpoints = loadings, textcex = 1, textpos = 2, textcol = NULL, ...){

  #check arguments
  loadings = as.matrix(loadings)
  textpoints = as.matrix(textpoints)
  dim <- ncol(loadings)

  if(!dim %in% 1:2){
    stop("The dimension of projection should be 1 or 2.")
  }


  if(is.null(label)){label = paste("V",1:nrow(loadings),sep = "")}else if(length(label) != nrow(loadings)){
    stop("The length of labels should be the same as the number of variables.")
  }

  if(is.null(xlab)){xlab = if(ncol(loadings) == 2)"Projection 1" else ""}
  if(is.null(ylab)){ylab = if(ncol(loadings) == 2)"Projection 2" else "Loading"}
  if(is.null(ylim)){ylim = if(ncol(loadings) == 2)c(min(loadings[,2])-0.05,
                                                    max(loadings[,2])+0.05) else c(min(loadings)-0.1,
                                                                                   max(loadings)+0.1)}

  if(dim == 1){

    barplot(loadings ~ label, xlab = xlab, ylab = ylab, ylim = ylim,main = main,...)
    abline(h = 0,col = "red")

  }else if(dim == 2){

    if(is.null(xlim)){xlim = c(min(loadings[,1])-0.2,max(loadings[,1])+0.2)}

    plot(loadings, pch = pch,xlim = xlim,ylim = ylim
         ,xlab =xlab,ylab = ylab, main = main,...)
    abline(h = 0,v = 0,lty = 4)
    text(textpoints,labels = label, cex = textcex, pos = textpos, col = textcol)

  }
}

