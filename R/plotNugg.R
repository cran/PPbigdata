##Function to plot projected data nuggets
library(dplyr)
library(weights)

plotNugg <- function(nuggproj,weight,qt = 0.8,pch = 16,cex = 0.5, hist = FALSE,
                     jitter = 0.1,freq = TRUE, breaks = 30, main = NULL, xlab = NULL, ylab = NULL,...){

    #check arguments

    nuggproj = as.matrix(nuggproj)
    dim <- ncol(nuggproj)
    weight = as.vector(weight)

    if(!dim %in% 1:2){
      stop("The dimension of projected values should be 1 or 2.")
    }


    if (length(weight) != nrow(nuggproj)){
      stop("The length of weights should be the same as the number of nuggets. ")
    }

    if(sum(weight <=0) > 0){
      stop("The nugget weights should be positive values.")
    }


    if(qt <= 0 | qt > 1 | !is.numeric(qt)){
      stop("The quantile probability for weights transformation qt
           should be a number between 0 and 1.")
    }

    #set up colors by nugget weights
    weight_trans = weight/sum(weight)
    weight_trans = weight_trans/quantile(weight_trans,qt)
    col=sapply(1:nrow(nuggproj), function(t){rgb(0,min(weight_trans[t],1),0)})


    #set up defalt title and axis labels
    if(is.null(main)){main = if(dim == 1 & hist)"Weighted Histogram of Projected Values \n considering data nugget weights"
                             else "Projected Data Nuggets"}
    if(is.null(xlab)){xlab = if(dim == 2)"Projection 1" else if(hist)"Projected Value" else "Data nugget"}
    if(is.null(ylab)){ylab = if(dim == 2)"Projection 2" else if(hist)"" else "Projected Value"}


    if(dim == 2){
      #plot for 2-d projection
       plot(nuggproj,col=col,pch=pch, cex=cex,xlab = xlab,
           ylab = ylab,main =main,...)
    }else if(dim == 1){
      if(!hist){
        #stripchart for 1-d projection
        V1=as.numeric(factor(col))
        V2=as.numeric(nuggproj)
        df = data.frame(V1=V1, V2=V2,col = col)
        df = df[order(nuggproj),]
        group_by(df,V1) %>% group_walk(~{stripchart(pull(.,V2), method="jitter",
                                                    col=I(pull(.,col)),
                                                    vertical=FALSE, pch=pch,
                                                    add=min(pull(.,V1))!=1,
                                                    cex=cex, jitter=jitter,
                                                    at=1,
                                                    xlab = xlab,ylab = ylab,
                                                    main = main,...)},
                                       .keep=TRUE)
      }else{
        #weighted histogram for 1-d projection
        wtd.hist(nuggproj, weight = weight, freq = freq, breaks = breaks,
                 xlab = xlab, ylab = ylab, main = main, ...)

      }

    }
}
