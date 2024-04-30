\name{plotLoadings}
\alias{plotLoadings}
\title{
Plot of variable loadings for 1-dim/2-dim projection
}
\description{
Draw a loading plot for 1-dim/2-dim projection
}
\usage{
plotLoadings(loadings, label = NULL,pch = 16,main = "Loadings",
             xlab = NULL,ylab = NULL, xlim = NULL, ylim = NULL,
             textpoints = loadings, textcex = 1, textpos = 2, textcol = NULL, ...)
}
\arguments{
  \item{loadings}{Loadings of original variables for projection in 1-dim/2-dim space, one column for each projection direction. Must be a data matrix (of class matrix, or data.frame) with one or two columns or a vector containing only entries of class numeric.}
  \item{label}{A character vector specifying the text to be written for the variables. Default is NULL, which uses the "V1, V2,..." as the labels of variables. If not NULL, the vector length should be the number of variables, i.e., nrow(loadings).
  }
  \item{pch}{For 2-dim projection loadings, either an integer specifying a symbol or a single character to be used as the default in plotting points. See \code{\link[graphics]{points}} for possible values and their interpretation. Default is 16.
  }
  \item{main}{The title for the loading plot.
  }
  \item{xlab}{The title for the x axis.
  }
  \item{ylab}{The title for the y axis.
  }
  \item{xlim}{The scale limits for the x axis.
  }
  \item{ylim}{The scale limits for the y axis.
  }
  \item{textpoints}{For 2-dim projection loadings, the x and y coordinates of the variable label text displayed on the loading plot. Must be a data matrix (of class matrix, or data.frame) with two columns. Default value is the 2-dim projection loading matrix.
  }
  \item{textcex}{For 2-dim projection loadings, numeric character expansion factor for the variable label text displayed on the loading plot. Default is 1. See \code{\link[graphics]{text}}.
  }
  \item{textpos}{For 2-dim projection loadings, a position specifier for the variable label text displayed on the loading plot. Default is 2, indicating the left of the coordinates. See \code{\link[graphics]{text}}.
  }
  \item{textcol}{For 2-dim projection loadings, the color font to be used for the variable label text displayed on the loading plot. See \code{\link[graphics]{text}}.
  }
  \item{...}{other arguments and graphical parameters passed to \code{\link{plot}} for 2-dim projection, or \code{\link[graphics]{barplot}} for 1-dim projection.
  }
}
\details{
This function draws a loading plot for 1-dim/2-dim projection. It plot a barplot for 1-dim projection loadings and a scatterplot with variable label text for 2-dim.
}
\value{
No return value, called for plotting.
}
\references{
Duan, Y., Cabrera, J., & Emir, B. (2023). A New Projection Pursuit Index for Big Data. ArXiv:2312.06465. https://doi.org/10.48550/arXiv.2312.06465

Cherasia, K. E., Cabrera, J., Fernholz, L. T., & Fernholz, R. (2022). Data Nuggets in Supervised Learning. \emph{In Robust and Multivariate Statistical Methods: Festschrift in Honor of David E. Tyler} (pp. 429-449). Cham: Springer International Publishing.

Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.
}
\author{
Yajie Duan, Javier Cabrera
}
\seealso{
\code{\link[PPbigdata]{PPnugg}}, \code{\link[PPbigdata]{faProj}}
}
\examples{
  require(datanugget)
  require(rstiefel)

  #4-d small example with visualization
  X = rbind.data.frame(matrix(rnorm(5*10^3, sd = 0.3), ncol = 4),
            matrix(rnorm(5*10^3, mean = 1, sd = 0.3), ncol = 4))


  #create data nuggets
  my.DN = create.DN(x = X,
                    R = 500,
                    delete.percent = .1,
                    DN.num1 = 500,
                    DN.num2 = 250,
                    no.cores = 0,
                    make.pbs = FALSE)


  #refine data nuggets
  my.DN2 = refine.DN(x = X,
                     DN = my.DN,
                     EV.tol = .9,
                     min.nugget.size = 2,
                     max.splits = 5,
                     no.cores = 0,
                     make.pbs = FALSE)

  #get nugget centers, weights, and scales
  nugg = my.DN2$`Data Nuggets`[,2:(ncol(X)+1)]
  weight = my.DN2$`Data Nuggets`$Weight
  scale = my.DN2$`Data Nuggets`$Scale


  #generate a random projection matrix to 2-dim space
  proj_2d = rustiefel(4, 2)

  #project data nugget centers into 2-dim space by the random projection matrix
  nuggproj_2d = as.matrix(nugg)\%*\%proj_2d

  #plot data nuggets in 2-dim space
  plotNugg(nuggproj_2d,weight)

  #plot loadings for the variables
  plotLoadings(proj_2d)


  #generate a random projection vector to 1-dim space
  proj_1d = rustiefel(4, 1)

  #project data nugget centers into 1-dim space by the random projection vector
  nuggproj_1d = as.matrix(nugg)\%*\%proj_1d

  #plot the weighted histogram for 1-dim projected data nuggets
  plotNugg(nuggproj_1d,weight,hist = TRUE)

  #plot loadings for the variables
  plotLoadings(proj_1d)

}


