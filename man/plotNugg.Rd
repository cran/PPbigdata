\name{plotNugg}
\alias{plotNugg}
\title{
Plot of projected 1-dim/2-dim data nuggets
}
\description{
Draw a scatterplot/stripchart/weighted histogram of projected 1-dim/2-dim data nuggets considering the weights of data nuggets.
}
\usage{
plotNugg(nuggproj, weight,
         qt = 0.8,pch = 16,cex = 0.5,
         hist = FALSE, jitter = 0.1,freq = TRUE, breaks = 30,
         main = NULL, xlab = NULL, ylab = NULL,...)
}
\arguments{

  \item{nuggproj}{Projected data nugget centers in 1-dim/2-dim space. Must be a data matrix (of class matrix, or data.frame) with two columns or a vector containing only entries of class numeric.
  }
  \item{weight}{Vector of the weight parameter for each data nugget. Its length should be the same as the number of data nuggets, i.e., nrow(nuggproj) for 2-dim projection/length(nuggproj) for 1-dim projection. Must be of class numeric or integer.
  }
  \item{qt}{A scalar with value in \code{[0,1]} indicating the probability used to obtain a sample quantile of the data nugget weights as the maximal value to transform weights to colors for the plot. Defaults to be 0.8.
  }
  \item{pch}{Either an integer specifying a symbol or a single character to be used as the default in plotting points. See \code{\link[graphics]{points}} for possible values and their interpretation.
  }
  \item{cex}{A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default. See \code{\link[graphics]{par}}.
  }
  \item{hist}{logical; If TRUE, a weighted histogram is plotted for 1-dim projected data nuggets considering the nugget weights. Otherwise, a stripchart of the points jittered is plotted for 1-dim projection with colors indicating nugget weights. Defalut to be FALSE. Ignorable for 2-dim projection.
  }
  \item{jitter}{If \code{hist = FALSE}(default), the amount of jittering applied to the stripchart of 1-dim projected data nuggets. Ignorable for 2-dim projection and 1-dim projection with \code{hist = TRUE}.
  }
  \item{freq}{If \code{hist = TRUE}, a logical value indicating whether to plot frequencies. or probability densities for the weighted histogram of 1-dim projected data nuggets. Ignorable for 2-dim projection and 1-dim projection with \code{hist = FALSE}.
  }
  \item{breaks}{If \code{hist = TRUE}, the \code{breaks} argument for the weighted histogram of 1-dim projected data nuggets. See details in \code{\link[weights]{wtd.hist}}. Ignorable for 2-dim projection and 1-dim projection with \code{hist = FALSE}.
  }
  \item{main}{an overall title for the plot: see \code{\link[graphics]{title}}.
  }
   \item{xlab}{a title for the x axis: see \code{\link[graphics]{title}}.
  }
   \item{ylab}{a title for the y axis: see \code{\link[graphics]{title}}.
  }
  \item{...}{further arguments and graphical parameters passed to \code{\link{plot}} for 2-dim projection, or \code{\link[graphics]{stripchart}} or \code{\link[weights]{wtd.hist}} for 1-dim projection.
  }
}
\details{
This function plots a scatterplot/stripchart/weighted histogram of projected 1-dim/2-dim data nuggets considering the weights of data nuggets.

Data nuggets are a representative sample meant to summarize Big Data by reducing a large dataset to a much smaller dataset by eliminating redundant points while also preserving the peripheries of the dataset. Each data nugget is defined by a center (location), weight (importance), and scale (internal variability). Data nuggets for a large dataset could be created and refined by functions \code{create.DN} or \code{refine.DN} in the package \code{datanugget}.

Based on the created and refined data nuggets, the projected 1-dim/2-dim data nuggets are plotted via a scattorplot (2d) or a stripchart with points jittered (1d) of data nugget centers, coloring with data nugget weights where lighter green represents larger weights. If \code{hist = TRUE}, a weighted histgram of 1-dim projected data nugget centers is plotted considering the data nugget weights.
}
\value{
If \code{hist = TRUE}, an object of class \code{histogram}. See \code{\link[weights]{wtd.hist}}. Otherwise, no returned values.
}
\references{
Cherasia, K. E., Cabrera, J., Fernholz, L. T., & Fernholz, R. (2022). Data Nuggets in Supervised Learning. \emph{In Robust and Multivariate Statistical Methods: Festschrift in Honor of David E. Tyler} (pp. 429-449). Cham: Springer International Publishing.

Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.
}
\author{
Yajie Duan, Javier Cabrera
}
\seealso{
\code{\link[datanugget]{datanugget-package}}, \code{\link[datanugget]{create.DN}}, \code{\link[datanugget]{refine.DN}},\code{\link[weights]{wtd.hist}}, \code{\link[graphics]{stripchart}}
}

\examples{
  require(datanugget)
  require(rstiefel)

  #2-d small example with visualization
  X = rbind.data.frame(matrix(rnorm(10^4, sd = 0.3), ncol = 2),
            matrix(rnorm(10^4, mean = 1, sd = 0.3), ncol = 2))


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

  #plot raw large dataset
  plot(X)

  #plot data nuggets in 2-dim space
  plotNugg(nugg,weight)

  #generate a random projection vector to 1-dim space
  proj_1d = rustiefel(2, 1)

  #project data nugget centers into 1-dim space by the random projection vector
  nuggproj_1d = as.matrix(nugg)\%*\%proj_1d

  #plot the stripchart for 1-dim projected data nuggets
  plotNugg(nuggproj_1d,weight)
  #plot the weighted histogram for 1-dim projected data nuggets
  plotNugg(nuggproj_1d,weight,hist = TRUE)

}


