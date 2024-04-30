\name{nuggKDE}
\alias{nuggKDE}
\title{
Density Estimation for projected 1-dim/2-dim big data based on data nuggets
}
\description{
This function estimates the density function of projected 1-dim/2-dim big data based on data nuggets and kernal density esimation.
}
\usage{
nuggKDE(nuggproj, weight, scale, h = NULL, gridn = 300,lims = NULL, gridnAd = TRUE)
}
\arguments{
  \item{nuggproj}{Projected data nugget centers in 1-dim/2-dim space. Must be a data matrix (of class matrix, or data.frame) with two columns or a vector containing only entries of class numeric.
  }
  \item{weight}{Vector of the weight parameter for each data nugget. Its length should be the same as the number of data nuggets, i.e., nrow(nuggproj) for 2-dim projection/length(nuggproj) for 1-dim projection. Must be of class numeric or integer.
  }
  \item{scale}{Vector of the scale parameter for each data nugget. Its length should be the same as the number of data nuggets, i.e., nrow(nuggproj) for 2-dim projection/length(nuggproj) for 1-dim projection. Must be of class numeric or integer.
  }
  \item{h}{Bandwidth in each direction that would be combined with data nuggets scales as the final bandwith for kernal density estimation of projected data nuggets. Defaults to normal reference bandwidth considering the nugget weights, i.e., \code{\link[MASS]{bandwidth.nrd}(rep(dat,weight))}. Can be scalar or a length-2 numeric vector. For 2-dim projection, a scalar value will be applied on both directions.
  }
  \item{gridn}{Number of grid points in each direction used for kernel density estimation of projected data. Can be scalar or a length-2 integer vector.
  }
  \item{lims}{The limits of each direction used for kernel density estimation of projected data. Must be a length-4 numeric vector as (xlow, xupper, ylow, yupper) for 2-dim projected data, or a length-2 numeric vector as (xlow, xupper) for 1-dim projected data. If NULL, defaulting to the range of each direction.
  }
  \item{gridnAd}{logical; if \code{TRUE} (default) and \code{gridn} is a scalar, for 2-dim projected data \code{rawproj}, \code{gridn} is used for x-direction, and the number of grid points in y-direction is adjusted by the limits of both directions, i.e., \code{round(gridn*diff(lims[3:4])/diff(lims[1:2]))}. Ignorable when \code{gridn} is a length-2 integer vector or projected data \code{rawproj} is 1-dim.
  }
}
\details{
  This function calculates the estimated density values of projected 1-dim/2-dim big based on data nuggets.

  Data nuggets are a representative sample meant to summarize Big Data by reducing a large dataset to a much smaller dataset by eliminating redundant points while also preserving the peripheries of the dataset. Each data nugget is defined by a center (location), weight (importance), and scale (internal variability). Data nuggets for a large dataset could be created and refined by functions \code{create.DN} or \code{refine.DN} in the package \code{datanugget}.

  Based on the created and refined data nuggets, the density of projected 1-dim/2-dim big data could be estimated via a revised version of kernal density estimation considering the data nugget centers, weightes and scales. For the estimation, the normal kernal is used with a bandwidth being a combination of pre-specified bandwidth and scales of data nuggets. By default, the pre-specified bandwidth in each direction is a normal reference bandwidth considering the nugget weights, i.e., \code{\link[MASS]{bandwidth.nrd}(rep(dat,weight))}.
}
\value{
A list containing the following components:

\item{x}{The coordinates of the grid points on x-direction.}

\item{y}{For 2-dim projection, the coordinates of the grid points on y-direction. Non-existing for 1-dim projection.}

\item{z}{For 2-dim projection, a \code{length(x)} by \code{length(y)} martix of the estimated density. For 1-dim projection, a vector of length \code{length(x)} of the estimated density values.}

}
\references{
  Duan, Y., Cabrera, J. & Emir, B. A New Projection Pursuit Index for Big Data. Under revision.

  Cherasia, K. E., Cabrera, J., Fernholz, L. T., & Fernholz, R. (2022). Data Nuggets in Supervised Learning. In Robust and Multivariate Statistical Methods: Festschrift in Honor of David E. Tyler (pp. 429-449). Cham: Springer International Publishing.

  Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.
}
\author{
  Yajie Duan, Javier Cabrera
}
\seealso{
\code{\link[PPbigdata]{NHnugg}},\code{\link[datanugget]{create.DN}}, \code{\link[datanugget]{refine.DN}}
}
\examples{

  require(datanugget)
  require(rstiefel)

  #4-dim small example
  X = cbind.data.frame(rnorm(5*10^3),
                       rnorm(5*10^3),
                       rnorm(5*10^3),
                       rnorm(5*10^3))


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

  #generate a random orthonormal matrix as a projection matrix to 2-dim space
  proj_2d = rustiefel(4, 2)

  #project data nugget centers into 2-dim space by the random projection matrix
  nuggproj_2d = as.matrix(nugg)\%*\%proj_2d

  #plot the projected data nuggets
  #lighter green represents larger weights
  plotNugg(nuggproj_2d, weight)

  #project raw large data into 2-dim space using the same projection matrix
  rawproj_2d = as.matrix(X)\%*\%proj_2d

  #plot projected raw large dataset
  plot(rawproj_2d)

  #estimated density for 2-dim projected data based on the data nuggets
  est_nugg = nuggKDE(nuggproj_2d, weight, scale)
  #plot the estimated density values
  image(est_nugg)
}

