\name{NHnugg}
\alias{NHnugg}
\title{
Natural Hermite Index for projected 1-dim/2-dim Big Data based on data nuggets
}
\description{
This function calculates the value of Nature Hermite index, a Projection Pursuit index proposed by Cook(1993) for projected 1-dim/2-dim big data based on data nuggets.
}
\usage{
NHnugg(nuggproj, weight, scale,
       bandwidth = NULL, gridn = 300,lims = NULL, gridnAd = TRUE)}
\arguments{
\item{nuggproj}{Projected data nugget centers in 1-dim/2-dim space. Must be a data matrix (of class matrix, or data.frame) with two columns or a vector containing only entries of class numeric.
  }
  \item{weight}{Vector of the weight parameter for each data nugget. Its length should be the same as the number of data nuggets, i.e., nrow(nuggproj) for 2-dim projection/length(nuggproj) for 1-dim projection. Must be of class numeric or integer.
  }
  \item{scale}{Vector of the scale parameter for each data nugget. Its length should be the same as the number of data nuggets, i.e., nrow(nuggproj) for 2-dim projection/length(nuggproj) for 1-dim projection. Must be of class numeric or integer.
  }
  \item{bandwidth}{Bandwidth in each direction that would be combined with data nuggets scales as the final bandwith for kernal density estimation of projected data nuggets. Defaults to normal reference bandwidth considering the weights. Can be scalar or a length-2 numeric vector. For 2-dim projection, a scalar value will be applied on both directions.
  }
  \item{gridn}{Number of grid points in each direction used for kernel density estimation of projected data. Can be scalar or a length-2 integer vector.
  }
  \item{lims}{The limits of each direction used for kernel density estimation of projected data. Must be a length-4 numeric vector as (xlow, xupper, ylow, yupper) for 2-dim projected data, or a length-2 numeric vector as (xlow, xupper) for 1-dim projected data. If NULL, defaulting to the range of each direction.
  }
  \item{gridnAd}{logical; if \code{TRUE} (default) and \code{gridn} is a scalar, for 2-dim projected data \code{rawproj}, \code{gridn} is used for x-direction, and the number of grid points in y-direction is adjusted by the limits of both directions, i.e., \code{round(gridn*diff(lims[3:4])/diff(lims[1:2]))}. Ignorable when \code{gridn} is a length-2 integer vector or projected data \code{rawproj} is 1-dim.
  }
}
\details{
  This function calculates the value of Nature Hermite index, a Projection Pursuit index proposed by Cook(1993) for projected 1-dim/2-dim Big Data based on data nuggets.

  Data nuggets are a representative sample meant to summarize Big Data by reducing a large dataset to a much smaller dataset by eliminating redundant points while also preserving the peripheries of the dataset. Each data nugget is defined by a center (location), weight (importance), and scale (internal variability). Data nuggets for a large dataset could be created and refined by functions \code{create.DN} or \code{refine.DN} in the package \code{datanugget}.

  Natural Hermite index is one kind of Projection Pursuit (PP) index, and it measures the distance between the density of projected data and the standard normal. Larger index values indicate a hidden structure of multivariate data, such as clustersing, outliers or other non-linear structures. However, it's computationally hard to calculate the index for big data because of the issue about density estimation of projected big data. A new PP index for big data was proposed by Duan(2023), which is based on the Natural Hermite index and data nuggets.

  In this function, the PP index value for projected 1-dim/2-dim big data is calculated based on created and refined data nuggets. Data nuggets are firstly created and refined for the big data. For Natural Hermite index, the data nugget centers need to be spherized considering nugget weights before projection. The projection is taken on the spherized data nugget centers to obtain projected data nuggets. The density values of projected big data are firstly estimated by \code{\link[PPbigdata]{nuggKDE}}. Based on it, the Natural Hermite index value is calculated via numerical integral by summation.
}
\value{
  A numeric value indicating Nature Hermite index value of the projected big data based on the data nuggets.
}
\references{
  Cook, D., Buja, A., & Cabrera, J. (1993). Projection pursuit indexes based on orthonormal function expansions. Journal of Computational and Graphical Statistics, 2(3), 225-250.

  Cook, D., Buja, A., Cabrera, J., & Hurley, C. (1995). Grand tour and projection pursuit. Journal of Computational and Graphical Statistics, 4(3), 155-172.

  Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.

  Duan, Y., Cabrera, J., & Emir, B. (2023). A New Projection Pursuit Index for Big Data. ArXiv:2312.06465. https://doi.org/10.48550/arXiv.2312.06465

}
\author{
  Yajie Duan, Javier Cabrera
}
\seealso{
\code{\link[PPbigdata]{nuggKDE}},\code{\link[datanugget]{create.DN}}, \code{\link[datanugget]{refine.DN}}
}
\examples{

  require(datanugget)
  require(rstiefel)

  #4-dim small example
  X = cbind.data.frame(rnorm(5*10^3),
                       rnorm(5*10^3,2,1),
                       rnorm(5*10^3,5,2),
                       rnorm(5*10^3))

  #raw data is recommended to be scaled firstly to generate data nuggets for Projection Pursuit
  X = as.data.frame(scale(X))

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

  #spherize the data nuggets with weights to calculate the PP index
  nugg_wsph <- wsph(nugg,weight)$data_wsph

  #generate a random orthonormal matrix as a projection matrix to 2-dim space
  proj_2d = rustiefel(4, 2)

  #project data nugget centers into 2-dim space by the random projection matrix
  nuggproj_2d = as.matrix(nugg_wsph)\%*\%proj_2d

  #plot the projected data nuggets
  #lighter green represents larger weights
  plotNugg(nuggproj_2d, weight)

  #calculate the Natural Hermite index for the projected 2-dim big data
  NHnugg(nuggproj_2d,weight,scale)

}

