\name{HoleNugg}
\alias{HoleNugg}
\title{
Hole Index for projected Big Data based on data nuggets
}
\description{
This function calculates the value of Hole index, a Projection Pursuit index, for projected big data based on data nuggets.
}
\usage{
HoleNugg(nuggproj,weight)
}
\arguments{
\item{nuggproj}{Projected data nugget centers. Must be a data matrix (of class matrix, or data.frame) or a vector containing only entries of class numeric.
  }
  \item{weight}{Vector of the weight parameter for each data nugget. Its length should be the same as the number of data nuggets, i.e., nrow(nuggproj) or length(nuggproj). Must be of class numeric or integer.
  }
}
\details{
  This function calculates the value of Hole index, a Projection Pursuit index for projected Big Data based on data nuggets.

  Data nuggets are a representative sample meant to summarize Big Data by reducing a large dataset to a much smaller dataset by eliminating redundant points while also preserving the peripheries of the dataset. Each data nugget is defined by a center (location), weight (importance), and scale (internal variability). Data nuggets for a large dataset could be created and refined by functions \code{create.DN} or \code{refine.DN} in the package \code{datanugget}.

  Hole index is a kind of Projection Pursuit (PP) index, and larger index values indicate a hole structure of multivariate data. However, it's computationally hard to calculate the index for big data because of the vector memory limit during calculation. To deal with big data, data nuggets could be used to calculate the index efficiently. In this function, based on the projected data nugget centers with data nugget weights, the Hole index is calculated via a weighted version of the original Hole index formula.
}
\value{
  A numeric value indicating Hole index value of the projected big data based on the data nuggets.
}
\references{
  Chen, C. H., Hardle, W. K., & Unwin, A. (Eds.). (2007). Handbook of data visualization. Springer Science & Business Media.

  Cook, D., Buja, A., Cabrera, J., & Hurley, C. (1995). Grand tour and projection pursuit. Journal of Computational and Graphical Statistics, 4(3), 155-172.

  Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.

  Duan, Y., Cabrera, J., & Emir, B. (2023). A New Projection Pursuit Index for Big Data. ArXiv:2312.06465. https://doi.org/10.48550/arXiv.2312.06465

}
\author{
  Yajie Duan, Javier Cabrera
}
\seealso{
\code{\link[PPbigdata]{CMNugg}}, \code{\link[PPbigdata]{NHnugg}},\code{\link[datanugget]{create.DN}}, \code{\link[datanugget]{refine.DN}}
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

  #calculate the Hole index for the projected 2-dim big data
  HoleNugg(nuggproj_2d,weight)

}

