\name{guidedTourNugg}
\alias{guidedTourNugg}
\title{
1-dim/2-dim Guided Tour for Big Data based on Data Nuggets
}
\description{
This function performs a 1-dim/2-dim guided tour path for big data based on constructed data nuggets with their weights and scales. The guided tour tries to find a projection  with a higher value of PP index than the current projection.
}
\usage{
guidedTourNugg(nugg, weight, scale, dim, index = c("NH","Hole","CM"), qt = 0.8,...)
}
\arguments{
  \item{nugg}{Data nugget centers obtained from raw data. Must be a data matrix (of class matrix, or data.frame) with at least two columns.
  }
  \item{weight}{Vector of the weight parameter for each data nugget. Its length should be the same as the number of data nuggets, i.e., nrow(nugg). Must be of class numeric or integer.
  }
  \item{scale}{Vector of the scale parameter for each data nugget. Its length should be the same as the number of data nuggets, i.e., nrow(nuggproj) for 2-dim projection/length(nuggproj) for 1-dim projection. Must be of class numeric or integer.
  }
  \item{dim}{A numerical value indicating the target dimensionality for the tour. It's either 1 or 2.}
  \item{index}{A character indicating the PP index function to be used to guide the tour:
  "NH" - Natural Hermite Index for data nuggets
  "Hole" - Hole Index for data nuggets
  "CM" - Central Mass Index for data nuggets}
  \item{qt}{For projected plots of 2-dim tour, a scalar with value in \code{[0,1]} indicating the probability used to obtain a sample quantile of the data nugget weights as the maximal value to transform weights to colors for the plot. Defaults to be 0.8. See \code{\link[PPbigdata]{plotNugg}}.
  }
  \item{...}{Other arguments sent to \code{\link[PPbigdata]{NHnugg}}, \code{\link[tourr]{guided_tour}}, \code{\link[tourr]{animate}}, \code{\link[tourr]{animate_xy}}, \code{\link[tourr]{display_dist}}, \code{\link[tourr]{display_xy}}, \code{\link[weights]{wtd.hist}}.
}
}
\details{
  This function performs a 1-dim/2-dim guided tour path for big data based on constructed data nuggets with their weights and scales. The guided tour tries to find a projection  with a higher value of PP index than the current projection.

  Data nuggets are a representative sample meant to summarize Big Data by reducing a large dataset to a much smaller dataset by eliminating redundant points while also preserving the peripheries of the dataset. Each data nugget is defined by a center (location), weight (importance), and scale (internal variability). Data nuggets for a large dataset could be created and refined by functions \code{create.DN} or \code{refine.DN} in the package \code{datanugget}.

  Based on the data nuggets from big data, a projection pursuit guided tour is performed to explore the multivariate data. Unlike walking randomly to discover projections in grand tour, the guided tour selects the next target basis by optimizing a projection pursuit index function defining interesting projections. Here the considered choices of PP indices include Nature Hermite index, Hole index and CM index for big data based on data nuggets. See details in \code{\link[PPbigdata]{NHnugg}}, \code{\link[PPbigdata]{HoleNugg}}, and \code{\link[PPbigdata]{CMNugg}}.

  This function for data nuggets is based on functions about guided tour in the package \code{tourr}. See details in \code{\link[tourr]{guided_tour}}, \code{\link[tourr]{animate}}, \code{\link[tourr]{animate_xy}}, \code{\link[tourr]{display_dist}}, and \code{\link[tourr]{display_xy}}. For 2-dim guided tour, the projected data nugget centers are plotted with colors based on their weights where lighter green represents larger weights. For 1-dim guided tour, a weighted density histgram of 1-dim projected data nugget centers is plotted considering the data nugget weights. See details in \code{\link[weights]{wtd.hist}}. The loadings of each variable for projections are also shown at each step.
}
\value{
A list containing the bases, index values, and other information during the tour.
}
\references{
  Cook, D., Buja, A., Cabrera, J., & Hurley, C. (1995). Grand tour and projection pursuit. Journal of Computational and Graphical Statistics, 4(3), 155-172.

  Cook, D., Buja, A., & Cabrera, J. (1993). Projection pursuit indexes based on orthonormal function expansions. Journal of Computational and Graphical Statistics, 2(3), 225-250.

  Wickham, H., Cook, D., Hofmann, H., & Buja, A. (2011). tourr: An R package for exploring multivariate data with projections. Journal of Statistical Software, 40, 1-18.

  Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.

  Duan, Y., Cabrera, J., & Emir, B. (2023). A New Projection Pursuit Index for Big Data. ArXiv:2312.06465. https://doi.org/10.48550/arXiv.2312.06465

}
\author{
  Yajie Duan, Javier Cabrera
}
\seealso{
\code{\link[PPbigdata]{PPnugg}}, \code{\link[PPbigdata]{NHnugg}},\code{\link[datanugget]{create.DN}}, \code{\link[datanugget]{refine.DN}}, \code{\link[tourr]{guided_tour}}, \code{\link[tourr]{animate}}
}
\examples{

  require(datanugget)

  #4-dim small example with cluster stuctures in V3 and V4
  X = cbind.data.frame(V1 = rnorm(5*10^3,mean = 5,sd = 2),
                       V2 = rnorm(5*10^3,mean = 5,sd = 1),
                       V3 = c(rnorm(3*10^3,sd = 0.3),
                              rnorm(2*10^3,mean = 2, sd = 0.3)),
                       V4 = c(rnorm(1*10^3,mean = -8, sd = 1),
                              rnorm(3*10^3,mean = 0,sd = 1),
                              rnorm(1*10^3,mean = 7, sd = 1.5)))

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

  #2-dim guided tour by Natural Hermite Index based on data nuggets
  guidedTourNugg(nugg,weight,scale,dim = 2,index = "NH",cex = 0.5,max.tries = 15)

  #1-dim guided tour by Hole Index based on data nuggets
  guidedTourNugg(nugg,weight,scale,dim = 1,index = "Hole",density_max = 4.5)

}

