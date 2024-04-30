\name{PPnuggOptim}
\alias{PPnuggOptim}
\title{
Optimize Projection Pursuit index for Big Data based on data nuggets
}
\description{
Optimize PP index for 1-dim/2-dim projection for big data based on data nuggets, using grand tour simulated annealing optimizationg method.
}
\usage{
PPnuggOptim(FUN, nugg_wsph, dimproj, tempInit = 1, cooling = 0.9, eps = 1e-3,
            tempMin = 0.01, maxiter = 1000, half = 10, tol = 1e-5, maxc = 15,
            seed = 3, initP = NULL, ...)
}
\arguments{
  \item{FUN}{The index function for data nuggets to optimize.
  }
  \item{nugg_wsph}{The data nugget centers spherized with nugget weights. Must be a data matrix (of class matrix, or data.frame) with at least two columns. See \code{\link[PPbigdata]{wsph}} for spherization with nugget weights.
  }
  \item{dimproj}{A numerical value indicating the dimension of the data projection. It's either 1 or 2.
  }
  \item{tempInit}{The initial temperature parameter for optimization, i.e., grand tour simulated annealing method. It's a numeric value between (0,1).
  }
  \item{cooling}{The cooling factor for optimization. Rate by which the temperature is reduced from one cycle to the next, i.e., the new temperature = cooling * old temperature. Default is 0.9.
  }
  \item{eps}{The approximation accuracy for cooling used in the interpolation between current projection and the target one. Default is 1e-3.
  }
  \item{tempMin}{The minimal temperature parameter to stop the optimization. It's a numeric value between (0,1). Default is 0.01.
  }
  \item{maxiter}{The maximal number of iterations for the optimization. Default is 1000.
  }
  \item{half}{The number of steps without incrementing the index before decreasing the temperature parameter by multiplying the cooling factor. Default is 10.
  }
  \item{tol}{The tolerance parameter for the index value during optimization. The increase of index value smaller than the tolerance would be ignored. Default is 1e-5.
  }
  \item{maxc}{The maximal number of temperature changes without incrementing the index value. The algorithm would stop if the number of temperature changes without index value increasing exceeds maxc paramter. Default is 15.
  }
  \item{seed}{A positive integer value indicating the seed for the optimization. Default is 3.
  }
  \item{initP}{The initial projection matrix to start the optimization. Must be a data matrix (of class matrix, or data.frame) with the dimension of ncol(nugg_wsph)*dimproj. Default is NULL, which uses a matrix of orthogonal initialization.
  }
  \item{...}{Other arguments passed to the index function FUN.
  }
}
\details{
  This function performs the optimization PP index for 1-dim/2-dim projection for big data based on data nuggets, using grand tour simulated annealing optimizationg method.

  Data nuggets are a representative sample meant to summarize Big Data by reducing a large dataset to a much smaller dataset by eliminating redundant points while also preserving the peripheries of the dataset. Each data nugget is defined by a center (location), weight (importance), and scale (internal variability). Data nuggets for a large dataset could be created and refined by functions \code{create.DN} or \code{refine.DN} in the package \code{datanugget}.

  After obtaining created and refined data nuggets for big data, data nugget centers needs to be spherized considering nugget weights before conducting projection pursuit. The optimal or most interested projection could be found by optimization of the projection pursuit index based on data nuggets. This function optimizes the PP index function for data nuggets by GTSA, i.e., the grand tour simulated annealing optimizationg method. The optimization starts with an initial projection and a initial temperature parameter. A target projection would be generated from the current projection plus the temperature times the initial base, from which a projection matrix is generated through interpolation between current and target projections. If the number of steps without incrementing the index exceeds the parameter half, the temperature parameter is decreased by multiplying the cooling factor. The increase of index value smaller than the tolerance would be ignored. The optimization would stop if either there are a maximal number of iterations, or the temperature has decreased to the minimal value, or there are a maximal number of temperature changes without incrementing the index value.
}
\value{
A list containing the following components:

\item{proj.nugg}{The projected data nugget centers under the optimal projection found.}

\item{proj.opt}{The optimal projection matrix found.}

\item{index}{A vector with the PP index values found in the optimization process.}

}
\references{
  Cook, D., Buja, A., & Cabrera, J. (1993). Projection pursuit indexes based on orthonormal function expansions. Journal of Computational and Graphical Statistics, 2(3), 225-250.

  Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.

  Duan, Y., Cabrera, J., & Emir, B. (2023). A New Projection Pursuit Index for Big Data. ArXiv:2312.06465. https://doi.org/10.48550/arXiv.2312.06465


}
\author{
  Yajie Duan, Javier Cabrera
}
\seealso{
\code{\link[PPbigdata]{PPnugg}}, \code{\link[PPbigdata]{NHnugg}},\code{\link[datanugget]{create.DN}}, \code{\link[datanugget]{refine.DN}}
}
\examples{
\donttest{

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
                    no.cores = 2,
                    make.pbs = FALSE)


  #refine data nuggets
  my.DN2 = refine.DN(x = X,
                     DN = my.DN,
                     EV.tol = .9,
                     min.nugget.size = 2,
                     max.splits = 5,
                     no.cores = 2,
                     make.pbs = FALSE)

  #get nugget centers, weights, and scales
  nugg = my.DN2$`Data Nuggets`[,2:(ncol(X)+1)]
  weight = my.DN2$`Data Nuggets`$Weight
  scale = my.DN2$`Data Nuggets`$Scale

  #spherize data nugget centers considering weights to conduct Projection Pursuit
  wsph.res = wsph(nugg,weight)
  nugg_wsph = wsph.res$data_wsph
  wsph_proj = wsph.res$wsph_proj

  #conduct Projection Pursuit in 2-dim by optimizing Natural Hermite index
  res = PPnuggOptim(NHnugg, nugg_wsph, dimproj = 2, weight = weight, scale = scale,
        tempMin = 0.05, maxiter = 1000, tol = 1e-5)

  #plot projected data nuggets
  plotNugg(nugg_wsph\%*\%res$proj.opt,weight,qt = 0.8)


  #conduct Projection Pursuit in 1-dim by optimizing Hole index
  res = PPnuggOptim(HoleNugg, nugg_wsph, dimproj = 1, weight = weight,
        tempMin = 0.05, maxiter = 1000, tol = 1e-5)

  #plot projected data nuggets
  plotNugg(nugg_wsph\%*\%res$proj.opt,weight)
}

}

