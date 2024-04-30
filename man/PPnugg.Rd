\name{PPnugg}
\alias{PPnugg}
\title{
Projection Pursuit for Big Data based on data nuggets
}
\description{
This function performs 1-dim/2-dim projection pursuit (PP) for big data based on data nuggets.
}
\usage{
PPnugg(data, index = c("NH","Hole","CM"), dim, h = NULL, fa = TRUE, den = TRUE,
       R = 5000, DN.num1 = 10^4, DN.num2 = 2000, max.splits = 5, seed_nugg = 5,
       cooling = 0.9, tempMin = 1e-3, maxiter = 2000, tol = 1e-6, seed_opt = 3,
       initP = NULL,qt = 0.8,label = colnames(data), ...)
}
\arguments{
  \item{data}{A big data matrix (of class matrix, or data.frame) containing only entries of class numeric. The data size is large in terms of the number of observations, i.e., nrow(data).
  }
  \item{index}{A character indicating the PP index function to be optimized:
  "NH" - Natural Hermite Index for data nuggets
  "Hole" - Hole Index for data nuggets
  "CM" - Central Mass Index for data nuggets
  }
  \item{dim}{A numerical value indicating the target dimensionality for the projection. It's either 1 or 2.}
  \item{h}{If \code{index == "NH"}, a scalar or a length-2 numeric vector indicating the bandwidth used in the calculation of Natural Hermite index for data nuggets. Defaults to NULL, which uses normal reference bandwidth considering the nugget weights. See details in \code{\link[PPbigdata]{NHnugg}}.
  }
  \item{fa}{Logical value indicating whether to perform factor rotation after obtaining the optimal projection. Default is TRUE. See details in \code{\link[PPbigdata]{faProj}}.
  }
  \item{den}{Logical value indicating whether to estimate the density function of the optimal projected data based on data nuggets. Default is TRUE. See details in \code{\link[PPbigdata]{nuggKDE}}.
  }
  \item{R}{The number of observations to sample from the data matrix when creating the initial data nugget centers. Must be of class numeric within [100,10000]. Default is 5000. See details in \code{\link[datanugget]{create.DN}}.
  }
  \item{DN.num1}{The number of initial data nugget centers to create. Must be of class numeric. Default is 10^4. See details in \code{\link[datanugget]{create.DN}}.
  }
  \item{DN.num2}{The number of data nuggets to create. Must be of class numeric. Default is 2000. See details in \code{\link[datanugget]{create.DN}}.
  }
  \item{max.splits}{A numeric value indicating the maximum amount of attempts that will be made to split data nuggets during the refining of data nuggets. Default is 5. See details in \code{\link[datanugget]{refine.DN}}.
  }
  \item{seed_nugg}{A numeric value indicating the random seed for replication of data nugget creation and refining. Default is 5. See details in \code{\link[datanugget]{create.DN}}, and \code{\link[datanugget]{refine.DN}}.
  }
  \item{cooling}{The cooling factor for optimization of the PP index. Default is 0.9. See details in \code{\link[PPbigdata]{PPnuggOptim}}.
  }
  \item{tempMin}{The minimal temperature parameter to stop the optimization of the PP index. It's a numeric value between (0,1). Default is 1e-3. See details in \code{\link[PPbigdata]{PPnuggOptim}}.
  }
  \item{maxiter}{The maximal number of iterations for optimization of the PP index. Default is 2000. See details in \code{\link[PPbigdata]{PPnuggOptim}}.
  }
  \item{tol}{The tolerance parameter for the PP index value during optimization. Default is 1e-6. See details in \code{\link[PPbigdata]{PPnuggOptim}}.
  }
  \item{seed_opt}{A positive integer value indicating the seed for optimization of the PP index. Default is 3.
  }
  \item{initP}{The initial projection matrix to start the optimization of PP index. Must be a data matrix (of class matrix, or data.frame) with the dimension of ncol(data)*dim. Default is NULL, which uses a matrix of orthogonal initialization. See details in \code{\link[PPbigdata]{PPnuggOptim}}.
  }
  \item{qt}{A scalar with value in \code{[0,1]} used in the plot of projected data nuggets. See details in \code{\link[PPbigdata]{plotNugg}}.
  }
  \item{label}{A character vector specifying the text to be written for the variables on the loading plot. Defaults to the column names of the bia data set. See details in \code{\link[PPbigdata]{plotLoadings}}.
  }
  \item{...}{Other arguments sent to \code{\link[datanugget]{create.DN}}, \code{\link[datanugget]{refine.DN}}, \code{\link[PPbigdata]{PPnuggOptim}}, \code{\link[PPbigdata]{NHnugg}}, \code{\link[PPbigdata]{faProj}}, \code{\link[PPbigdata]{plotNugg}}, \code{\link[PPbigdata]{plotLoadings}}, and \code{\link[PPbigdata]{nuggKDE}}.
}
}
\details{
  This function performs 1-dim/2-dim projection pursuit (PP) for big data based on data nuggets.

  Projection Pursuit (PP) is a tool for high-dim data to find low-dim projections indicating hidden structures such as clusters, outliers, and other non-linear structures. The interesting low dimensional projections of high-dimensional data could be found by optimizing PP index function. PP Index function numerically measures features of low-dimensional projections. Higher values of PP indices correspond to more interesting structure, such as point mass, holes, clusters, and other non-linear structures.

  Data nuggets are a representative sample meant to summarize Big Data by reducing a large dataset to a much smaller dataset by eliminating redundant points while also preserving the peripheries of the dataset. Each data nugget is defined by a center (location), weight (importance), and scale (internal variability). Data nuggets for a large dataset could be created and refined by functions \code{create.DN} or \code{refine.DN} in the package \code{datanugget}.

  The projection pursuit for big data with a large number of observations could be performed based on data nuggets. Before creating data nuggets for big data, the raw data is recommended to be standardized first for Projection Pursuit, which is performed by \code{\link{scale}} in the function. After obtaining created and refined data nuggets for big data, data nugget centers needs to be spherized considering nugget weights before conducting projection pursuit. The optimal or interested projection found by projection pursuit would be taken on the spherized nugget centers. The optimization of PP index is performed by \code{\link[PPbigdata]{PPnuggOptim}} for 1-dim/2-dim projection using grand tour simulated annealing method. The PP index for data nuggets could be Natural Hermite index, Hole index and CM index. See details in \code{\link[PPbigdata]{NHnugg}}, \code{\link[PPbigdata]{HoleNugg}}, and \code{\link[PPbigdata]{CMNugg}}.

  After obtaining the optimal 1-dim/2-dim projection by maximizing PP index for data nuggets, a factor analysis could be performed on the projection to back to the original variables. See details in \code{\link[PPbigdata]{faProj}}. The rotation is taken on the overall transformation matrix for the raw data nuggets, which is a combination of spherization matrix and optimial projection matrix, to back to the original variables.

  The density of projected data can be estimated, which is performed by \code{\link[PPbigdata]{nuggKDE}} in the function. Moreover, the same centering, spherization and the optimal projection found for the data nuggets is also performed on the standardized raw data to output the projected raw big data found by Projection Pursuit based on data nuggets, i.e., \code{dataproj}. Based on the 1-dim/2-dim projection found, the hidden structures such as clusters, outliers, and other non-linear structures inside the big data set could be explored.
}
\value{
A list containing the following components:

\item{DN}{An object of class datanugget obtained from the standardized raw big data after data nugget creation and refining. It contains a data frame including data nugget centers, weights, and scales, and a vector of length nrow(data) indicating the data nugget assignments of each observation in the raw data.}

\item{nuggproj}{If \code{fa == TRUE}, the rotated projected data nugget centers after conducting factor ratation on the optimal projection found. If \code{fa == FALSE}, the projected data nugget centers under the optimal projection found.}

\item{loadings}{A matrix of loadings for original variables, one column for each projection direction. If \code{fa == TRUE}, it's the rotated transformation matrix to obtain the projected data nugget centers after factor rotation. If \code{fa == FALSE}, it's an overall transformation matrix for the centered data nuggets, which is a combination of spherization matrix and the optimal projection matrix found, to back to the original variables. In either case, \code{nuggproj} can be obtained by multiplying the centered data nuggets \code{nugg_wcen} with this loading matrix \code{loadings}.
}

\item{nugg_wcen}{The centered data nugget centers that has a zero weighted mean for each column considering nugget weights. It's obtained by extracting the weighted mean from the original data nugget centers.}

\item{dataproj}{The corresponding projected raw data. It's obtained by performing the same projection found on the standardized raw big data, i.e., \code{dataproj = data_cen \%*\% loadings.}}

\item{data_cen}{The standardized raw data centered by the weighted mean of data nugget centers.}

\item{index.opt}{The optimal PP index value found.}

\item{density}{If \code{den == TRUE}, a list containing the estimated density for projected data. See details in \code{\link[PPbigdata]{nuggKDE}}.}

}
\references{
  Cook, D., Buja, A., & Cabrera, J. (1993). Projection pursuit indexes based on orthonormal function expansions. Journal of Computational and Graphical Statistics, 2(3), 225-250.

  Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.

  Duan, Y., Cabrera, J., & Emir, B. (2023). A New Projection Pursuit Index for Big Data. ArXiv:2312.06465. https://doi.org/10.48550/arXiv.2312.06465

  Cabrera, J., & McDougall, A. (2002). Statistical consulting. Springer Science & Business Media.

  Horst, P. (1965). Factor Analysis of Data Matrices. Holt, Rinehart and Winston. Chapter 10.

  Kaiser, H. F. (1958). The varimax criterion for analytic rotation in factor analysis. Psychometrika, 23(3), 187-200.

}
\author{
  Yajie Duan, Javier Cabrera
}
\seealso{
\code{\link[PPbigdata]{PPnuggOptim}}, \code{\link[PPbigdata]{NHnugg}},  \code{\link[PPbigdata]{HoleNugg}}, \code{\link[datanugget]{create.DN}}, \code{\link[datanugget]{refine.DN}}
}
\examples{
\donttest{

  require(datanugget)

  #4-dim small example with cluster stuctures in V3 and V4
  X = cbind.data.frame(V1 = rnorm(5*10^4,mean = 5,sd = 2),
                       V2 = rnorm(5*10^4,mean = 5,sd = 1),
                       V3 = c(rnorm(3*10^4,sd = 0.3),
                              rnorm(2*10^4,mean = 2, sd = 0.3)),
                       V4 = c(rnorm(1*10^4,mean = -8, sd = 1),
                              rnorm(3*10^4,mean = 0,sd = 1),
                              rnorm(1*10^4,mean = 7, sd = 1.5)))


  #perform 2-dim Projection Pursuit for the big data
  #based on Hole index for data nuggets
  res = PPnugg(X, index = "Hole", dim = 2, R = 5000, DN.num1 = 1*10^4, DN.num2 = 2000,
  no.cores = 2, tempMin = 0.05, maxiter = 1000, tol = 1e-4)

  #data nuggets created and refined from the standardized raw data
  nugg = res$DN$`Data Nuggets`

  #data nugget assignments of each observation in the raw data
  nugg_assign = res$DN$`Data Nugget Assignments`

  #plot projected data nuggets
  plotNugg(res$nuggproj,nugg$Weight,qt = 0.8)

  #plot the corresponding projected raw big data
  plot(res$dataproj,cex = 0.5,main = "Projected Raw Data")

  #plot the estimated density of the projected data
  image(res$density)

  #plot loadings of original variables
  #V3 and V4 have large loadings, same as the simulation setting.
  plotLoadings(res$loadings)



  #perform 1-dim Projection Pursuit for the big data
  #based on Natural Hermite index for data nuggets
  res = PPnugg(X, index = "NH", dim = 1, R = 5000, DN.num1 = 1*10^4, DN.num2 = 2000,
  no.cores = 2, tempMin = 0.05, maxiter = 1000, tol = 1e-5)

  #data nuggets created and refined from the standardized raw data
  nugg = res$DN$`Data Nuggets`

  #data nugget assignments of each observation in the raw data
  nugg_assign = res$DN$`Data Nugget Assignments`

  #plot projected data nuggets
  plotNugg(res$nuggproj,nugg$Weight,qt = 0.8,hist = TRUE)

  #plot the corresponding projected raw big data
  hist(res$dataproj,breaks = 100)
}

}

