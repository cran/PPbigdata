\name{faProj}
\alias{faProj}
\title{
Factor rotation for projected Big Data in multi-dimensional space based on data nuggets
}
\description{
This function performs the factor rotation for projected big data in multi-dimensional space based on data nuggets.
}
\usage{
faProj(nugg, weight, wsph_proj = NULL, proj, method = c("varimax","promax"))
}
\arguments{
  \item{nugg}{Data nugget centers obtained from raw data. Must be a data matrix (of class matrix, or data.frame) with at least two columns.
  }
  \item{weight}{Vector of the weight parameter for each data nugget. Its length should be the same as the number of data nuggets, i.e., nrow(nugg). Must be of class numeric or integer.
  }
  \item{wsph_proj}{Matrix of size ncol(nugg) by ncol(nugg). It's the sphering/whitening matrix considering nugget weights for the transformation. The projection is on the spherized data nugget centers considering weights, which is obtained by multiplying the centered data nuggets with weights by this sphering/whitening matrix. Default is NULL, which would be obtained by function \code{\link[PPbigdata]{wsph}}. Must be a data matrix containing only entries of class numeric.
  }
  \item{proj}{Matrix of size ncol(nugg) by projection dimenstion. It's the orthonormal projection matrix that would be taken on the spherized data nugget centers considering weights, to obtain projected data nuggets. Must be a data matrix containing only entries of class numeric.
  }
  \item{method}{A character indicating the rotation method used for factor analysis. The default method "\code{varimax}" uses function \code{\link[stats]{varimax}} to take rotation; the alternative "\code{promax}" uses function \code{\link[stats]{promax}}. The rotation is taken on the overall transformation matrix for the raw data nuggets, which is a combination of spherization matrix and projection matrix, to back to the original variables.
  }
}
\details{
  This function performs the factor rotation for projected big data in multi-dimensional space based on data nuggets.

  Data nuggets are a representative sample meant to summarize Big Data by reducing a large dataset to a much smaller dataset by eliminating redundant points while also preserving the peripheries of the dataset. Each data nugget is defined by a center (location), weight (importance), and scale (internal variability). Data nuggets for a large dataset could be created and refined by functions \code{create.DN} or \code{refine.DN} in the package \code{datanugget}.

  After obtaining created and refined data nuggets for big data, data nugget centers needs to be spherized considering nugget weights before conducting projection pursuit. The optimal or interested projection found by projection pursuit would be taken on the spherized nugget centers. This function conducts the factor analysis for the projected data nugget centers.  The default rotation method "\code{varimax}" uses function \code{\link[stats]{varimax}} to take rotation; the alternative "\code{promax}" uses function \code{\link[stats]{promax}}. The rotation is taken on the overall transformation matrix for the raw data nuggets, which is a combination of spherization matrix and projection matrix, to back to the original variables.
}
\value{
A list containing the following components:

\item{nuggproj_rotat}{The rotated projected data nugget centers after conducting factor ratation. It's obtained by multiplying the centered data nuggets \code{nugg_wcen} with the rotated matrix \code{loadings}.}

\item{loadings}{A matrix of loadings for original variables, one column for each projection direction. It's the rotated transformation matrix to obtain updated projected data nugget centers.}

\item{nugg_wcen}{The centered data nugget centers that has a zero weighted mean for each column considering nugget weights. It's obtained by extracting the weighted mean from the original data nugget centers.}

}
\references{
  Cook, D., Buja, A., & Cabrera, J. (1993). Projection pursuit indexes based on orthonormal function expansions. Journal of Computational and Graphical Statistics, 2(3), 225-250.

  Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.

  Duan, Y., Cabrera, J., & Emir, B. (2023). A New Projection Pursuit Index for Big Data. ArXiv:2312.06465. https://doi.org/10.48550/arXiv.2312.06465

  Hendrickson, A. E., & White, P. O. (1964). Promax: A quick method for rotation to oblique simple structure. British journal of statistical psychology, 17(1), 65-70.

  Horst, P. (1965). Factor Analysis of Data Matrices. Holt, Rinehart and Winston. Chapter 10.

  Kaiser, H. F. (1958). The varimax criterion for analytic rotation in factor analysis. Psychometrika, 23(3), 187-200.

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
  require(rstiefel)

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

  #spherize data nugget centers considering weightsn to conduct Projection Pursuit
  wsph.res = wsph(nugg,weight)
  nugg_wsph = wsph.res$data_wsph
  wsph_proj = wsph.res$wsph_proj

  #conduct the same spherization projection on the standardized raw data
  X_cen = X- as.matrix(rep(1,nrow(X)))\%*\%wsph.res$wmean
  X_sph = as.matrix(X_cen)\%*\%wsph_proj

  #conduct Projection Pursuit in 2-dim by optimizing Natural Hermite index
  res = PPnuggOptim(NHnugg, nugg_wsph, dimproj = 2, weight = weight, scale = scale)

  #optimal projection matrix obtained
  proj_opt = res$proj.opt

  #plot projected data nuggets
  plotNugg(nugg_wsph\%*\%proj_opt,weight,qt = 0.8)

  #conduct varimax rotation for projection
  fa = faProj(nugg,weight,proj = proj_opt)

  #obtain rotated projected data nuggets and
  #corresponding loadings of original variables
  nuggproj_rotat = fa$nuggproj_rotat
  loadings = fa$loadings

  #plot rotated projected data nuggets after varimax rotation
  plotNugg(nuggproj_rotat,weight,qt = 0.8)

  #plot corresponding projected raw big data after factor roation
  X_proj = as.matrix(X_cen)\%*\%loadings
  plot(X_proj,cex = 0.5)

  #plot loadings of original variables
  #V3 and V4 have large loadings, same as the simulation setting.
  plotLoadings(loadings)
}
}

