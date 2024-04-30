\name{PPbigdata-package}
\alias{PPbigdata-package}

\title{
Projection Pursuit for Big Data Based on Data Nuggets
}

\description{
This package contains functions to perform 1-dim/2-dim projection pursuit for big data based on data nuggets. It includes PP indices for data nuggets, static PP for big data by optimization of PP index, grand tour and guided tour for big data based on data nuggets, and visialization functions for projections and variable loadings.
}

\author{
Yajie Duan, Javier Cabrera
}

\references{
  Duan, Y., Cabrera, J., & Emir, B. (2023). A New Projection Pursuit Index for Big Data. ArXiv:2312.06465. https://doi.org/10.48550/arXiv.2312.06465

  Beavers, T. E., Cheng, G., Duan, Y., Cabrera, J., Lubomirski, M., Amaratunga, D., & Teigler, J. E. (2024). Data Nuggets: A Method for Reducing Big Data While Preserving Data Structure. Journal of Computational and Graphical Statistics, (just-accepted), 1-21.

  Cherasia, K. E., Cabrera, J., Fernholz, L. T., & Fernholz, R. (2022). Data Nuggets in Supervised Learning. \emph{In Robust and Multivariate Statistical Methods: Festschrift in Honor of David E. Tyler} (pp. 429-449). Cham: Springer International Publishing.

  Cook, D., Buja, A., & Cabrera, J. (1993). Projection pursuit indexes based on orthonormal function expansions. Journal of Computational and Graphical Statistics, 2(3), 225-250.

  Cook, D., Buja, A., Cabrera, J., & Hurley, C. (1995). Grand tour and projection pursuit. Journal of Computational and Graphical Statistics, 4(3), 155-172.

  Wickham, H., Cook, D., Hofmann, H., & Buja, A. (2011). tourr: An R package for exploring multivariate data with projections. Journal of Statistical Software, 40, 1-18.

  Cabrera, J., & McDougall, A. (2002). Statistical consulting. Springer Science & Business Media.

  Horst, P. (1965). Factor Analysis of Data Matrices. Holt, Rinehart and Winston. Chapter 10.

  Kaiser, H. F. (1958). The varimax criterion for analytic rotation in factor analysis. Psychometrika, 23(3), 187-200.
}

\seealso{
\code{\link[PPbigdata]{PPnugg}}, \code{\link[PPbigdata]{NHnugg}},  \code{\link[PPbigdata]{HoleNugg}}, \code{\link[PPbigdata]{plotNugg}}, \code{\link[PPbigdata]{faProj}}, \code{\link[PPbigdata]{guidedTourNugg}}, \code{\link[PPbigdata]{grandTourNugg}}, \code{\link[datanugget]{create.DN}}, \code{\link[datanugget]{refine.DN}}
}
