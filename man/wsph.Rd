\name{wsph}
\alias{wsph}
\title{
Spherize/Whiten data with observational weights
}
\description{
This function performs PCA sphering/whitening transformation on data with observational weights.
}
\usage{
wsph(data,weight)
}
\arguments{

  \item{data}{
A data matrix (of class matrix, or data.frame) containing only entries of class numeric.
}

  \item{weight}{Vector of length nrow(data) of weights for each observation in the dataset. Must be of class numeric or integer or table. If NULL, the default value is a vector of 1 with length nrow(data), i.e., weights equal 1 for all observations.
}
}
\details{
This function performs PCA sphering/whitening transformation on data with observational weights. Specifically, weighted sample mean and weighted sample covariance matrix are firstly calculated. Next, data are centered with weights, and spectral decomposition of the weighted covariance matrix is conducted to obtain its eigenvalues and eigenvectors. Based on them, the PCA sphering/whitening transformation is performed to obtain a spherized data matrix considering the observational weights. The spherized data matrix has a zero weighted mean for each column, and a weighted covariance that equals identity matrix.
}
\value{
A list containing the following components:

\item{data_wsph}{The spherized data matrix that has a zero weighted mean for each column, and a weighted covariance that equals identity matrix.}

\item{data_wcen}{The centered data matrix that has a zero weighted mean for each column. It's obtained by extracting the weighted sample mean from the original data matrix.}

\item{wmean}{Vector of length ncol(data). It's the weighted sample mean of the original data matrix.}

\item{wcov}{Matrix of size ncol(data) by ncol(data). It's the weighted sample covariance matrix of the original data matrix.}

\item{wsph_proj}{Matrix of size ncol(data) by ncol(data). It's the sphering/whitening matrix for the transformation. The spherized data matrix \code{data_wsph} is obtained by multiplying the centered data matrix with weights \code{data_wcen} with this sphering/whitening matrix. }

}
\references{
Cabrera, J., & McDougall, A. (2002). Statistical consulting. Springer Science & Business Media.
}
\author{
Yajie Duan, Javier Cabrera
}

\examples{

 dataset = matrix(rnorm(300),100,3)

 # assign random weights to observations
 weight = sample(1:20,100,replace = TRUE)

 # spherize the dataset with observational weights
 res = wsph(dataset,weight)

 # spherized data matrix considering the observation weights
 res$data_wsph

 # The spherzied data matrix has a zero weighted sample mean for each column
 1/sum(weight)*t(as.matrix(weight))\%*\%as.matrix(res$data_wsph)

 # The spherzied data matrix has a weighted covariance that equals identity matrix
 1/sum(weight)*t(as.matrix(res$data_wsph))\%*\%diag(weight)\%*\%as.matrix(res$data_wsph)

}

