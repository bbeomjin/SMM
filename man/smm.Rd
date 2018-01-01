\name{smm}
\alias{smm}
\alias{smm.default}
\alias{smm.formula}
\alias{predict.smm}

\title{
Support Matrix Machine
}
\description{
smm is used to train a support matrix machine. It can be used to carry out classification, A formula interface is provided.
}
\usage{
smm(X, ...)

\method{smm}{default}(X, y, p, q, C, tau, max_iter, inner_iter, eps, rho, eta, ...)

\method{smm}{formula}(formula, data = list(), p, q, C, tau, max_iter, inner_iter, eps, rho, eta, ...)

\method{predict}{smm}(object, newdata = NULL, type = "class", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a symbolic description of the model to be fit.
  }
  \item{data}{
an optional data frame contatining the variables in the model. By default the variables are taken from the encironment which 'smm' is called from.
  }
  \item{X}{
  a data matrix or a sparse matrix
  }
  \item{y}{
  a response vector with one label for each row/component of X. It must be 1 or -1
  }
  \item{p}{
  the number of row of an input matrix. 
  }
  \item{q}{
  the number of column of an input matrix.
  }
  \item{C}{
  cost of constraints violation
  }
  \item{tau}{
  cost of regularization term. see "Detail"
  }
  \item{max_iter}{
  Maximum number of iterations for ADMM algorithm
  }
  \item{inner_iter}{
  Maximum number of iterations for SMO algorithm
  }
  \item{eps}{
  Tightness of stopping conditions
  }
  \item{rho}{
  ADMM algorithm parameter
  }
  \item{eta}{
  ADMM algorithm parameter
  }
  \item{object}{
  an object of class \code{"smm"}, i.e., a fitted model.
  }
  \item{newdata}{
  An object containing the new input data
  }
  \item{type}{
  one of "class", "raw" indicating the type of output
  }
  \item{\dots}{
  currently not used.
  }
}
\value{
An object of class 'smm'containing the fitted model, including:
  \item{w_k}{The corresponding coefficients.}
  \item{b}{The intercept.}
  \item{rk}{The rank of coefficients matrix.}
  \item{stop_iter}{The number of iteration the algorithm needed.}
}

\details{
This code is builded by regerring to the SMM paper. Detailed theory is included in the SMM paper.
}

\references{
  Luo, L., Xie, Y., Zhang, Z., and Li, W.-J. (2015). Support matrix machine. \emph{Proceedings of the 32nd International Conference on Machine Learning}, Lille, France.
}
\author{
Park Beomjin(based on Matlab code by Vojtech Franc)
}

\seealso{
\link{gsmo}
}
\examples{
require(SMM)
data(nottingham)

X = as.matrix(nottingham[, -1])
y = ifelse(nottingham[, 1] == 1, 1, -1)

smm_model = smm(X, y, 200, 200, 100, 1, 1e+4, 1e+5, 1e-6, 20, 0.999)
table(y, predict(smm_model, newdata = X))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
