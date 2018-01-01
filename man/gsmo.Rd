\name{gsmo}
\alias{gsmo}
\title{
General sequential minimal optimization
}
\description{
General sequential minimal optimization(GSMO) solves the quadratic programming problem: \cr
\eqn{min(f'*x + 1/2*x'*H*x)}\cr
subject to \cr 
\eqn{A*x = b} \cr
\eqn{l <= x <= u}
}
\usage{
gsmo(H, f, A, b, LB, UB, x0, MaxIter, TolKKT, verb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H}{
  Matrix appearing in the quadratic function to be minimized.
}
  \item{f}{
  Vector or one column matrix appearing in the quadratic function to be minimized.
  }
  \item{A}{
  Vector or one column matrix defining the constrains under which we minimize the quadratic function.
  }
  \item{b}{
  Scalar defining the constrains.
  }
  \item{LB}{
  Lower bound vector or one column matrix.
  }
  \item{UB}{
  Upper bound vector or one column matrix.
  }
  \item{x0}{
  Initial vector. Default is zero vector.
  }
  \item{MaxIter}{
  Maximum number of iterations.
  }
  \item{TolKKT}{
  Tightness of KKT stopping conditions. Default = 1e-9
  }
  \item{verb}{
  Display convergence information. Default = 0
  }
}
\details{
Sequential minimal optimization(SMO) is an algorithm for solvind the quadratic programming(QP) problem that arises during the training of support vector machines.
}
\value{
  An S3 object with the following slots
  \item{alpha}{Vector containing the primal solution of the quadratic problem.}
  \item{obj_value}{Primal objective value.}
  \item{exitflag}{Indicates which stopping condition was used: \cr 
                  -1 ... not enough memory \cr 
                   0 ... maximal number of iterations reached: nIter >= MaxIter \cr 
                   4 ... Relaxed KKT conditions satisfied.}
  \item{nIter}{Number of iterations.}
}
\references{
S.S. Keerthi, E.G. Gilbert. Convergence of a generalized SMO algorithm \cr
     for SVM classier design. Technical Report CD-00-01, Control Division, \cr 
     Dept. of Mechanical and Production Engineering, National University \cr
     of Singapore, 2000. \cr
     http://citeseer.ist.psu.edu/keerthi00convergence.html \cr
}
\author{
Park Beomjin, Moon Sangjun(based on Matlab code by Vojtech Franc)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(iris)
iris = subset(iris, subset = Species != "versicolor")
X = as.matrix(iris[, 1:4])
y = ifelse(iris$Species == "setosa", 1, -1)
C = 100

K = (y \%*\% t(y)) * (X \%*\% t(X))

# Dual solution for linear SVM
alpha = gsmo(H = K, f = -rep(1, nrow(K)), A = t(y), b = 0,
             LB = rep(0, nrow(K)), UB = rep(C, nrow(K)), MaxIter = 1e+5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
