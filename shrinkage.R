shrinkage = function (X, tau) {
  X = as.matrix(X)
  SVD = svd(X, nu = nrow(X), nv = ncol(X))
  U = SVD$u
  V = SVD$v
  D = SVD$d
  d1 = pmax(0, D - tau)
  d = U %*% diag(d1, nrow = nrow(X), ncol = ncol(X)) %*% t(V)
  return(d)
}
