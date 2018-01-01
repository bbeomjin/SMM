gsmo = function(H, f, A, b, LB, UB, x0, MaxIter, TolKKT, verb)
{
  H = as.matrix(H)
  f = as.numeric(f)
  A = as.numeric(A)
  
  if (missing(x0) == TRUE) {x0 = rep(0, ncol(H))}
  if (missing(TolKKT) == TRUE) {TolKKT = 1e-9}
  if (missing(verb) == TRUE) {verb = 0}
  
  alpha = .Call("gsmo", H, f, A, b, LB, UB, x0, MaxIter, TolKKT)
  
  if (verb == TRUE) {
    cat("Settings of GSMO solver \n")
    cat(paste("MaxIter : ", MaxIter, "\n", sep = ""))
    cat(paste("TolKKT : ", TolKKT, "\n", sep = ""))
    cat(paste("nVar : ", ncol(H), "\n", sep = ""))
  }
  
  return(alpha)
}

