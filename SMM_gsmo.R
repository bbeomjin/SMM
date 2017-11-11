# matrix multiplication with cpp
#####################################################################################################
#####################################################################################################
# using gurobi packages
dyn.load("~/beomfile/gsmo/libqp_gsmo_r.so")
source("~/beomfile/Rscript/shrinkage.R")


fastADMM_gsmo = function(X, y, p, q, C, tau, max_iter, inner_iter, eps, rho, eta) {
  if (class(X) != "matrix") {X = as.matrix(X)}
  K = tcrossprod(X) * tcrossprod(y) / (rho + 1)
  n = nrow(X)
  d = ncol(X)
  
  s_km1 = s_hatk = rep(0, d)
  lambda_km1 = lambda_hatk = rep(1, d)
  t_k = 1
  c_km1 = 0
  
  recent_number = 50
  recent_idx = 0
  obj_recent = rep(0, recent_number)
  
  LB = matrix(rep(0, n))
  UB = matrix(C * rep(1, n))
  
  obj = function(w, p, q, b, X, y, C, tau, svd_v) {
      obj = 0.5 * crossprod(w) + C * sum(pmax(0, 1 - y * (X %*% w + b))) + tau * sum(svd_v$d)
  }   
 
  for (k in 1:max_iter) {
    f = 1 - (tcrossprod(X, t(lambda_hatk + rho * s_hatk)) * y) / (rho + 1)
   
    alpha = .Call("gsmo", K, -f, t(y), 0, LB, UB, rep(0, ncol(K)), inner_iter, eps/100, 1)
    w_k = (lambda_hatk + rho * s_hatk + crossprod(X, (alpha * y))) / (rho + 1)
    
    sel = (alpha > 0) & (alpha < C)
    b = as.numeric(crossprod(sel, (y - tcrossprod(X, t(w_k)))) / sum(sel))
    
    W_k = array(w_k, dim = c(p, q))
    Lambda_k = array(lambda_hatk, dim = c(p, q))
    S = shrinkage(rho * W_k - Lambda_k, tau) / rho
    s_k = as.vector(S) 
    lambda_k = lambda_hatk - rho * (w_k - s_k)
    c_k = crossprod((lambda_k - lambda_hatk)) / rho + rho * crossprod((s_k - s_hatk))
    
    if (c_k < eta * c_km1) {
      t_kp1 = 0.5 * (1 + sqrt(1 + 4 * t_k^2))
      s_hatkp1 = s_k + (t_k - 1) / t_kp1 * (s_k - s_km1)
      lambda_hatkp1 = lambda_k + (t_k - 1) / t_kp1 * (lambda_k - lambda_km1)
      restart = FALSE 
    } else {
      t_kp1 = 1
      s_hatkp1 = s_km1
      lambda_hatkp1 = lambda_km1
      c_k = c_km1 / eta
      restart = TRUE
    }
    s_hatk = s_hatkp1
    lambda_hatk = lambda_hatkp1
    c_km1 = c_k
    s_km1 = s_k
    lambda_km1 = lambda_k
    t_k = t_kp1
    
    svd_v = svd(W_k, nu = 0, nv = 0)
    
    obj_k = obj(w_k, p, q, b, X, y, C, tau, svd_v)
    recent_idx = recent_idx + 1
    obj_recent[recent_idx] = obj_k
    rk = sum(svd_v$d > 1e-6)
    if (recent_idx == recent_number) {
      recent_idx = 0
    }
    if ((k %% 1000) == 0) {
      sprintf('k = %d, obj = %f, restart = %d, rank = %d', k, obj_k, restart, rk)
    }
    if (abs(obj_k - mean(obj_recent, na.rm = T)) / abs(mean(obj_recent, na.rm = T)) < eps & k > recent_number) break;
  }
  stop_iter = k
  return(list(w_k = w_k, b = b, rk = rk, stop_iter = stop_iter))
}
