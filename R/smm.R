smm_core = function(X, y, p, q, C, tau, max_iter, inner_iter, eps, rho, eta)
{
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
   
    alpha = gsmo(H = K, f = -f, A = t(y), b = 0, 
                LB = LB, UB = UB, x0 = rep(0, ncol(K)), MaxIter = inner_iter, TolKKT = eps/100)$alpha
    
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

shrinkage = function(X, tau) 
{
  X = as.matrix(X)
  SVD = svd(X, nu = nrow(X), nv = ncol(X))
  U = SVD$u
  V = SVD$v
  D = SVD$d
  d1 = pmax(0, D - tau)
  d = U %*% diag(d1, nrow = nrow(X), ncol = ncol(X)) %*% t(V)
  return(d)
}


smm = function(X, ...) {UseMethod("smm")}

smm.default = function(X, y, p, q, C, tau, max_iter = 1e+4, inner_iter = 1e+5, eps = 1e-6, rho = 20, eta = 0.999, ...) 
{
  X = as.matrix(X)
  y = as.numeric(y)
  
  if(ncol(X) != (p * q)) {warning("ncol(X) != (p * q)")}
  
  smm_sol = smm_core(X, y, p, q, C, tau, max_iter, inner_iter, eps, rho, eta)
  smm_sol$fitted.values = as.vector(X %*% smm_sol$w_k + smm_sol$b)
  smm_sol$call = match.call()
  
  class(smm_sol) = "smm"
  
  return(smm_sol)
}

smm.formula = function(formula, data = list(), p, q, C, tau, max_iter = 1e+4, inner_iter = 1e+5, eps = 1e-6, rho = 20, eta = 0.999, ...) 
{
  mf = model.frame(formula = formula, data = data)
  X = as.matrix(mf[, -1])
  y = mf[, 1]
  
  if(ncol(X) != (p * q)) {warning("ncol(X) != (p * q)")}
    
  smm_sol = smm.default(X, y, p, q, C, tau, max_iter, inner_iter, eps, rho, eta)
  smm_sol$fitted.values = as.vector(X %*% smm_sol$w_k + smm_sol$b)
  smm_sol$call = match.call()
  smm_sol$formula = formula
  
  return(smm_sol)
}

predict.smm = function(object, newdata = NULL, type = "class", ncores = 1, ...)
{
  if(is.null(newdata)) {
    y = fitted(object)
  } else {
    X = newdata
    y = X %*% object$w_k + object$b
  }
  
  if (type == "class") {res = ifelse(y > 0, 1, -1)} else {res = y}
  
  return(res)
}

kfold_cv = function(X, y, p, q, cost_range, tau_range, nfolds, optModel = TRUE, ncores = 1, ...)
{
  params = expand.grid(cost = cost_range, tau = tau_range)
  
  fold_list = createFolds(y, k = nfolds, list = FALSE)
  valid_err_mat = matrix(NA, nrow = nfolds, ncol = nrow(params))
  
  for (i in 1:nfolds) {
    cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
    fold = which(fold_list == i)
    y_fold = y[-fold]
    x_fold = X[-fold, , drop = FALSE]
    y_valid = y[fold]
    x_valid = X[fold, , drop = FALSE]
    
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          smm_fit = smm.default(X = x_fold, y = y_fold, p, q, C = params$cost[j], tau = params$tau[j],
                                                ...)
                          pred_val = predict.smm(smm_fit, newdata = x_valid)
                          acc = sum(y_valid == pred_val) / length(y_valid)
                          return(acc)
                        }, mc.cores = ncores)
    
    valid_err_mat[i, ] = sapply(fold_err, "[[", 1)
  }
  valid_err = colMeans(valid_err_mat, na.rm = TRUE)
  opt_ind = max(which(valid_err == min(valid_err)))
  opt_param = params[opt_ind, ]
  opt_valid_err = min(valid_err)
  
  out = list()
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  
  if (optModel) {
    opt_model = smm.default(X = X, y = y, p, q, C = opt_param$cost, tau = opt_param$tau, ...)
    out$opt_model = opt_model
  }
  return(out)
}
