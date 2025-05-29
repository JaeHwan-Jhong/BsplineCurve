##########################################
##### ADMM for Bspline curve fitting #####
##########################################

# 24.07.28
#-----------------------------------------------------------------------------#
bspline.curve.admm_lambdas = function(y, D, B, 
                                      lambdas = NULL, 
                                      lam_max = 100, 
                                      lam_min = 1e-05, 
                                      n_lambda = 100, 
                                      max_iter = 500, 
                                      epsilon = 1e-5,
                                      eta_c = 0.1)
{
  results = list()
  
  # lam_1 > ... > lam_r
  if(is.null(lambdas))
    lambdas = exp(seq(log(lam_min), log(lam_max), length = n_lambda))
  n_lambda = length(lambdas)
  
  # fixed values
  n = nrow(as.matrix(y))
  m = ncol(as.matrix(y))
  Bty = t(B) %*% y
  BtB = t(B) %*% B
  DtD = t(D) %*% D
  #eta = 1
  #inv = solve((BtB + eta * DtD))
  #eta_D = eta * t(D)
  K = nrow(D)
  J = ncol(D)
  
  # initial values
  alpha = matrix(0, K, m)
  xi = matrix(0, J, m)
  u = matrix(0, K, m)
  aic = bic = rep(NA, length = n_lambda)
  
  for(k in 1:n_lambda)
  {
    lambda = lambdas[k]
     eta = lambda * eta_c
    lambda.eta = lambda / eta
     eta_D = eta * t(D)
     inv = solve(BtB + eta * DtD)
    rss_pen_store = Inf
    for(iter in 1:max_iter)
    {
      # xi update
      xi = inv %*% (Bty + eta_D %*% (alpha + u))
      Dxi = D %*% xi
      # alpha update
      alpha = BST(Dxi - u, lambda.eta)
      
      #alpha[k, ] = BST(D[k, ]%*%beta-v[k, ], lambda/xi)
      
      # z update
      u = u + (alpha - Dxi)
      # check convergence
      pen = lambda * sum(sqrt(rowSums(Dxi^2)))
      rss = 0.5 * sum((y - B %*% xi)^2)
      rss_pen = rss + pen
      if (abs(rss_pen - rss_pen_store) < epsilon)
        break;
      rss_pen_store = rss_pen
    }
    p = sum(rowSums(alpha) != 0)
    aic[k] = n * log(2 * rss / n) + p * 2
    bic[k] = n * log(2 * rss / n) + p * log(n)
    
    # cat("iteration = ", iteration, "\r")

    results[[k]] = list()
    results[[k]]$xi = xi
    results[[k]]$Bb = B %*% xi
    results[[k]]$alpha = alpha
    results[[k]]$eta = eta
    results[[k]]$lambda = lambda
    results[[k]]$u = u
    results[[k]]$iteration = iter
    results[[k]]$obj_f = rss_pen
    results[[k]]$aic = aic[k]
    results[[k]]$bic = bic[k]
  }
  
  best_aic = which.min(aic)
  best_bic = which.min(bic)
  
  results$best_aic = results[[best_aic]]
  results$best_bic = results[[best_bic]]
  results$aic = aic
  results$bic = bic
  results$r = n_lambda
  
  return(results)
}

# prox.lasso
# ST
ST = function(v, lambda)
{
  p = length(v)
  s = rep(0, p)
  for (i in 1 : p)
  {
    if (v[i] < -lambda)
      s[i] = v[i] + lambda
    else if (v[i] > +lambda)
      s[i] = v[i] - lambda
    else
      s[i] = 0
  }
  s
}
#-----------------------------------------------------------------------------#
## BST (Group Proximal) ##
BST = function(v, lambda)
{
   # v : K X M
   p = matrix(0, nrow(v), ncol(v))
   for(m in 1:nrow(v))
   {
      row.norm = sqrt(sum(v[m, ]^2))
      if(row.norm > lambda)
         p[m, ] = (1 - lambda / row.norm) * v[m, ]
   }
   return(p)
}
