kernelMatrix = function(x, y, kernel = "gaussian", kparam = 1.0, kparam2 = 1) {
  
  x = as.matrix(x)
  y = as.matrix(y)
  p = ncol(x)
  
  if (NCOL(x) == 0) {
    x = matrix(0, nrow = nrow(x), ncol = 1)
  }
  
  if (NCOL(y) == 0) {
    y = matrix(0, nrow = nrow(y), ncol = 1)
  }
  
  if (kernel == "poly") {
    K = (x %*% t(y) + kparam2)^kparam
  } else if(kernel == "gaussian" | kernel == "gaussian2") {
    normx = rowSums(x^2)
    normy = rowSums(y^2)
    temp = x %*% t(y)
    temp = (-2.0 * temp) + outer(normx, rep(1.0, nrow(y)), "*") + outer(rep(1.0, nrow(x)), normy, "*")
    K = exp(-temp * kparam)
    # obj = kernelMatrix(rbfdot(sigma = kparam), x, y)
  } else if (kernel == "spline") {
    K = 0
    for (d in 1:p) {
      K_temp = spline_kernel(x[, d, drop = FALSE], y[, d, drop = FALSE])
      K = K + K_temp$K1 + K_temp$K2
    }
  } else if (kernel == "linear") {
    K = tcrossprod(x, y)
  } else if (kernel == "anova_gaussian") {
    K = 0
    for (d in 1:p) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "gaussian", kparam = kparam)
      K = K + K_temp
    }
  } else {
    K = NULL
  }
  return(K)
}


spline_kernel = function(x, y)
{
  x = as.matrix(x)
  y = as.matrix(y)
  K1x = (x - 1 / 2)
  K1y = (y - 1 / 2)
  K2x = (K1x^2 - 1 / 12) / 2
  K2y = (K1y^2 - 1 / 12) / 2
  ax = x %x% matrix(1, 1, nrow(y)) 
  ay = y %x% matrix(1, 1, nrow(x))
  b = abs(ax - t(ay))
  K1 = K1x %x% t(K1y)
  K2 = K2x %x% t(K2y) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
  list(K1 = K1, K2 = K2)
}

dgaussian = function(X, xp, kparam = 1)
{
  gamma = kparam
  # diff_mat = sweep(X, xp, MARGIN = 2)
  # diff_mat = t(t(X) - xp)
  np = dim(X)
  diff_mat = X - matrix(xp, nrow = np[1], ncol = np[2], byrow = TRUE)
  # diff_mat = X - xp[col(X)]
  # tmp = exp(-rowSums(diff_mat^2) / (2 * sigma)) * (diff_mat / sigma)
  tmp = exp(-rowSums(diff_mat^2) * gamma) * (2 * gamma * diff_mat)
  # dgd = drop(crossprod(alpha, tmp))
  return(tmp)
}

ddgaussian = function(X, xp, kparam = 1, comb_set)
{
  gamma = kparam
  np = dim(X)
  diff_mat = X - matrix(xp, nrow = np[1], ncol = np[2], byrow = TRUE)
  # diff_mat2 = diff_mat[, comb_set]
  
  diff_mat2 = sapply(1:ncol(comb_set), function(i) diff_mat[, comb_set[, i][1]] * diff_mat[, comb_set[, i][2]])
  
  # tmp = exp(-rowSums(diff_mat^2) * gamma) * 4 * gamma^2 * (diff_mat[, ind1]) * (diff_mat[, ind2])
  tmp = exp(-rowSums(diff_mat^2) * gamma) * 4 * gamma^2 * diff_mat2
  # dgd = drop(crossprod(alpha, tmp))
  return(tmp)
}



dlinear = function(X, xp, kparam = 1)
{
  return(X)
}


dpoly = function(X, xp, kparam = 1)
{
  # degree = kernel_par$degree
  # scale = kernel_par$scale
  # offset = kernel_par$offset
  degree = kparam
  scale = 1
  offset = 0
  tmp = degree * (scale * drop(X %*% xp) + offset)^{degree - 1} * scale * X
  return(tmp)
}

ddspline = function(x, y)
{
  m1 = (x - (1 / 2)) + (1 / 2) * ((x - (1 / 2))^2 - (1 / 12)) * (y - (1 / 2))
  m2 = ifelse(x > y, -(1 / 6) * (x - y - (1 / 2))^3 + (1 / 24) * (x - y - (1 / 2)), 
              (4 * (y - x - (1 / 2))^3 - (y - x - (1 / 2))) / 24)
  return(m1 - m2)
}


pderiv = function(alpha, x, y, kernel = c("linear", "poly", "gaussian", "spline", "anova_gaussian"), kparam = 1)
{
  n = length(y)
  k = length(unique(y))
  
  dkernel = switch(kernel,
                   linear = dlinear,
                   poly = dpoly,
                   gaussian = dgaussian,
                   spline = dspline,
                   anova_gaussian = dgaussian)
  
  grad_mat = 0
  for (i in 1:n) {
    dK_sq = crossprod(alpha, dkernel(x, x[i, ], kparam))^2
    grad_mat = grad_mat + dK_sq / n
  }
  
  res = colSums(grad_mat)
  return(res)
}



pderiv_so = function(alpha, x, y, kernel = c("linear", "poly", "gaussian"), kparam = 1, active_set = NULL)
{
  n = length(y)
  k = length(unique(y))
  
  ddkernel = switch(kernel,
                    linear = ddlinear,
                    poly = ddpoly,
                    gaussian = ddgaussian)
  
  comb_set = combn(active_set, 2)
  
  grad_mat = 0
  for (i in 1:n) {
    dK_sq = crossprod(alpha, ddkernel(x, x[i, ], kparam, comb_set))^2
    grad_mat = grad_mat + dK_sq / n
  }
  
  res = colSums(grad_mat)
  return(res)
}

data_split = function(y, nfolds, seed = length(y))
{
  # k: the number of classes
  y = as.factor(y)
  n_data = length(y)
  n_class = length(levels(y))
  class_size = table(y)
  classname = names(class_size)
  ran = rep(0, n_data) 
  if ((min(class_size) < nfolds) & (nfolds != n_data))
  {
    warning('The given fold is bigger than the smallest class size. \n Only a fold size smaller than the minimum class size \n or the same as the sample size (LOOCV) is supported.\n')
    return(NULL)
  }
  
  if (min(class_size) >= nfolds) {
    set.seed(seed)
    for (j in 1:n_class) {  
      ran[y == classname[j]] = ceiling(sample(class_size[j]) / (class_size[j] + 1) * nfolds) 
    }
  }
  else if (nfolds == n_data) {
    ran = 1:n_data
  }
  return(ran)
}



strong_heredity = function(main_effect, interaction)
{
  if (sum(interaction) != 0) {
    p = length(main_effect)
    comb_mat = combn(1:p, 2)
    ind_mat = comb_mat[, interaction == 1, drop = FALSE]
    for (i in 1:nrow(ind_mat)) {
      ind = ind_mat[i, ]
      main_effect[ind] = 1
    }
  }
  res = c(main_effect, interaction)
  return(res)
}


# interaction_graph = function(comb, p, min = 3)
# {
#   int_mat = Matrix::Matrix(0, nrow = p, ncol = p)
#   int_mat[t(comb)] = 1
#   int_mat = Matrix::t(int_mat) + int_mat
#   g = graph_from_adjacency_matrix(int_mat, mode = "undirected")
#   cliques_list = max_cliques(g, min = 3)
#   return(cliques_list)
# }

interaction_kernel = function(x, u, kernel, kparam, active_set, interaction_set)
{
  if (!is.matrix(x)) {
    x = as.matrix(x)
  } else {
    x = as.matrix(x)
  }
  u = as.matrix(u)
  dimx = ncol(x)
  
  scaled_kernel = function(x, u, kernel, kparam, active_set, index)
  {
    X1 = matrix(rowSums(x[, active_set, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u))
    U1 = matrix(rowSums(u[, active_set, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u), byrow = TRUE)
    X2 = matrix(rowSums(x[, index, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u))
    U2 = matrix(rowSums(u[, index, drop = FALSE]^2), nrow = nrow(x), ncol = nrow(u), byrow = TRUE)
    K = exp(-kparam * ((X1 + U1) - (X2 + U2)))
    K1 = exp(-kparam * (X1 + U1))
    K2 = exp(-kparam * (X2 + U2))
    K_mat = kernelMatrix(x[, index, drop = FALSE], u[, index, drop = FALSE], kernel = kernel, kparam = kparam)
    res = K * K_mat - K1
    return(list(res = res, K = K, K1 = K1, K2 = K2))
  }
  
  main_effects = vector(mode = "list", dimx)
  
  for (j in 1:length(active_set)) {
    temp_kernel = scaled_kernel(x, u, kernel = kernel, kparam = kparam, active_set = active_set, index = active_set[j])
    main_effects[[active_set[j]]] = temp_kernel$res
  }
  
  interaction_kernel = 0
  if (ncol(interaction_set) != 0) {
    
    for (i in 1:ncol(interaction_set)) {
      ind = interaction_set[, i]
      interaction_kernel = interaction_kernel + ((main_effects[[ind[1]]]) * (main_effects[[ind[2]]])) / temp_kernel$K1
    }
  }
  
  K = temp_kernel$K1 + Reduce("+", main_effects[active_set]) + interaction_kernel
  
  return(K)
}

code_ramsvm = function(y)
{
  n_class = length(unique(y))
  n = length(y)
  yyi = Y_matrix_gen(k = n_class, nobs = n, y = y)
  W = XI_gen(n_class)
  
  y_index = cbind(1:n, y)
  index_mat = matrix(-1, nrow = n, ncol = n_class)
  index_mat[y_index] = 1
  
  Hmatj = list()
  Lmatj = list()
  for (j in 1:(n_class - 1)) {
    Hmatj_temp = NULL
    Lmatj_temp = NULL
    for (i in 1:n_class) {
      temp = diag(n) * W[j, i]
      diag(temp) = diag(temp) * index_mat[, i]
      Hmatj_temp = rbind(Hmatj_temp, temp)
      Lmatj_temp = c(Lmatj_temp, diag(temp))
    }
    Hmatj[[j]] = Hmatj_temp
    Lmatj[[j]] = Lmatj_temp
  }
  return(list(yyi = yyi, W = W, Hmatj = Hmatj, Lmatj = Lmatj, y_index = y_index))
}

find_nonzero = function(Amat)
{
  nr = nrow(Amat)
  nc = ncol(Amat)
  Amat_compact = matrix(0, nr, nc)
  Aind = matrix(0, nr + 1, nc)
  for (j in 1:nc) {
    index = (1:nr)[Amat[, j] != 0]
    number = length(index)
    Amat_compact[1:number, j] = Amat[index, j]
    Aind[1, j] = number
    Aind[2:(number+1), j] = index
  }
  max_number = max(Aind[1, ])
  Amat_compact = Amat_compact[1:max_number, ]
  Aind = Aind[1:(max_number + 1), ]
  return(list(Amat_compact = Amat_compact, Aind = Aind))
}

fixit = function(A, epsilon = .Machine$double.eps, is_diag = FALSE)
{
  if (is_diag) {
    d = diag(A)
    tol = epsilon
    eps = max(tol * max(d), 0)
    d[d < eps] = eps
    Q = diag(d)
  } else {
    eig = eigen(A, symmetric = TRUE)
    tol = epsilon
    eps = max(tol * abs(eig$values[1]), 0)
    eig$values[eig$values < eps] = eps
    Q = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  }
  return(Q)
}

