make_D = function(dimension, order)
{
  D = diag(rep(-1, dimension))
  for(j in 1:(dimension - 1))
    D[j, (j + 1)] = 1
  if(order > 2)
    for(k in 1:(order - 2))
      D = D %*% D
  return(D)
}

c.basic.interval = function(t, knots, order)
{
   dimension = length(knots) - order
   (t - knots[2:dimension]) / 
      (knots[(2 + order - 1) : (dimension + order - 1)] - knots[2:dimension])
}

# t : point-wise
spline.fit = function(t, knots, order, xi)
{
   dimension = nrow(xi)
   if(order > 1)
   {
      # c_j vector
      c = c.basic.interval(t, knots, order)
      xi_update = c * xi[2:dimension, ] + (1 - c) * xi[1:(dimension - 1), ]
      return(spline.fit(t, knots = knots[-c(1, length(knots))], order = order - 1, xi = xi_update))
   }
   else
   {
      # constant fit
      return(as.vector(bsplines(t, knots, order) %*% xi))
   }
}

########################################################################################################
# spline box R-version 1.6

# b-spline basis for one predcitor
bspline = function(x, knots, order = 1, derivative = 0, j = 1)
{
   bspline_vector = rep(0, length(x))
   support = knots[j] <= x & x < knots[j + order]
   if (derivative == 0)
      bspline_vector[support] = bspline0(x[support], knots, order, j)
   else
      bspline_vector[support] = bspline_derivative(x[support], knots, order,
                                                   derivative, j)
   return(bspline_vector)
}

# collection(matrix) of b-spline basis(each column)
bsplines = function(x, knots, order = 2, derivative = 0)
{
   dimension = length(knots) - order
   bspline_matrix = matrix(0, length(x), dimension)
   for (j in 1 : dimension)
      bspline_matrix[, j] = bspline(x, knots, order, derivative, j)
   return(bspline_matrix)
}

# b-spline basis
bspline0 = function(x, knots, order, j)
{
   if (order > 1)
   {
      if (knots[j + order - 1] > knots[j])
         a = (x - knots[j]) / (knots[j + order - 1] - knots[j])
      else
         a = 0
      if (knots[j + order] > knots[j + 1])
         b = (knots[j + order] - x) / (knots[j + order] - knots[j + 1])
      else
         b = 0
      return(a * bspline0(x, knots, order - 1, j) +
                b * bspline0(x, knots, order - 1, j + 1))
   }
   else
   {
      bspline = rep(0, length(x))
      bspline[knots[j] <= x & x < knots[j + 1]] = 1
      return(bspline)
   }
}

# derivative of b-spline basis
bspline_derivative = function(x, knots, order, derivative, j)
{
   if (derivative > 0)
   {
      if (knots[j + order - 1] > knots[j])
         a = (order - 1) / (knots[j + order - 1] - knots[j])
      else
         a = 0
      if (knots[j + order] > knots[j + 1])
         b = (order - 1) / (knots[j + order] - knots[j + 1])
      else
         b = 0
      return(a * bspline_derivative(x, knots, order - 1, derivative - 1, j) -
                b * bspline_derivative(x, knots, order - 1, derivative - 1, j + 1))
   }
   else
      return(bspline0(x, knots, order, j))
}

# compute jump matrix for b-spline
bspline_jump = function(knots, order, transpose = FALSE)
{
   dimension = length(knots) - order
   number_interior_knots = dimension - order
   midpoint_between_knots = (knots[order : (length(knots) - order)]
                             + knots[(order + 1) : (length(knots) - order + 1)]) / 2
   # transpose of theoritical jump matrix
   jump_matrix = matrix(nrow = dimension, ncol = number_interior_knots)
   for (j in 1 : dimension)
   {
      derivatives = bspline(midpoint_between_knots, knots, order, order - 1, j)
      for (l in 1 : number_interior_knots)
         jump_matrix[j, l] = derivatives[l + 1] - derivatives[l]
   }
   # re-transpose the transposed jump matirx
   if (transpose == FALSE)
      jump_matrix = t(jump_matrix)
   return(jump_matrix)
}


bspline_jump_xi = function(knots, order, transpose = FALSE)
{
  dimension = length(knots) - order
  number_interior_knots = dimension - order
  midpoint_between_knots = (knots[order : (length(knots) - order)]
                            + knots[(order + 1) : (length(knots) - order + 1)]) / 2
  # transpose of theoritical jump matrix
  jump_matrix = matrix(nrow = dimension, ncol = number_interior_knots)
  for (j in 1 : dimension)
  {
    derivatives = bspline(midpoint_between_knots, knots, order, order - 1, j)
    for (l in 1 : number_interior_knots)
      jump_matrix[j, l] = derivatives[l + 1] - derivatives[l]
  }
  # re-transpose the transposed jump matirx
  if (transpose == FALSE)
    jump_matrix = t(jump_matrix)
  return(jump_matrix)
}


# setting knots by quantile of data
knots_quantile = function(x, dimension, order = 2, type = "bs")
{
   if (type == "bs")
   {
      dimension = max(dimension, order)
      number_interior_knots = dimension - order
      if (number_interior_knots > 0)
         probs = (1 : number_interior_knots) / (number_interior_knots + 1)
      else
         probs = NULL
   }
   else
   {
      dimension = max(dimension, 2)
      probs = seq(0, 1, length = dimension)
      order = 4
   }
   interior_knots = quantile(x, probs, type = 1)
   knots = add_boundary_knots(x, interior_knots, order)
   return(knots)
}

# adding boundary knots to interior knots
add_boundary_knots = function(x, interior_knots, order = 1, tiny = 1e-5)
{
   #knots = c(to = seq(min(x) - tiny, by = tiny, length = order), interior_knots, seq(from = max(x) + tiny, by = tiny, length = order))
   knots = c(rep(min(x) - tiny, order), interior_knots, rep(max(x) + tiny, order))
   return(knots)
}
