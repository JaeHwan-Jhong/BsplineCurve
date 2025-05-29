# bspline curve version 1

# bspline curves

alpha = function(s, v, w)
{
   theta = dist(v, w)
   R(1 - s, theta) %p% v + R(s, theta) %p% w
}

ebc = function(t, xi, knots, order)
{
   n = length(t)
   m = ncol(xi)
   gamma = matrix(0, n, m)
   dimension = length(knots) - order
   for (l in order : dimension)
   {
      i = knots[l] <= t & t < knots[l + 1]
      gamma[i, ] = ebc.one(t[i], xi, knots, order, l)
   }
   gamma
}

ebc.one = function(t, xi, knots, k, l, r = 1)
{
   if (r < k)
   {
      c = (t - knots[l]) / (knots[l + r] - knots[l])
      c * ebc.one(t, xi, knots, k, l, r + 1) + 
         (1 - c) * ebc.one(t, xi, knots, k, l - 1, r + 1)
   }
   else
      Rep(xi[l, ], length(t))
}

ebc.derivative = function(t, xi, knots, order, m)
{
   n = length(t)
   dgamma = matrix(0, n, 2)
   dimension = length(knots) - order
   for (l in order : dimension)
   {
      i = knots[l] <= t & t < knots[l + 1]
      if(sum(i) > 0)
        dgamma[i, ] = ebc.derivative.one(t[i], xi, knots, order, m, l)
   }
   dgamma
}

ebc.derivative.one = function(t, xi, knots, order, m, l)
{
   if (m >= order)
      return(Rep(c(0, 0), length(t)))
   dm.xi = dxi(xi, knots, order, m, l)
   ebc.one(t, dm.xi, knots, order - m, l)
}

dxi = function(xi, knots, order, m, l)
{
   if (m >= order)
      stop("m must less than order of splines")
   if (m > 0)
   {
      xi = dxi(xi, knots, order, m - 1, l)
      dm.xi = diff(xi[(l - order + m) : l, ]) / (diff(knots[(l - order + 1 + m) : (l + order - m)], lag = order - m) / (order - m))
      xi[(l - order + 1 + m) : l, ] = dm.xi
      xi
   }
   else
      xi
}












ibc = function(t, xi, knots, order)
{
   n = length(t)
   gamma = matrix(0, n, 3)
   dimension = length(knots) - order
   for (l in order : dimension)
   {
      i = knots[l] <= t & t < knots[l + 1]
      gamma[i, ] = ibc.one(t[i], xi, knots, order, l)
   }
   gamma
}

ibc.one = function(t, xi, knots, k, l, r = 1)
{
   if (r < k)
   {
      #s = (t - knots[l]) / (knots[l + r] - knots[l])
      xi.lm1 = ibc.one(t, xi, knots, k, l - 1, r + 1)
      xi.l = ibc.one(t, xi, knots, k, l, r + 1) 
      s = (t - knots[l]) / (knots[l + r] - knots[l])
      alpha(s, xi.lm1, xi.l)
   }
   else
      Rep(xi[l, ], length(t))
}

knots.quantile = function(t, m, order, tiny = 1e-10)
{
   m = max(m, order)
   n = m - order
   if (n > 0)
      probs = (1 : n) / (n + 1)
   else
      probs = NULL
   interior = as.vector(quantile(t, probs, type = 1))
   c(rep(min(t) - tiny, order), interior, rep(max(t) + tiny, order))
}

# spherical operations

c.to.s = function(xyz)
{
   # convert cartesian coordinates to spherical coordinates
   xyz = row.matrix(xyz)
   theta = Acos(xyz[, 3])
   phi = Atan(xyz[, 2], xyz[, 1])
   theta.phi = cbind(theta, phi)
   row.names(theta.phi) = NULL
   theta.phi
}

dist = function(v, w)
{
   # spherical distance of two points v and w on the sphere
   v.dot.w = dot(v, w)
   v.dot.w = pmin(pmax(-1, dot(v, w)), +1)
   acos(v.dot.w)
}

Exp = function(x, v)
{
   # the exponential map on the unit sphere given a base point x and a vector v
   u = norm2(v)
   cos(u) * x + v / A(u)
}

Log = function(x, y)
{
   theta = dist(x, y)
   A(theta) * proj(x, y)
}

proj = function(x, v, normalize = FALSE)
{
   xv = row.match(x, v)
   x = xv[[1]]
   v = xv[[2]]
   v - dot(x, v) * x
}

proj.n = function(x, v)
{
   normalize(proj(x, v))
}

reflection = function(v, w, tau)
{
   theta = dist(v, w)
   if (theta > pi)
      stop("reflection can not be done")
   a = ((tau[3] - tau[2]) / (tau[2] - tau[1])) * theta
   (cos(a) + cos(theta) * sin(a) / sin(theta)) * w - (sin(a) / sin(theta)) * v
}

s.to.c = function(theta.phi)
{
   # convert spherical coordinates to cartesian coordinates
   theta.phi = row.matrix(theta.phi, 2)
   theta = theta.phi[, 1]
   phi = theta.phi[, 2]
   x = sin(theta) * cos(phi)
   y = sin(theta) * sin(phi)
   z = cos(theta)
   xyz = cbind(x, y, z)
   row.names(xyz) = NULL
   xyz
}

# array operations

backward = function(from, to)
{
   if (from >= to)
      from : to
   else
      NULL
}


cross = function(v, w)
{
   v = row.matrix(v)
   w = row.matrix(w)
   cbind(v[, 2] * w[, 3] - v[, 3] * w[, 2],
         v[, 3] * w[, 1] - v[, 1] * w[, 3],
         v[, 1] * w[, 2] - v[, 2] * w[, 1])
}

cross.n = function(v, w)
{
   normalize(cross(v, w))
}

dot = function(v, w) 
{
   vw = row.match(v, w)
   v = vw[[1]]
   w = vw[[2]]
   as.vector(rowSums(v * w))
}

eye = function(m)
{
   diag(1, m, m)
}

forward = function(from, to)
{
   if (from <= to)
      from : to
   else
      NULL
}

norm2 = function(a)
{
   sqrt(dot(a, a))
}

normalize = function(a) 
{
   a / norm2(a)
}

plus = function(a)
{
   pmax(a, 0)
}

prox.norm2 = function(a, lambda = 0)
{
   if (lambda > 0)
      plus(1 - lambda / norm2(a)) * a
   else
      a
}

Rep = function(v, m)
{
   # repeat a vector m times in row storage mode
   t(replicate(m, as.vector(v)))
}

row.match = function(v, w)
{
   # match the nrows of v and w
   v = row.matrix(v)
   w = row.matrix(w)
   if (nrow(w) > nrow(v))
      v = Rep(v, nrow(w))
   else if (nrow(v) > nrow(w))
      w = Rep(w, nrow(v))
   vw = list()
   vw[[1]] = v
   vw[[2]] = w
   vw
}

row.matrix = function(a, m = 3)
{
   matrix(a, ncol = m)
}

Sign = function(z)
{
   z.sign = sign(z)
   z.sign[z.sign == 0] = 1
   z.sign
}

# new operator

"%ni%" = function(k, ind) 
{
   !(k %in% ind)
}

"%p%" = function(x, v) 
{
   n = length(x)
   v = row.matrix(v)
   if (nrow(v) == 1)
      v = Rep(v, n)
   else if (nrow(v) != n)
      stop("x and v can not be multiplied")
   x * v
}

# special functions

A = function(z, tiny = 1e-5)
{
   x = truncate.below(z, tiny)
   x / sin(x)
}

B = function(z)   
{
   # B(z) = z cos(z) / sin(z) = A(z) cos(z)
   A(z) * cos(z)
}

C = function(z, tiny = 1e-5)
{
   x = truncate.below(z, tiny)
   (sin(x) - x * cos(x)) / sin(x)^3
}

Q = function(s, z, tiny = 1e-5)
{
   x = truncate.below(z, tiny)
   (sin(s * x) * cos(x) - s * cos(s * x) * sin(x)) / sin(x)^3
}

R = function(s, z, tiny = 1e-5)
{
   x = truncate.below(z, tiny)
   sin(s * x) / sin(x)
}

truncate.below = function(z, tiny = 1e-100)
{
   Sign(z) * pmax(abs(z), tiny)
}
