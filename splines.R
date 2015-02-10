#load ('to.splines.RData')

require (plyr)
require (geomorph)
require (expm)
require (numDeriv)

Norm <- function (x) return (sqrt (sum (x * x)))

interpol.TPS <- function (x, y)
      {
        distance <- Norm (x - y)
        if (length (x) == 2)
          {
            if (distance < .Machine $ double.eps)
              return (0)
            else
              return (- log (distance ^ 2) * (distance ^ 2))
          }
        if (length (x) == 3)
          return (distance)
      }

interpol.L1.TPS <- function (y, X)
  return (aaply (X, 1, interpol.TPS, y = y))

interpol.L2.TPS <- function (X, Y)
  return (aaply (Y, 1, interpol.L1.TPS, X = X))

theta.p <- function (p)
    return (b + (A %*% p) + t (W) %*% interpol.L1.TPS (p, Q))

## interpol.EBS <- function (x, y, nu = 0.25)
##       {
##         p <- x - y
##         distance <- Norm (p)
##         alpha <- 12 * (1 - nu) - 1
##         Gp <- alpha * distance ^ 2
##         Gp <- Gp * diag (3)
##         Gp <- Gp - (3 * outer (p, p))
##         Gp <- Gp * distance
##         return (Gp)
##       }

## interpol.L1.EBS <- function (y, X, ...)
##   return (aaply (X, 1, interpol.EBS, y = y, ...))

## interpol.L2.EBS <- function (X, Y, ...)
##   return (aaply (Y, 1, interpol.L1.EBS, X = X, ...))

TPS <- function (target.shape, reference.shape)
  {
    ## calculates thin plate splines deformations between reference and target shapes
    theta.p <- function (p)
      return (b + (A %*% p) + t (W) %*% interpol.L1.TPS (p, Q))
    Q <- reference.shape
    P <- target.shape
    k <- nrow (Q)
    m <- ncol (Q)
    R <- interpol.L2.TPS (Q, Q)
    # print (R)
    one.k <- array (1, c (k, 1))
    zero.11 <- array (0, c (1, 1))
    zero.m1 <- array (0, c (m, 1))
    zero.mm <- array (0, c (m, m))
    L <- rbind (R, t (one.k), t (Q))
    L <- cbind (L, rbind (one.k, zero.11, zero.m1))
    L <- cbind (L, rbind (Q, t (zero.m1), zero.mm))
    P.zero <- rbind (P, t (zero.m1), zero.mm)
    spline <- solve (L, P.zero)
    W <- spline [1:k, ]
        b <- t (t (spline [1+k, ]))
    A <- t (spline [(2:(m+1))+k, ])
    return (list ('W' = W, 'A' = A,
                  'b' = b, 'Q' = Q, 'P' = P, 'theta.p' = theta.p))
  }

## EBS <- function (target.shape, reference.shape, nu = 0.25)
##     {
##         ## calculates elastic body splines between target and reference shapes (Davis et. al. 1995)
##         ## note: notation switched to match TPS (Q: reference, P: target)
##         Q <- reference.shape
##         P <- target.shape
##         N <- nrow (Q)
##         Y <- t (P - Q)
##         dim (Y) = NULL
##         Y <- c (Y, rep (0, times = 12))
##         K <- interpol.L2.EBS (Q, Q, nu = nu)
##         K <- aperm (K, c (3, 1, 4, 2), resize = TRUE)
##         dim (K) <- c (3 * N, 3 * N)
##         Qmat <- cbind (Q %x% diag (3), rep (1, times = N) %x% diag (3))
##         O <- array (0, c (12, 12))
##         L <- cbind (K, Qmat)
##         L <- rbind (L, cbind (t (Qmat), O))
##         W <- solve (L, Y)
##         C <- t (array (W [1:(3*N)], c(3, N)))
##         A <- t (array (W [(3*N)+(1:9)], c (3, 3)))
##         b <- W [(3*N)+(10:12)]
##         theta.p <- function (p)
##             {
##                 Gp.s <- interpol.L1.EBS (p, Q, nu = nu)
##                 sum.ux <- array (0, dim (C))
##                 for (i in 1:dim (Gp.s) [1])
##                     sum.ux [i, ] <- Gp.s [i, , ] %*% as.vector (C [i, ])
##                 d.x <- colSums (sum.ux) + as.vector (A %*% p) + as.vector (b)
##                 return (d.x)
##             }
##         return (list ('C' = C, 'A' = A, 'b' = b,
##                       'Q' = Q, 'theta.p' = theta.p))
##     }

Rotate2MidlineMatrix <- function (X, midline)
  {
    ## returns the rotation matrix that aligns a specimen saggital line
    ## to plane y = 0 (2D) or z = 0 (3D)
    ncl <- ncol (X) 
    Xm <- na.omit (X [midline, ])
    Mm <- matrix (apply (Xm, 2, mean), byrow = TRUE, nr = nrow (X), nc = ncl)
    Xc <- X - Mm 
    W <- na.omit (Xc [midline, ])
    RM <-svd (var (W))$v
    return (RM)
  }

JacobianArray <- function (spline, tesselation, ...)
    {
        ## calculates jacobians for a given interpolation in a set of points
        ## determined from tesselation (as centroids of each tetrahedron defined, for now...)
        with (spline,
              {
                  Q.tetra <- Q [tesselation, ]
                  dim (Q.tetra) <- c (dim (tesselation), ncol (Q))
                  Q.centroids <- apply (Q.tetra, 1, colMeans)
                  jacobs <- apply (Q.centroids, 2, jacobian, func = theta.p, ...)
                  dim (jacobs) <- c (ncol (Q), ncol (Q), ncol (Q.centroids))
                  return (jacobs)
              })
    }

Center2MeanJacobian <- function (jacobArray, max.steps = 100)
    {
        ### calculates mean jacobian matrix for a set of jacobian matrices
        ### describing a local aspect of shape deformation for a given set of volumes, 
        ### returning log determinants of deviations from mean jacobian (Woods, 2003).
        logm.single <- function (Ai, inv.Mk) return (logm (Ai %*% inv.Mk))
        A <- jacobArray
        N <- dim (A) [3]
        Mk <- diag (nrow (A))
        i <- 1
        repeat
            {
              inv.Mk <- solve (Mk)
              centered.now <- apply (A, 3, logm.single, inv.Mk = inv.Mk)
              o <- array (- rowMeans (centered.now), c(nrow (A), nrow (A)))
              if (all (abs (o) < .Machine$double.eps))
                {
                    break
                }
              Mk <- expm (- o) %*% Mk
              i <- i + 1
              if (i == max.steps)
                {
                  # print ('Did not converge.')
                  stop ('Convergence has not been achieved in number of steps.')
                }
            }
        centered.now <- array (centered.now, c(nrow (A), nrow (A), N))
        log.det <- apply (centered.now, 3, function (x) return (sum (diag (x))))
        return (log.det)
    }


wrapMarquez <- function (OTU, tesselation, which = 'sym')
    {
      with (OTU,
            {
                if (which == 'sym')
                    {
                        gpa <- gpagen (sym, ShowPlot = FALSE) $ coords
                        mshape <- mshape (gpa)
                    }
                else if (which == 'ass')
                    {
                        gpa <- gpagen (ass, ShowPlot = FALSE) $ coords
                        mshape <- mshape (gpa)
                    }
                dimnames (mshape) <- dimnames (sym.mean)
                tps <- alply (gpa, 3, TPS,
                              reference.shape = mshape, .parallel = TRUE)
                print ('tps done')
                jacobs <- laply (tps, JacobianArray, tesselation = tesselation,
                                 .parallel = TRUE)
                jacobs <- aperm (jacobs, c(2, 3, 1, 4), resize = TRUE)
                print ('jacobs done')
                local <- aaply (jacobs, 4, Center2MeanJacobian, .parallel = TRUE)
                local <- t (local)
                return (list ('tps' = tps,
                              'jacobians' = jacobs,
                              'local' = local,
                              'reference' = mshape,
                              'cs' = cs,
                              'info' = info))
            })
  }
