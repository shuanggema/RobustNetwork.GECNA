
dyn.load("quicksort.so")

weighted.median <- function(x1,x2) 
{  

          x3 = cbind(x1,x2)
          n = length(x1)
          w = 1:n
          fit<-.C("quicksort",as.double(x1), as.double(w), as.integer(n))
          x = x3[fit[[2]],]
          weight = x[,2]
          cum.weight = cumsum(x[,2])
          norm.weight = cumsum(x[,2])/sum(x[,2])
          p.1 = min(which(norm.weight>=0.5))
          w.m = x[p.1]
          return(w.m)
}



# adjacency matrix

N.3 = function(X)
{
   n = nrow(X)
   p = ncol(X)
   r = cor(X); r[which(r==1)] = 1 - 0.01
   z = 0.5*log((1+r[upper.tri(r)])/(1-r[upper.tri(r)]))
   c0 = mean(sqrt(n-3)*z) + 2*sd(sqrt(n-3)*z)
   cutoff = (exp(2*c0/sqrt(n-3))-1)/(exp(2*c0/sqrt(n-3))+1)
   r = cor(X)
   A = (r)^5*(abs(r)>cutoff)
   diag(A) = 0
   A
}


dmcp <- function(x, gamma , lambda) {
        lambda * (1 - abs(x)/(gamma*lambda))*(abs(x) <= gamma*lambda)
    }


gALAL <- function (xx, y, tau, beta=0.9995, eps, wpp) 
{
    n <- length(y)
    p <- ncol(xx)
    if (n != nrow(xx)) 
        stop("xx and y don't match n")
    lambda =   wpp * n
    if (length(lambda) != p) 
        stop(paste("lambda must be either of length ", p, " or length one"))
    if (any(lambda < 0)) 
        stop("negative lambdas disallowed")
    if (tau < eps || tau > 1 - eps) 
        stop("No parametric Frisch-Newton method.  Set tau in (0,1)")
    R <- matrix(0,length(lambda),length(lambda))
    diag(R) <- lambda
    index <- which(lambda != 0)
    len = length(index)
  
    if (len>0)
       {
    R <- R[index, ]
    if (len==1){R=as.matrix(t(R))}
    r <- rep(0, len)
    X <- rbind(xx, R)
    Y <- c(y, r)
    N <- length(Y)
    rhs <- (1 - tau) * apply(xx, 2, sum) + 0.5 * apply(R, 2, sum)
    d <- rep(1, N)
    u <- rep(1, N)
    wn <- rep(0, 10 * N)
    wn[1:N] <- c(rep(1 - tau, n), rep(0.5, nrow(R)))
    z <- .Fortran("rqfnb", as.integer(N), as.integer(p), a = as.double(t(as.matrix(X))), 
        c = as.double(-Y), rhs = as.double(rhs), d = as.double(d), 
        as.double(u), beta = as.double(beta), eps = as.double(eps), 
        wn = as.double(wn), wp = double((p + 3) * p), aa = double(p * p), 
        it.count = integer(3), info = integer(1), PACKAGE = "quantreg")
    if (z$info != 0) 
        stop(paste("Error info = ", z$info, "in stepy2: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(xx)[[2]]
    residuals <- y - xx %*% coefficients
      }
    list(coefficients = coefficients, residuals = residuals)
}

