# generate simulated data


Block = function(p, d, rho.out, rho.in)
{
     sigma = matrix(rho.out, ncol=p, nrow=p)
     for (i in 1:(p/d))
     {
         sigma[((i-1)*d+1) : (i*d), ((i-1)*d+1) : (i*d)] = rho.in
     }
     diag(sigma) = 1
     sigma
}


AR = function(p, rho)
{
   sigma = outer(1:p, 1:p, FUN = function(x, y) rho^(abs(x - y)))
   sigma
}


Blar = function(p, d, rho.out, rho.in)
{
    sigma = matrix(rho.out, ncol=p, nrow=p)
    for (i in 1:(p/d))
    {
         sigma[((i-1)*d+1) : (i*d), ((i-1)*d+1) : (i*d)] = AR(d, rho.in)
    }
    diag(sigma) = 1
    sigma
}

generator = function(m, n, p, bet, Vstr, Vrho)
{
    require(mvtnorm)
    if (Vstr == "Block") { 
       V= Block(p, 5, 0, Vrho)
    } else if (Vstr == "AR") {
       V = AR(p, Vrho)
    } else if (Vstr == "Blar") {
       V = Blar(p, 5, 0, Vrho)
    }


      U = runif(n)
      X = matrix(0,n,p)
    for(i in 1:n){
        if(U[i]<.6){
            X[i,] = rmvnorm(1, mean=rep(0, p), sigma = V)
         }else{
            X[i,] = rmvnorm(1, mean=rep(0, p), sigma = 0.7*V)
        }
 
      }

    # normalize X
    Center <- colMeans(X)
    X.c <- sweep(X, 2, Center)
    Scale <- sqrt(apply(X.c,2,crossprod)/n)
    XX <- sweep(X.c, 2, Scale,"/")
    X <- XX
     bet.n = seq(1:100)
     bet.p = c(seq(6,10),seq(56,60))
     bet.n[bet.p] = 0
     bet.n = bet.n[ bet.n != 0 ] 

    Z <- matrix(0,n,p)    
    Z[,bet.n] <- X[,bet.n]

     U = runif(n)
     Error = matrix(0,n,p)
    for(i in 1:n){
        if(U[i]< 0.85){
            Error[i,] = rnorm(m)*0.5
        }else{   Error[i,] <- rcauchy(m)*0.05 }
 
      }

    Y = X %*% bet + Error + apply(Z, 2, function(x) 2*sin(x*pi/2))
    list(X=X, Y=Y)
}


