
galasso<-function(norm,lambda)
{
       if (norm <= lambda)
         res = 0
       else 
         res = norm - lambda
       return(res)
}   

gmcp<-function(a,b,c)
{
     if (abs(a) <= b) return(0)
     else
         {
           if(a > b) return((a-b)/c)     
           else      return((a+b)/c)
         }
       
}   


solve.new <- function(mat)
{
  k = nrow(mat)
  tiny = matrix(0,k,k)
  diag(tiny) = 0.1*rep(1,k)
  p = svd(mat+tiny)
  f = p$d
  f[f <= 0.001] = 0.001
  m = matrix(0,k,k)
  diag(m) = f^(-1)
  inv = p$v%*%m%*%t(p$u)
  return(inv)
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
