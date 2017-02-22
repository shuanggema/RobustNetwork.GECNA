#  Created by Cen Wu on 07/2014
#  Updated on 09/2014
#  Copyright 2014-2017 Yale University. All rights reserved.



A3 <- function(Y,bet,X,Y.test,X.test,tol,eps,m,n,p) 
{     
      
          lambda1 =   seq(53, 58, by = 1)/300   ; nlambda1 = length(lambda1);
          lambda3 =  10^(seq(-2,-1,by=1))  ; nlambda3 = length(lambda3)

          n = dim(Y)[1]
          p = dim(bet)[1]
          m = dim(bet)[2]
          nt = 100*n; 
          kn = degree = 2; 
          q = kn+degree+1;
          gamma = 3; 
          tau = 0.5;

          coef.p = array(0,dim=c(p,m,nlambda1,nlambda3))
          beta0.p = matrix(0,p,m)                             # the parametric components
          design.n = NULL
          cv.test = array(0, dim=c(nlambda1, nlambda3))
          A = N.3(X); 
          res = matrix(0,n,m)    # the residual matrix

            for (i in 1:m)
            {
              xx = X
              y = as.matrix(Y[,i])      
              lasso.cv <- cv.glmnet(xx,y,alpha=1,nfolds=5)   
              alpha <- lasso.cv$lambda.min/3
              lasso.fit1 <- glmnet(xx,y,family="gaussian",alpha=1,nlambda=100,intercept=FALSE)
              coef0 <- as.matrix(as.vector(predict(lasso.fit1, s=alpha, type="coefficients"))[-1])  # sum((y[,1]-xx%*%coef0)^2)/n  
              beta0.p[,i] = as.matrix(coef0)
              res[,i] = y - as.matrix(X) %*% as.matrix(beta0.p[,i])
            }
            
            res0 = res     

          for (i1 in 1:nlambda1)
             {       
               for (i3 in 1:nlambda3)
                {   
                   lam1 = lambda1[i1]   
                   lam3 = lambda3[i3]    
                   beta.p = beta0.p
                   res = res0
                   diff.in = diff.out = 1  
                   step.in = step.out = 0

                   while(diff.out>=10^{-4} & step.out<=25) 
                      {
                            beta.p.out = beta.p
                         for(k in 1:m) 
                              {
                                   y = as.matrix(Y[,k]); 
                                   for(i in 1:p) 
                                     {
                                        if (i<p) # update the parametric coefficients
                                           {
                                               old_beta.p = beta.p[i,k]
                                               u.m = matrix(0,n,1)
                                               w.m = matrix(0,n,1)                                             
                                             
                                               u.m = (y - X %*% as.matrix(beta.p[,k]) + X[,i]*beta.p[i,k])/X[,i]
                                               for (s in 1:n)
                                                   {
                                                      w.m[s] = abs(X[s,i]*(tau - ((u.m[s]*X[s,i])< 0)))/n
                                                   }

                                               u.m = matrix(c(u.m,rep(0,1+p-i)),n+1+p-i,1)
                                               w.m = matrix(c(w.m,rep(0,1+p-i)),n+1+p-i,1)
                                               w.m[n+1] = dmcp(abs(beta.p[i,k]),gamma,lam1)

                                               for (s in (n+2):(n+1+p-i))
                                                   {
                                                      st = s-n-1;
                                                      u.m[s] = abs(beta.p[(i+st),k])
                                                      w.m[s] = lam3*abs(A[i,(i+st)])
                                                   }

                                               beta.p[i,k] = weighted.median(u.m,w.m) 
                                          
                                            }               

                                        if (i==p) 
                                            {

                                               old_beta.p = beta.p[i,k]
                                               u.m = matrix(0,n,1)
                                               w.m = matrix(0,n,1)                                             
                                               u.m = (y - X %*% as.matrix(beta.p[,k]) + X[,i]*beta.p[i,k])/X[,i]
                                               for (s in 1:n)
                                                   {
                                                      w.m[s] = abs(X[s,i]*(tau - ((u.m[s]*X[s,i])< 0)))/n
                                                   }
                                               u.m = matrix(c(u.m,0),n+1,1)
                                               w.m = matrix(c(w.m,dmcp(abs(beta.p[i,k]),gamma,lam1)),n+1,1)                                              
                                               beta.p[i,k] = weighted.median(u.m,w.m) 
                                            }   

                                  }                                   
                               } 
                          diff.out=mean((beta.p-beta.p.out)^2); 
                          step.out=step.out+1
                       } #  end of while loop 
              
                   ####  validation on an independent dataset 
 
                   res.test = matrix(0,nt,m)    
                   res.test = Y.test - as.matrix(X.test) %*% as.matrix(beta.p) 
                   for (i.test in 1:m) { cv.test[i1,i3] = cv.test[i1,i3] + sum(res.test[,i.test]*(tau-(res.test[,i.test]<0)))}
                   coef.p[,,i1,i3] = beta.p         
             } 
           } 

          cv = cv.test/nt        
          i1.cv = which(cv == min(cv), arr.ind = TRUE)[1,1]
          i3.cv = which(cv == min(cv), arr.ind = TRUE)[1,2]
          flag.lambda1 = lambda1[i1.cv]          
          flag.lambda3 = lambda3[i3.cv]

          beta2.p = coef.p[,,i1.cv,i3.cv]  
          bet.n = seq(1:200)
          bet.p = c(seq(6,10),seq(56,60),seq(101,200))
          bet.n[bet.p] = 0

          TP1 = sum((which(bet!=0) %in% which(beta2.p!=0))*1) 
          TP2 = sum((which(bet.n!=0) %in% which(diag(beta2.p)!=0))*1) 
          TP = sum((which(bet!=0) %in% which(beta2.p!=0))*1) + sum((which(bet.n!=0) %in% which(diag(beta2.p)!=0))*1) 


          FP1 = sum((which(beta2.p!=0) %in% which(bet==0))*1)
          FP2 = sum( (which(diag(beta2.p)!=0) %in% which(bet.n==0))*1)
          FP = sum((which(beta2.p!=0) %in% which(bet==0))*1) + sum( (which(diag(beta2.p)!=0)  %in% which(bet.n==0))*1)         
               
tun = list(beta.p=beta.p,TP=TP, FP=FP,TP1=TP1, FP1 = FP1,TP2=TP2, FP2 = FP2)
return(tun)
}  




