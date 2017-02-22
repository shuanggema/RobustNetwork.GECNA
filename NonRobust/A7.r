#  Created by Cen Wu on 06/2014
#  Copyright 2014-2017 Yale University. All rights reserved.


A7 <- function(Y,bet,X,Y.test,X.test,tol,eps) 
{     
          lambda1 =  seq(100, 140, by = 10)/300  ; nlambda1 = length(lambda1);
          lambda3 =  10^(seq(-2,-1,by=1))  ; nlambda3 = length(lambda3);
          n = dim(Y)[1]
          p = dim(bet)[1]
          m = dim(bet)[2]
          gamma = 3; 
          kn = degree = 2; 
          q = kn+degree+1;

          coef.n = array(0,dim=c(q,m,nlambda1,nlambda3))
          coef.p = array(0,dim=c(p,m,nlambda1,nlambda3))
          beta0.p = matrix(0,p,m)   # the parametric components
          beta0.n = matrix(0,q,m)  # the non parametric components 
          design.n = NULL
          res = matrix(0,n,m)  
          cv.test = array(0, dim=c(nlambda1, nlambda3))
          A = N.3(X); 

           for (i in 1:m)
            {
              xx = X
              y = as.matrix(Y[,i])        
              lasso.cv <- cv.glmnet(xx,y,alpha=1,nfolds=5)   
              alpha <- lasso.cv$lambda.min
              lasso.fit1 <- glmnet(xx,y,family="gaussian",alpha=1,nlambda=100,intercept=FALSE)
              coef0 <- as.matrix(as.vector(predict(lasso.fit1, s=alpha, type="coefficients"))[-1])    
              beta0.p[,i] = as.matrix(coef0)
              res[,i] = y - as.matrix(X) %*% as.matrix(beta0.p[,i])
            }
            
            res0 = res


          for (i1 in 1:nlambda1)
             {       
               for (i3 in 1:nlambda3)
                {   
                   lam1 = lambda1[i1]    
                   lam2 = lambda1[i1]*q^(1/2) 
                   lam3 = lambda3[i3]        
                   beta.p = beta0.p
                   res = res0;
                   diff.in = diff.out = 1  
                   step.in = step.out = 0

                   while(diff.out>=10^{-4} & step.out<=30) 
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
                                           if (abs(beta.p[i,k])< gamma*lam1) {
                                                u = 1 -1/gamma + lam3*sum(abs(A[i,((i+1):p)]))  
                                             } else  { 
                                                u = 1 + lam3*sum(abs(A[i,((i+1):p)])) } 
                                             if (abs(beta.p[i,k])< gamma*lam1)
                                                 {
                                                     v = lam1 
                                                 } else  { v = 0 }
                                              
                                               w = t(as.matrix(X[,i]))%*%(res[,k]+as.matrix(X[,i])*beta.p[i,k])/n + lam3*sum(A[i,((i+1):p)]*beta.p[((i+1):p),k])                                
                                               beta.p[i,k] = gmcp(w,v,u); 
                                               res[,k] = res[,k]- as.matrix(X[,i])*(beta.p[i,k]- old_beta.p)
                                           }                

                                      if (i==p) # update the parametric coefficients
                                        {
                                             old_beta.p = beta.p[i,k]
                                           if (abs(beta.p[i,k])< gamma*lam1) {
                                                u = 1 -1/gamma  
                                             } else  { 
                                                u = 1  } 
                                             if (abs(beta.p[i,k])< gamma*lam1)
                                                 {
                                                     v = lam1 
                                                 } else  { v = 0 }
                                               w = t(as.matrix(X[,i]))%*%(res[,k]+as.matrix(X[,i])*beta.p[i,k])/n                                 
                                               beta.p[i,k] = gmcp(w,v,u); 
                                               res[,k] = res[,k]- as.matrix(X[,i])*(beta.p[i,k]- old_beta.p)
                                           }         
                                  }                               
                               } 
                          diff.out=mean((beta.p-beta.p.out)^2)
                          step.out=step.out+1
                       } #  end of while loop 

                   ###  cross validation on an independent dataset 
                   res.test = matrix(0,n,m)    # the residual matrix
                   res.test = Y.test - as.matrix(X.test) %*% as.matrix(beta.p)     
                   cv.test[i1,i3] = sum(res.test^2/n)            
                   coef.p[,,i1,i3] = beta.p           
             } 
           } 

          cv = cv.test
          i1.cv = which(cv == min(cv), arr.ind = TRUE)[1,1]
          i3.cv = which(cv == min(cv), arr.ind = TRUE)[1,2]
          beta2.p = coef.p[,,i1.cv,i3.cv]  
          beta3.p = beta2.p;
          diag(beta3.p) = 0;  
          bet.n = seq(1:200)
          bet.p = c(seq(6,10),seq(56,60),seq(101,200))
          bet.n[bet.p] = 0
          TP1 = sum((which(bet!=0) %in% which(beta3.p!=0))*1) 
          TP2 = sum((which(bet.n!=0) %in% which(diag(beta2.p)!=0))*1) 
          TP = TP1+TP2

          FP1 = sum((which(beta3.p!=0) %in% which(bet==0))*1)
          FP2 = sum( (which(diag(beta2.p)!=0) %in% which(bet.n==0))*1)
          FP = FP1 + FP2 


        
tun = list(beta.p=beta2.p,TP=TP, FP = FP, TP1=TP1, FP1=FP1,TP2=TP2, FP2=FP2)
return(tun)
}  




