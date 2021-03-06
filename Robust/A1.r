#  Created by Cen Wu on 06/2014
#  updated on  08/2014
#  Copyright 2014-2017 Yale University. All rights reserved.



A1 <- function(Y,bet,X,Y.test,X.test,tol,eps) 
{     


          lambda1 =   seq(58, 62, by = 1)/600   ; nlambda1 = length(lambda1);
          lambda3 =  10^(seq(-2,-1,by=1))  ; nlambda3 = length(lambda3);
    
          n = dim(Y)[1]
          p = dim(bet)[1]
          m = dim(bet)[2]
          nt=100*n;
          gamma = 3; tau = 0.5;
          kn = degree = 2; 
          q = kn+degree+1;

          coef.p = array(0,dim=c(p,m,nlambda1,nlambda3))
          coef.n = array(0,dim=c(q,m,nlambda1,nlambda3))
          beta0.p = matrix(0,p,m)                           # the parametric components
          beta0.n = matrix(0,q,m)                           # the non parametric components 
          design.n = design.n.test = NULL
          res = matrix(0,n,m)   
          cv.test = array(0, dim=c(nlambda1, nlambda3))
          A = N.3(X); 

            for (i in 1:m)
            {
              u = as.matrix(X.test[,i])
              u[which(u>max(X[,i]))]= max(X[,i])
              u[which(u<min(X[,i]))]= min(X[,i])
              u.k = seq(0, 1, length=kn+2)[-c(1,kn+2)]
              Knots = as.numeric(quantile(u, u.k))  
              pi.u = bs(u, knots=Knots, intercept=TRUE, degree=degree)[,1:q]
              design.n.test = as.matrix(cbind(design.n.test,pi.u))
            }


            for (i in 1:m)
            {
              u = as.matrix(X[,i])
              u.k = seq(0, 1, length=kn+2)[-c(1,kn+2)]
              Knots = as.numeric(quantile(u, u.k))  
              pi.u = bs(u, knots=Knots, intercept=TRUE, degree=degree)[,1:(q)]
              xx = cbind(pi.u,as.matrix(X[,-i]))
              design.n = as.matrix(cbind(design.n,pi.u))
              y = as.matrix(Y[,i])        
              lasso.cv <- cv.glmnet(xx,y,alpha=1,nfolds=5)   
              alpha <- lasso.cv$lambda.min/3
              lasso.fit1 <- glmnet(xx,y,family="gaussian",alpha=1,nlambda=100,intercept=FALSE)
              coef0 <- as.matrix(as.vector(predict(lasso.fit1, s=alpha, type="coefficients"))[-1])  # sum((y[,1]-xx%*%coef0)^2)/n  
              beta0.n[,i] = coef0[1:q]
              beta0.p[,i][-i] = as.matrix(coef0[(q+1):(p+q-1)])
              res[,i] = y - as.matrix(X[,-i]) %*% as.matrix(beta0.p[,i][-i])- pi.u %*%as.matrix(beta0.n[,i]) 
            }

            res0 = res;

  

          for (i1 in 1:nlambda1)
             {       
               for (i3 in 1:nlambda3)
                { 
                   lam1 = lambda1[i1]    
                   lam2 = lambda1[i1]*q^(1/2)   
                   lam3 = lambda3[i3]          
                   beta.n = beta0.n
                   beta.p = beta0.p
                   res = res0;
                   diff.in1 = diff.out1 = diff.out2 = 1  
                   step.in1 = step.out1 = step.out2 = 0

                   while(diff.out1>=10^{-4} & step.out1<=30) 
                      {
                            beta.n.out1 = beta.n
                         for(k in 1:m) 
                              {
                                    y = as.matrix(Y[,k]); 
                                    sub = ((k-1)*q+1):(k*q)   
                                    old_beta.n = beta.n[sub]  
                                    y.sub = y-design.n[,-sub]%*%beta.n[-sub]
                                    x.sub = design.n[,sub]
                                    betak = beta.n[sub]
                                    normk = sqrt(sum((betak)^2)) 
                                    gwt2 = (max(normk,eps/100))^(-1)
                                    gwt = normk*sqrt(q)
                                    if (dmcp(gwt,gamma,lam2)>0)
                                      {
                                        wpp = rep(gwt2,q)
                                        beta.n[sub] = gALAL(x.sub, y.sub, tau, beta=0.9995, eps, wpp)$coefficients 
                                       }                                                          
                             } 
                          diff.out1=mean((beta.n-beta.n.out1)^2); 
                          step.out1=step.out1+1
                       } #  end of while loop 
                       beta.n[which(abs(beta.n)<eps)]=0

                    while(diff.out2>=10^{-4} & step.out2<=20) 
                      {
                            beta.p.out2 = beta.p
                         for(k in 1:m) 
                              {
                                   y = as.matrix(Y[,k]);
                                   y = y- design.n[,c(((k-1)*q+1):(k*q))]%*%beta.n[,k]
                                   for(i in 1:p) 
                                     {
                                        if (k != i) # update the parametric coefficients
                                        {

                                          if (i<p) 
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
                                                } # end of i<p


                                          if (i==p) # update parametric coefficients w.o. network structure
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
                               } 
                          diff.out2=mean((beta.p-beta.p.out2)^2); 
                          step.out2=step.out2+1
                       } #  end of while loop 
 
                   ### validation on an independent dataset 
                   res.test = matrix(0,nt,m)    
                   for (i.test in 1:m)
                        {
                            res.test[,i.test] = design.n.test[,c(((i.test-1)*q+1):(i.test*q))]%*%beta.n[,i.test]
                        }
                   res.test = Y.test - as.matrix(X.test) %*% as.matrix(beta.p) - res.test
                   for (i.test in 1:m) { cv.test[i1,i3] = cv.test[i1,i3] + sum(res.test[,i.test]*(tau-(res.test[,i.test]<0)))}
                   coef.p[,,i1,i3] = beta.p          
                   coef.n[,,i1,i3] = beta.n   
              } 
           } 


          cv = cv.test/nt

          i1.cv = which(cv == min(cv), arr.ind = TRUE)[1,1]
          i3.cv = which(cv == min(cv), arr.ind = TRUE)[1,2]
          flag.lambda1 = lambda1[i1.cv]          
          flag.lambda3 = lambda3[i3.cv]
       
          beta2.p = coef.p[,,i1.cv,i3.cv]  
          beta2.n = coef.n[,,i1.cv,i3.cv]

          bet.n = seq(1:200)
          bet.p = c(seq(6,10),seq(56,60),seq(101,200))
          bet.n[bet.p] = 0

          TP1 = sum((which(bet!=0) %in% which(beta2.p!=0))*1) 
          TP2 = sum((which(bet.n!=0) %in% which(colSums(abs(beta2.n))>tol))*1)
          TP = sum((which(bet!=0) %in% which(beta2.p!=0))*1) + sum((which(bet.n!=0) %in% which(colSums(abs(beta2.n))>tol))*1) 

          FP1 = sum((which(beta2.p!=0) %in% which(bet==0))*1)
          FP2 = sum( (which(colSums(abs(beta2.n))>tol) %in% which(bet.n==0))*1)
          FP = sum((which(beta2.p!=0) %in% which(bet==0))*1) + sum( (which(colSums(abs(beta2.n))>tol) %in% which(bet.n==0))*1)


        
tun = list(beta.p=beta2.p,beta.n=beta2.n,TP=TP, FP = FP, TP1=TP1, FP1=FP1,TP2=TP2, FP2=FP2)
return(tun)
}  




