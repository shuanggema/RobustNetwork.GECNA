
          rm(list=ls(all=TRUE))
          # setwd("/home2/work")


          source("generator.r")
          source("functions.r")
          source("A5.r")
          source("A6.r")
          source("A7.r")
          source("A8.r")
          library(mvtnorm)      
          library(splines)
          library(glmnet)
          library(quantreg)



          n = 200;tol=1e-04; eps=1e-06;
          m = p = 200;  
          reps = 100; nt=n;
          bet <- as.matrix(read.csv("beta.csv",sep=",",header=F))
          Vstr = "Block"; Vrho = 0.9; 

          tp.n1 = fp.n1 = tp.n2 = fp.n2 =  tp.n3 = fp.n3 = tp.n4 = fp.n4 =rep(0,reps)
          tp1.n1 = fp1.n1 = tp1.n2 = fp1.n2 = tp1.n3 = fp1.n3 = tp1.n4 = fp1.n4 = rep(0,reps)
          tp2.n1 = fp2.n1 = tp2.n2 = fp2.n2 = tp2.n3 = fp2.n3 = tp2.n4 = fp2.n4 = rep(0,reps)
          

          for (i in 1:reps)
                {
                  data = generator(m, n, p, bet, Vstr, Vrho)
                  Y = data$Y; X = data$X;  

                  # generate an independent testing set
                  data.test = generator(m, n, p, bet, Vstr, Vrho)
                  Y.test = data.test$Y; X.test = data.test$X; 
 
                   res.n1 = A5(Y,bet,X,Y.test,X.test,tol,eps)   # A5
                   res.n2 = A6(Y,bet,X,Y.test,X.test,tol,eps)   # A6
                   res.n3 =  A7(Y,bet,X,Y.test,X.test,tol,eps)   # A7
                   res.n4 =  A8(Y,bet,X,Y.test,X.test,tol,eps)   # A8

                     ## A5
                     tp.n1[i] = res.n1$TP
                     fp.n1[i] = res.n1$FP
                     tp1.n1[i] = res.n1$TP1
                     fp1.n1[i] = res.n1$FP1
                     tp2.n1[i] = res.n1$TP2
                     fp2.n1[i] = res.n1$FP2


                     ## A6
                     tp.n2[i] = res.n2$TP
                     fp.n2[i] = res.n2$FP 
                     tp1.n2[i] = res.n2$TP1
                     fp1.n2[i] = res.n2$FP1 
                     tp2.n2[i] = res.n2$TP2
                     fp2.n2[i] = res.n2$FP2 


                     ## A7
                     tp.n3[i] = res.n3$TP
                     fp.n3[i] = res.n3$FP 
                     tp1.n3[i] = res.n3$TP1
                     fp1.n3[i] = res.n3$FP1 
                     tp2.n3[i] = res.n3$TP2
                     fp2.n3[i] = res.n3$FP2 
           
                     ## A8
                     tp.n4[i] = res.n4$TP
                     fp.n4[i] = res.n4$FP 
                     tp1.n4[i] = res.n4$TP1
                     fp1.n4[i] = res.n4$FP1 
                     tp2.n4[i] = res.n4$TP2
                     fp2.n4[i] = res.n4$FP2 

                     cat(i, "\n")

                }
                   




