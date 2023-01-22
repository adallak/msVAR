###########################################
## Model 4 example as in the manuscript
###########################################

rm(list = ls())

###########################################
## Please set your directory here
###########################################

#directory = ""
#setwd(directory)

source("msVAR_function.r")



genDiagA <- function(K, p = 3, block = 5) {
      ##This function generates coefficient 
      ##matrix A for model 4, as described  
      ##in the manuscript
      stopifnot(K%%block == 0)
      A = matrix(0, K, p*K)
      block_size = K / block
      c = 0
      for (lag in seq(p)) {
            for (k in seq(block)) {
                  left = c * K + 1
                  A[,left : (lag * K)][((k-1) * block_size + 1):(k * block_size),
                                       ((k-1)*block_size + 1):(k * block_size)] <- 
                        matrix(runif(block_size^2, -0.4, 0.4),
                               block_size, block_size)
            }
            c = c + 1
      }
      A = A/ (max(abs((svd(A)$d))))
      return(A)
}

K = 25   ## Data dimension
block = 5 ## Number of blocks
p = 3    ## Number of lags

###########################################
## Generate true coefficient matrix
###########################################

A = genDiagA(K = K, p = p, block = block)

## True adjacency matrix
adj = abs(A[,1:K])
adj[adj > 0] = 1


###########################################
## msVAR inputs
###########################################

size = 512
Sigma = diag(K)
burn = 300
halfWindowLength     = 12
lambda = 0.1
p.seq      		     = c(1,2,3)
p.ub       		     = max(p.seq)

###########################################
## Generate data
###########################################

set.seed(1223)
dta = simulateVAR(coefMatr=A, intercept =rep(0,K),
                  size=size, burn=burn,
                  Sigma=Sigma, error.dist="normal")$simData

##########################################
## Select lambda using BIC
##########################################

lambda = selectlambda(dta, gam = 0.5, thresh = 1e-6,
                      halfWindowLength = 18, rho = 10, 
                      Max_Iter = 100, diag = FALSE, 
                      verbose = TRUE,rho.flex = TRUE, 
                      alpha = 1.5, trim_max = 0.8, 
                      trim_min = 0.01,
                      freq = NULL, smooth = FALSE, 
                      criteria = "BIC")$lambda

#########################################
## Run msvar
#########################################

msVAR.result = msVAR(dta = dta, p.seq = p.seq, lambda = lambda,
                     halfWindowLength = halfWindowLength, 
                     stage1.showStatus = FALSE, 
                     stage2.showStatus = FALSE,           
                     standardize = FALSE, 
                     stage1.info = "bic", stage2.info = "bic",
                     iteMax = 300, reltol = 1e-3,
                     rho.flex = TRUE, rho = 10) 
msvar_A = msVAR.result$stage2$estA

#########################################
## Extract estimated adjacency matrix
#########################################

msvarAadj = matrix(1, K, K)
compAadj = abs(msvar_A)
lag = dim(msvar_A)[2] / K
if (lag > 1)
{
      for( i in 1 : K)
      {
            for( j in 1 : K)
            {
                  if (compAadj[i,j] == compAadj[(i), (j + K)])
                  {
                        if(compAadj[i,j] == 0)
                        {
                              msvarAadj[i, j] = 0
                        }
                  }
            }
      }
}else{
      msvarAadj = abs(msvar_A)
      msvarAadj[msvarAadj > 0] = 1
}

####################################
## Report TPR, FPR, TDR
####################################

compA = compareA(msvarAadj, adj)
print(compA)

####################################
## Plot the adjacency matrices
####################################

par(mfrow= c(1,2))
plot_mat(adj, main = "True")
plot_mat(msvarAadj, main = "msVAR")
par(mfrow= c(1,1))
