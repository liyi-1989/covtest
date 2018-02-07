library(CVXR)
 
GL=crossprod( fd(p) )

data = mvrnorm(n, mu = rep(0,p), Sigma = Sigma1)
Q=cov(data)

l1=0.05; l2=1000; alpha=100
#S <- Semidef(p)
S=Variable(p,p)
obj <- CVXR::norm(S-Q*(Sigma0==0),"F")+l1*norm1(S)#+l2*norm_nuc(S)##+  l1*matrix_trace(S%*%GL%*%t(S))
#constr <- list(matrix_trace(t(S)%*%GL%*%S) <= alpha)
prob <- Problem(Minimize(obj))
result <- solve(prob)

Sstar=result$getValue(S)
printM(Sstar)
printM(Sstar*(abs(Sstar)>1e-10))


base::norm(Q-Sigma1,"F")
base::norm(Sstar+Q*(Sigma0!=0)-Sigma1,"F")



