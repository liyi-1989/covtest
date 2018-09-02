source("covtest_simu_multi_fun.R")

cl=makeCluster(1)
clusterEvalQ(cl,{library(mvtnorm);source("utils_lscm.R")})
clusterExport(cl,c("paras","df")) #clusterCall(cl,print,a+b)
clusterApplyLB(cl, 1:nrow(df), simu, paras,df)
stopCluster(cl)



for(i in 1:nrow(df)){
  print(i)
  simu(i,paras,df)
}