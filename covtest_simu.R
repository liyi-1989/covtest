source("covtest_simu_fun.R")

cl=makeCluster(288)
clusterEvalQ(cl,{library(mvtnorm);source("utils_lscm.R")})
clusterExport(cl,c("paras","df")) #clusterCall(cl,print,a+b)
clusterApplyLB(cl, 1:nrow(df), simu, paras,df)
stopCluster(cl)
