source("covtest_simu_diffa_fun.R")

cl=makeCluster(599)
clusterEvalQ(cl,{library(mvtnorm);source("utils_lscm.R");NULL})
#clusterExport(cl,c("paras","df")) #clusterCall(cl,print,a+b)
clusterApplyLB(cl, 1:nrow(df), simu)
stopCluster(cl)
