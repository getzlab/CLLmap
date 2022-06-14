# Parallel computing utils
# Author: Binyamin Knisbacher

# Rlib_dir = "/private/common/Software/R/Rlibs"
# .libPaths(Rlib_dir)
# install.packages("doParallel")
install_repos = "http://cran.us.r-project.org"
if("doParallel" %in% rownames(installed.packages()) == FALSE) {install.packages("doParallel", repos=install_repos , quiet=T)}


library(doParallel, quietly = T) #for backend registering
library(foreach, quietly = T) #for parallel for loops
# library(parallel) #for parallel

parallel.registerCluster <- function(ncores2use=8, retCL=F){
  # Find out how many cores are available (if you don't already know)
  total_cores = detectCores()
  # Create cluster with desired number of cores
  cl <- makeCluster(min(ncores2use, total_cores))
  # Register cluster
  registerDoParallel(cl)
  # Find out how many cores are being used
  # getDoParWorkers()
  if(retCL){
    return(cl)
  }
}

#To remove all registered clusters
parallel.unregisterClusters <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
parallel.unregisterClusters() #unregister clusters upon activation
