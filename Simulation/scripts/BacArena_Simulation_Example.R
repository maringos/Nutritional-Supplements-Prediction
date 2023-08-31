#THIS IS AN EXAMPLE OF BACARENA SIMULATIONS
library(BacArena)
library(parallel)
setwd("Simulation/rawdata/")
source("source/source.R")
#DIET
NGM <- read.csv("data/NGM2020_gapseq.csv")
NGM$compounds <- paste0("EX_",NGM$compounds)
NGM$compounds <- paste0(NGM$compounds,"_e0")
NGM[c(35),3] <- 0.003 #enrich metals
NGM[c(12,13),3] <- 0.005 #enrich iron
#MODELS
models <- readRDS("data/cembio-20201113_NGM-adapted.RDS")
reactions <- find.react_al_vec(models) # CANDIDATE COMPOUNDS
#ORIGINAL PARALLEL SIMULATION BACARENA
replicates <- 15
cores <- 15
cl <- makeCluster(cores, type = "PSOCK")
clusterExport(cl, c("models","NGM","reactions"))
clusterEvalQ(cl, sink(paste0("*****", Sys.getpid(), ".txt")))
  simlist <- parLapply(cl, 1:(replicates), function(i){
  library(BacArena)
  library(parallel)
#BACARENA SIMULATION SETTINGS   
arena <- Arena(n = 30, m = 30)
    bac1 <- Bac(models[[1]])
  bac2 <- Bac(models[[2]])
arena <- addOrg(object = arena, specI = bac1,amount = 20)
arena <- addOrg(object = arena, specI = bac2,amount = 20)
  arena <- addSubs(object = arena, smax = 0.001*NGM$maxFlux, 
                   mediac = NGM$compounds, unit  = "mM")
arena@tstep <- 1
sim <- simEnv(arena, time = 12, sec_obj = "mtf")
})
saveRDS(simlist, paste0('****/sim_biolog_original.RDS'))
stopCluster(cl) 

#SUPPLEMENTATION SIMULATION BACARENA
for (j in 1:20) { # EXAMPLE: ADDING REACTIONS FROM 1 TO 20
  cl <- makeCluster(cores, type = "PSOCK")
  clusterExport(cl, c("models","NGM","reactions","j"))
  clusterEvalQ(cl, sink(paste0("****", Sys.getpid(), ".txt")))
  simlist <- parLapply(cl, 1:(replicates), function(i){
  library(BacArena)
  library(parallel)
  #BACARENA SIMULATION SETTINGS   
  arena <- Arena(n = 30, m = 30)
    bac1 <- Bac(models[[1]])
  bac2 <- Bac(models[[2]])
  arena <- addOrg(object = arena, specI = bac1,amount = 20)
  arena <- addOrg(object = arena, specI = bac2,amount = 20)
  arena <- addSubs(object = arena, smax = 0.001*NGM$maxFlux,
                   mediac = NGM$compounds, unit = "mM")
  arena1 <- arena
  arena1 <- addSubs(arena1,smax = 0.001*10, add = T, mediac = reactions[[j]],
                    unit = "mM")
  arena1 <- simEnv(arena1, time = 12, sec_obj = "mtf")})
  saveRDS(simlist, paste0('*****sim_biolog_',reactions[[j]],'.RDS'))
  stopCluster(cl)
  }

