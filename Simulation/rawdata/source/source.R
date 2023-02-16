library(sybil)
library(sybilSBML)
library(cplexAPI)
library(MicrobiomeGS2)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
#source("https://raw.githubusercontent.com/Waschina/sybil.diagnostics/master/getModelDiagnostics.R")
#The object iJO1366 is not related with model "http://bigg.ucsd.edu/models/iJO1366", apart from that
#it is only used to give the notation that the object is a model

reconstrain.models.biomass <- function(models, minerals){
  for (e in 1:length(models)) {
    iJO1366 <- models[[e]]
    reactions_viJO1366 <- grep(pattern = "EX_",
                               iJO1366@react_id, value = T, fixed = T)
    reactions_piJO1366 <- grep(pattern = "EX_",
                               iJO1366@react_id, value = F, fixed = T)
    if (is.null(minerals) == F) {
    min_p <- minerals$compounds %in% reactions_viJO1366
    
    iJO1366 <- changeBounds(model = iJO1366, lb = 0,
                            react = iJO1366@react_id[reactions_piJO1366])
    iJO1366 <- changeBounds(model = iJO1366, lb = -minerals$maxFlux[min_p],
                            react = minerals$compounds[min_p])}
    
    iJO1366@obj_coef <- rep(0, length(iJO1366@react_id))
    
    iJO1366@obj_coef[grep("EX_cpd11416_c0",iJO1366@react_id)] <- 1
    
    models[[e]] <- iJO1366
    iJO1366 <- NA
  }
  return(models)
}


find.react_al_vec <- function(models){
  react_al_vec <- vector()
  for (e in 1:length(models)) {
    iJO1366 <- models[[e]]
    reactions_viJO1366 <- grep(pattern = "EX_",
                               iJO1366@react_id, value = T, fixed = T)
    react_al_vec <- unique(c(react_al_vec,reactions_viJO1366))}
  return(react_al_vec)}

make.df.react <- function(models){
  react_al_vec <- find.react_al_vec(models)
  df.opt <- data.frame("model" = NA, "original_lp" = NA,stringsAsFactors = F)
  w <- 1
  df.opt[1:length(names(models)),w] <- names(models)
  w <- 2
  for (r in 1:length(react_al_vec)) {
    w <- w + 1 
    df.opt[,w] <- NA
    colnames(df.opt)[w] <- react_al_vec[[r]]}
  return(df.opt)
}


find.reactions.id.name <- function(models){
  reac <- find.react_al_vec(models)
dd <- data.frame(row = 1:length(reac), id = NA, name = NA)
for (d in 1:length(reac)) {
  h <- reac[[d]]
for (k in 1:length(models)) {
  if (length(grep(pattern = h,x = models[[k]]@react_id)) > 0) {
    p1 <- models[[k]]@react_name[grep(pattern = h,x = models[[k]]@react_id)]
  }
}
  dd[d,2] <- h
  dd[d,3] <- p1}
return(dd)}

joined.fba.bio.plus <- function(models,plus, name.join, df.exist=NULL) {
  if (is.null(df.exist) == T) {
    react_al_vec <- find.react_al_vec(models)
    df.opt <- make.df.react(models)}
  df.opt[3,] <- NA
  df.opt[3,1] <- name.join
  k = 2
  mod.joined <- join_mult_models(model.list = models)
  
  mod.joined$modj@obj_coef <- rep(0, length(mod.joined$modj@react_id))
  mod.joined$modj@obj_coef[grep("M[0-9]+_EX_cpd11416_c0",mod.joined$modj@react_id)] <- 1
  
  mod.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11416",mod.joined$modj@react_id)] <- 1000
  mod.joined$modj@lowbnd[grep("M[0-9]+_EX_cpd11416",mod.joined$modj@react_id)] <- 0

  # block quinon exchanges
  mod.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11606",mod.joined$modj@react_id)] <- 0.5
  mod.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11451",mod.joined$modj@react_id)] <- 0.5
  
  coupling <- get_coupling_constraints_mult(mod.joined$modj, cpl_c = 400,
                                            cpl_u = 0.01)
  
  modj_warm <- sysBiolAlg(mod.joined$modj,
                          algorithm = "mtfEasyConstraint2",
                          easyConstraint = coupling,
                          pFBAcoeff = 1e-6,
                          scaling = 1)
  mod.joined$solj <- optimizeProb(modj_warm)
  if (mod.joined$solj$stat == 1) {
    df.opt[1:2,k] <- round(mod.joined$solj$fluxes[grep("M[0-9]+_EX_cpd11416_c0",mod.joined$modj@react_id,value = F)],6)
    df.opt[3,k]  <-  mod.joined$solj$obj  
  } else {df.opt[,k] <- c(0,0,0)}
  
  for (k in 3:length(colnames(df.opt))) {
    print(paste0(k,"/",length(colnames(df.opt))))
    PLUS <- models
    for (e in 1:length(PLUS)) {
      check_model <- PLUS[[e]]
      R1 <- which((check_model@react_id == colnames(df.opt)[[k]]) == T)
      if (length(R1) == 0) {next()}
      if (length(R1) > 0) {
        
        iJO1366_new <- changeBounds(model = check_model,
                                    lb = check_model@lowbnd[R1] - plus,
                                    ub = 1000,
                                    react = check_model@react_id[R1])
        PLUS[[e]] <- iJO1366_new
      }
    }
    
    mod.joined <- join_mult_models(model.list = PLUS)
    
    mod.joined$modj@obj_coef <- rep(0, length(mod.joined$modj@react_id))
    mod.joined$modj@obj_coef[grep("M[0-9]+_EX_cpd11416_c0",mod.joined$modj@react_id)] <- 1
    
    mod.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11416",mod.joined$modj@react_id)] <- 1000
    mod.joined$modj@lowbnd[grep("M[0-9]+_EX_cpd11416",mod.joined$modj@react_id)] <- 0
    
    # block quinon exchanges
    mod.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11606",mod.joined$modj@react_id)] <- 0.5
    mod.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11451",mod.joined$modj@react_id)] <- 0.5
    
    coupling <- get_coupling_constraints_mult(mod.joined$modj, cpl_c = 400,
                                              cpl_u = 0.01)
    
    modj_warm <- sysBiolAlg(mod.joined$modj,
                            algorithm = "mtfEasyConstraint2",
                            easyConstraint = coupling,
                            pFBAcoeff = 1e-6,
                            scaling = 1)
    mod.joined$solj <- optimizeProb(modj_warm)
    if (mod.joined$solj$stat == 1) {
      df.opt[1:2,k] <- round(mod.joined$solj$fluxes[grep("M[0-9]+_EX_cpd11416_c0",mod.joined$modj@react_id,value = F)],6)
      df.opt[3,k]  <-  mod.joined$solj$obj  
    } else {df.opt[,k] <- c(0,0,0)}
  } 
  return(df.opt)
}
  
