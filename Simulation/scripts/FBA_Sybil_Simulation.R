#LOADING BACKGROUND FUNCTIONS
source("Simulation/rawdata/source/source.R")

#LOADING DIET
NGM <- read.csv("Simulation/rawdata/data/NGM2020_gapseq.csv")
NGM$compounds <- paste0("EX_",NGM$compounds)
NGM$compounds <- paste0(NGM$compounds,"_e0")
NGM[c(35),3] <-  0.003  #enrich metals
NGM[c(12,13),3] <- 0.005 #enrich iron

##### FIRST PART: SINGLE-MODEL SUPPLEMENTATION #####
models <- readRDS("Simulation/rawdata/data/cembio-20201113_NGM-adapted.RDS")
models.biomass <- reconstrain.models.biomass(models = models, minerals = NGM)
checklist <- find.reactions.id.name(models.biomass)

#PREPARING THE DATAFRAME OF RESULTS
data.frame.change <- data.frame("react_id" = checklist$id,"Myb11_rel" = NA,
                                "Myb71_rel" = NA, "react_name" = checklist$name)
data.frame.change[195, 1] <- "ORIGINAL_LP"
data.frame.change[195, 2] <- optimizeProb(object = models.biomass[[1]],
                                          algorithm = "fba", retOptSol = F)$obj
data.frame.change[195, 3] <- optimizeProb(object = models.biomass[[2]],
                                          algorithm = "fba", retOptSol = F)$obj
data.frame.change$react_name <- gsub(pattern = "EX_cpd03561_e0",
                                     replacement = "2-Hydroxybutyrate-e0 Exchange",
                                     x = data.frame.change$react_name, fixed = T)
data.frame.change$react_name <- gsub(pattern = "EX_cpd00122_e0",
                                     replacement = "N-Acetyl-D-glucosamine-e0 Exchange",
                                     x = data.frame.change$react_name, fixed = T)
data.frame.change$react_name <- gsub(pattern = "EX_cpd00094_e0",
                                     replacement = "2-Oxobutyrate-e0 Exchange",
                                     x = data.frame.change$react_name, fixed = T)
data.frame.change$react_name <- gsub(pattern = "EX_cpd00380_e0",
                                     replacement = "Itaconate Exchange",
                                     x = data.frame.change$react_name, fixed = T)
data.frame.change$react_name <- gsub(pattern = "EX_cpd02351_e0",
                                     replacement = "Glucosaminate Exchange",
                                     x = data.frame.change$react_name, fixed = T)
data.frame.change$react_name <- gsub(pattern = "EX cpd11416 c0",
                                     replacement = "Biomass",
                                     x = data.frame.change$react_name, fixed = T)
data.frame.change$react_name <- gsub(pattern = "TRHL-e0 Exchange",
                                     replacement = "D-Trehalose-e0 Exchange",
                                     x = data.frame.change$react_name, fixed = T)
data.frame.change$react_name <- gsub(pattern = "PAN-e0 Exchange",
                                     replacement = "D-Pantothenate-e0 Exchange",
                                     x = data.frame.change$react_name, fixed = T)

#SUPPLEMENTATION AND SAVING OF THE RESULT
for (r in 1:2) {
  for (t in 1:length(checklist$id)) { 
    original <- models.biomass[[r]]
    lb_original <-  models.biomass[[r]]@lowbnd[grep(pattern = checklist$id[t],
                                                    x = models.biomass[[r]]@react_id)]
    if (length(lb_original) > 0) {
      changedbounds <- changeBounds(model = original, lb = lb_original - 10,
                                    react = checklist$id[t])
      new <- optimizeProb(object = changedbounds, algorithm = "fba", retOptSol = F)$obj
      rchange <- (round(new,6) - round(data.frame.change[195, r + 1],6))/(round(data.frame.change[195, r + 1],6))
      data.frame.change[t, r + 1] <- rchange}
  }
}
#save data frame for the supplementary file S1
#write.table(data.frame.change, "S1FBA.csv", sep = ",", row.names = F) #SEE S3 TAB
data.frame.change <- data.frame.change[which(is.na(data.frame.change[,4]) == F),] # REMOVE ORIGINAL OPTIMISATION
onlymyb11 <- data.frame.change[which(data.frame.change[,2] > 0.01),] # POSITIVE EFFECT ON MYB11
onlymyb11 <- onlymyb11[which(onlymyb11[,2] > onlymyb11[,3]),] # EFFECT ON MYB11 GREATER THAN ON MYB71

onlymyb11$rank <- rank(onlymyb11$Myb11_rel,ties.method = "average") # RANKING

saveRDS(onlymyb11,"Simulation/rawdata/3_1.RDS")


