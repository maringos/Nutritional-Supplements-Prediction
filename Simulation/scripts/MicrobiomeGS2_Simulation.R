#LOADING BACKGROUND FUNCTIONS
source("Simulation/rawdata/source/source.R")

#LOADING DIET
NGM <- read.csv("Simulation/rawdata/data/NGM2020_gapseq.csv")
NGM$compounds <- paste0("EX_",NGM$compounds)
NGM$compounds <- paste0(NGM$compounds,"_e0")
NGM[c(35),3] <-  0.003  #enrich metals
NGM[c(12,13),3] <- 0.005 #enrich iron

###### SECOND PART: MULTI-MODEL SUPPLEMENTATION #####
#LOADING MODELS
models <- readRDS("Simulation/rawdata/data/cembio-20201113_NGM-adapted.RDS")
models.biomass <- reconstrain.models.biomass(models = models, minerals = NGM)
checklist <- find.reactions.id.name(models.biomass)
JOIN.BIO_biolog <- joined.fba.bio.plus(models = models.biomass, plus = 10,name.join = "MYb71 x MYb11",
                                       df.exist = NULL) # SUPPLEMENTATION EXPERIMENTS

#NORMALIZATION SO THAT TOTAL OF GROWTH IS ALWAYS 1
for (k in 2:length(colnames(JOIN.BIO_biolog))) {
  JOIN.BIO_biolog[1,k] <- JOIN.BIO_biolog[1,k]/JOIN.BIO_biolog[3,k]
  JOIN.BIO_biolog[2,k] <- JOIN.BIO_biolog[2,k]/JOIN.BIO_biolog[3,k]
  JOIN.BIO_biolog[3,k] <- JOIN.BIO_biolog[3,k]/JOIN.BIO_biolog[3,k]
  
}

#save data frame for the supplementary file S4
write.table(JOIN.BIO_biolog, "S1MGS2.csv", sep = ",", row.names = F)

#EXTRACTING RESULTS 
#FORMATING
data.frame.change1 <- data.frame("react_id" = colnames(JOIN.BIO_biolog)[3:196],
                                 "Myb11_rel" = NA, "Myb71_rel" = NA,
                                 "react_name" = NA)
for (k in 3:length(colnames(JOIN.BIO_biolog))) {
  data.frame.change1$react_name[[k - 2]] <- checklist$name[grep(pattern = colnames(JOIN.BIO_biolog)[[k]],
                                                                x = checklist$id,value = F)[1]]
}

#COMPOUNDS THAT SUPPORT MYB11
for (k in 3:length(colnames(JOIN.BIO_biolog))) {
  #RECORD CHANGES DUE TO SUPPLEMENTATION
  rchange <- (round(JOIN.BIO_biolog[1,k],6) - round(JOIN.BIO_biolog[1,2],6))/round((JOIN.BIO_biolog[1,2]),6)
  data.frame.change1$Myb11_rel[[k - 2]] <- rchange
  rchange <- (round(JOIN.BIO_biolog[2,k],6) - round(JOIN.BIO_biolog[2,2],6))/round((JOIN.BIO_biolog[2,2]),6)
  data.frame.change1$Myb71_rel[[k - 2]] <- rchange
}

# REMOVE NA IF ANY
data.frame.change1 <- data.frame.change1[which(is.na(data.frame.change1[,4]) == F),]

onlymyb11comb <- data.frame.change1[which(data.frame.change1[,2] > 0.01),] # POSITIVE EFFECT ON MYB11
onlymyb11comb <- onlymyb11comb[which(onlymyb11comb[,2] > onlymyb11comb[,3]),] # EFFECT ON MYB11 GREATER THAN ON MYB71

onlymyb11comb$rank <- rank(onlymyb11comb$Myb11_rel,ties.method = "average") # RANKING
saveRDS(onlymyb11comb,"Simulation/rawdata/3_2.RDS")
saveRDS(JOIN.BIO_biolog,"Simulation/rawdata/3_2_alternative.RDS")
