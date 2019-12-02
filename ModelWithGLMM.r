#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1]: LogReg formated lineages, DB_directory 

library(MuMIn)
library(glmmTMB)

CV <- read.csv(args[1], header=FALSE, na.strings=c("","NA"))
CV$s <- ifelse(as.character(CV$V8)!=as.character(CV$V15) & !is.na(CV$V8) & !is.na(CV$V15), 1, 0)
CV$g <- ifelse(as.character(CV$V7)!=as.character(CV$V14) & !is.na(CV$V7) & !is.na(CV$V14), 1, 0)
CV$f <- ifelse(as.character(CV$V6)!=as.character(CV$V13) & !is.na(CV$V6) & !is.na(CV$V13), 1, 0)
CV$o <- ifelse(as.character(CV$V5)!=as.character(CV$V12) & !is.na(CV$V5) & !is.na(CV$V12), 1, 0)
CV$c <- ifelse(as.character(CV$V4)!=as.character(CV$V11) & !is.na(CV$V4) & !is.na(CV$V11), 1, 0)
#CV$p <- ifelse(as.character(CV$V3)!=as.character(CV$V10) & !is.na(CV$V3) & !is.na(CV$V10), 1, 0)
#CV$k <- ifelse(as.character(CV$V2)!=as.character(CV$V9) & !is.na(CV$V2) & !is.na(CV$V9), 1, 0)

# running kingdom modelling
#cat("\n")
#print('producing kingdom model')
#cat("\n")
#CV <- subset(CV, CV$V18 >= 10)
#Kmod <- glmmTMB(k ~ V18 + V17 + (1|V2), data = CV, family = binomial()) 
#save(Kmod, file = paste(args[2], '/', "KingdomGLMM.rda", sep=''))
#print(summary(Kmod))
#print(r.squaredGLMM(Kmod))

# running phylum modelling
#print('producing phylum model')
#cat("\n")
#CV <- subset(CV, CV$V18 >= 10)
#Pmod <- glmmTMB(p ~ V18 + V17 + (1|V2/V3), data = CV, family = binomial()) 
#save(Pmod, file = paste(args[2], '/', "PhylumGLMM.rda", sep=''))
#print(summary(Pmod))
#print(r.squaredGLMM(Pmod))

# running class modelling
print('producing class model')
cat("\n")
CV <- subset(CV, CV$V18 >= 20)
Cmod <- glmmTMB(c ~ V18 + V17 + (1|V10/V11), data = CV, family = binomial()) # nest to class
save(Cmod, file = paste(args[2], '/', "ClassGLMM.rda", sep=''))
print(summary(Cmod))
print(r.squaredGLMM(Cmod))

# running order modelling
print('producing order model')
cat("\n")
CV <- subset(CV, CV$V18 >= 75)
Omod <- glmmTMB(o ~ V18 + V17 + (1|V11/V12), data = CV, family = binomial()) # nest to order
save(Omod, file = paste(args[2], '/', "OrderGLMM.rda", sep=''))
print(summary(Omod))
print(r.squaredGLMM(Omod))

# running family modelling
print('producing family model')
cat("\n")
CV <- subset(CV, CV$V18 >= 88)
Fmod <- glmmTMB(f ~ V18 + V17 + (1|V11/V12/V13), data = CV, family = binomial()) # nest glmm down to family
save(Fmod, file = paste(args[2], '/', "FamilyGLMM.rda", sep=''))
print(summary(Fmod))
print(r.squaredGLMM(Fmod))

# running Genus modelling
print('producing genus model')
cat("\n")
CV <- subset(CV, CV$V18 >= 92)
Gmod <- glmmTMB(g ~ V18 + V17 + (1|V11/V12/V13), data = CV, family = binomial()) 
save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep=''))
print(summary(Gmod))
print(r.squaredGLMM(Gmod))

# running species modelling
print('producing species model')
cat("\n\n")
CV <- subset(CV, CV$V18 >= 95)
Smod <- glmmTMB(s ~ V18 + V17 + (1|V11/V12/V13), data = CV, family = binomial())
save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep=''))
print(summary(Smod))
print(r.squaredGLMM(Smod))
