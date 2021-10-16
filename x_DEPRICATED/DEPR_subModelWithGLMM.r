#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1]: LogReg formated lineages, DB_directory, lowest reLevel, highest rank 

# try evaluating models with performance_accurcacy() AUC analysis from Performance package
# print some model estimates, need to look at those and think about whether we're using the most logical specification

library(MuMIn)
library(glmmTMB)

CV <- read.csv(args[1], header=FALSE, na.strings=c("","NA"))
CV$s <- ifelse(as.character(CV$V8)!=as.character(CV$V15) & !is.na(CV$V8) & !is.na(CV$V15), 1, 0)
CV$g <- ifelse(as.character(CV$V7)!=as.character(CV$V14) & !is.na(CV$V7) & !is.na(CV$V14), 1, 0)
CV$f <- ifelse(as.character(CV$V6)!=as.character(CV$V13) & !is.na(CV$V6) & !is.na(CV$V13), 1, 0)
CV$o <- ifelse(as.character(CV$V5)!=as.character(CV$V12) & !is.na(CV$V5) & !is.na(CV$V12), 1, 0)
if (args[4] == "Class") {
	CV$c <- ifelse(as.character(CV$V4)!=as.character(CV$V11) & !is.na(CV$V4) & !is.na(CV$V11), 1, 0)
}
if (args[4] == "Phylum") {
	CV$p <- ifelse(as.character(CV$V3)!=as.character(CV$V10) & !is.na(CV$V3) & !is.na(CV$V10), 1, 0)
}
if (args[4] == "Kingdom") {
	CV$k <- ifelse(as.character(CV$V2)!=as.character(CV$V9) & !is.na(CV$V2) & !is.na(CV$V9), 1, 0)
}

#g,k, p, c, o, f,g,s,i,l
#1,2, 3, 4, 5, 6,7,8,9,10
#     p, c, o, f, , 
#    10,11,12,13
# 17: alignLength, 18: alignID

CV$V18 <- CV$V18^(1/2)

#print(head(CV))

# running kingdom modelling
if (args[4] == "Kingdom") {
	cat("\n"); print('producing kingdom model'); cat("\n")
	CV <- subset(CV, CV$V18 >= 10)
	Kmod <- glmmTMB(k ~ V18 + V17 + (1|V2), data = CV, family = binomial()) 
	save(Kmod, file = paste(args[2], '/', "KingdomGLMM.rda", sep=''))
	print(summary(Kmod))
	print(r.squaredGLMM(Kmod))
}

# running phylum modelling
if (args[4] == "Phylum" | args[4] == "Kingdom") {
	cat("\n"); print('producing phylum model'); cat("\n")
	CV <- subset(CV, CV$V18 >= 10)
	Pmod <- glmmTMB(p ~ V18 + V17 + (1|V2/V3), data = CV, family = binomial()) 
	save(Pmod, file = paste(args[2], '/', "PhylumGLMM.rda", sep=''))
	print(summary(Pmod))
	print(r.squaredGLMM(Pmod))
}

# running class modelling
if (args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom") {
	cat("\n"); print('producing class model'); cat("\n")
	CV <- subset(CV, CV$V18 >= 20^(1/2))
	Cmod <- glmmTMB(c ~ V18 + V17 + (1|V3/V4), data = CV, family = binomial()) # nest to class
	save(Cmod, file = paste(args[2], '/', "ClassGLMM.rda", sep=''))
	print(summary(Cmod))
	print(r.squaredGLMM(Cmod))
}

# running order modelling
cat("\n"); print('producing order model'); cat("\n")
CV <- subset(CV, CV$V18 >= 60^(1/2))
Omod <- glmmTMB(o ~ V18 + (1|V4/V5), data = CV, family = binomial()) # nest to order ### + V17
save(Omod, file = paste(args[2], '/', "OrderGLMM.rda", sep=''))
print(summary(Omod))
print(r.squaredGLMM(Omod))

# running family modelling
cat("\n"); print('producing family model'); cat("\n")
CV <- subset(CV, CV$V18 >= 70^(1/2))
Fmod <- glmmTMB(f ~ V18 + (1|V4/V5/V6), data = CV, family = binomial()) # nest glmm down to family
save(Fmod, file = paste(args[2], '/', "FamilyGLMM.rda", sep=''))
print(summary(Fmod))
print(r.squaredGLMM(Fmod))

# running Genus modelling
if (args[3] == "Family") {
	cat("\n"); print('producing genus model w/ family as lowest random intercept term'); cat("\n")
	CV <- subset(CV, CV$V18 >= 80^(1/2))
	Gmod <- glmmTMB(g ~ V18 + (1|V4/V5/V6), data = CV, family = binomial()) 
	save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep=''))
	print(summary(Gmod))
	print(r.squaredGLMM(Gmod))
}
if (args[3] == "Genus" | args[3] == "Species") {
	cat("\n"); print('producing genus model w/ genus as lowest random intercept term'); cat("\n")
	CV <- subset(CV, CV$V18 >= 80^(1/2))
	Gmod <- glmmTMB(g ~ V18 + (1|V4/V5/V6/V7), data = CV, family = binomial())
	save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep=''))
	print(summary(Gmod))
	print(r.squaredGLMM(Gmod))
}

# running species modelling
if (args[3] == "Family") {
	cat("\n"); print('producing species model w/ family as lowest random intercept term'); cat("\n\n")
	CV <- subset(CV, CV$V18 >= 85^(1/2))
	Smod <- glmmTMB(s ~ V18 + (1|V4/V5/V6), data = CV, family = binomial())
	save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep=''))
	print(summary(Smod))
	print(r.squaredGLMM(Smod))
}
if (args[3] == "Genus") {
	cat("\n"); print('producing species model w/ genus as lowest random intercept term'); cat("\n\n")
	CV <- subset(CV, CV$V18 >= 85^(1/2))
	Smod <- glmmTMB(s ~ V18 + (1|V4/V5/V6/V7), data = CV, family = binomial())
	save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep=''))
	print(summary(Smod))
	print(r.squaredGLMM(Smod))
}
if (args[3] == "Species") {
	cat("\n"); print('producing species model w/ species as lowest random intercept term'); cat("\n\n")
	CV <- subset(CV, CV$V18 >= 85^(1/2))
	Smod <- glmmTMB(s ~ V18 + (1|V5/V6/V7/V8), data = CV, family = binomial())
	save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep=''))
	print(summary(Smod))
	print(r.squaredGLMM(Smod))
}

## default to family but include options for genus or species level random effect implementations
