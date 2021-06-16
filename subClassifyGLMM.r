#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1]: Metaxa2 formated tax lineages, args[2] dabatase path, args[3] outfile, args[4] highest rank analyzed 

library(MuMIn)
library(glmmTMB)

if (args[4] == "Kingdom") {
	load(paste(args[2], "KingdomGLMM.rda", sep=''))
	load(paste(args[2], "PhylumGLMM.rda", sep='')) }
if (args[4] == "Phylum") {
        load(paste(args[2], "PhylumGLMM.rda", sep='')) }
if (args[4] == "Class") {
	load(paste(args[2], "ClassGLMM.rda", sep='')) }
load(paste(args[2], "OrderGLMM.rda", sep=''))
load(paste(args[2], "FamilyGLMM.rda", sep=''))
load(paste(args[2], "GenusGLMM.rda", sep=''))
load(paste(args[2], "srGenusGLMM.rda", sep=''))
load(paste(args[2], "SpeciesGLMM.rda", sep=''))
load(paste(args[2], "srSpeciesGLMM.rda", sep=''))

df <- read.csv(args[1], header=FALSE, na.strings=c("","NA"))
colnames(df) <- c('GI','kReal','pReal','cReal','oReal','fReal','gReal','sReal','na','TopID','Length','oScndID','fScndID','gScndID','sScndID','oTpHts','fTpHts','gTpHts','sTpHts')

# Transform TopID if specified
df$TopID <- df$TopID^(1/2)

# mark hits with only one hit to the top taxon
df$sTpHts <- ifelse(df$sTpHts>1&df$sTpHts<4,2,df$sTpHts); df$sTpHts <- ifelse(df$sTpHts>3,3,df$sTpHts); df$sTpHts <- as.factor(df$sTpHts)
df$gTpHts <- ifelse(df$gTpHts>1&df$gTpHts<4,2,df$gTpHts); df$gTpHts <- ifelse(df$gTpHts>3,3,df$gTpHts); df$gTpHts <- as.factor(df$gTpHts)
df$fTpHts <- ifelse(df$fTpHts>1&df$fTpHts<4,2,df$fTpHts); df$fTpHts <- ifelse(df$fTpHts>3,3,df$fTpHts); df$fTpHts <- as.factor(df$fTpHts)
df$oTpHts <- ifelse(df$oTpHts>1&df$oTpHts<4,2,df$oTpHts); df$oTpHts <- ifelse(df$oTpHts>3,3,df$oTpHts); df$oTpHts <- as.factor(df$oTpHts)

# if second hit option is RandomEffect, bin ScndID distances
df$sDist <- (df$TopID - df$sScndID)
df$gDist <- (df$TopID - df$gScndID)
df$fDist <- (df$TopID - df$fScndID)
df$oDist <- (df$TopID - df$oScndID)
df$sDist <- ifelse(df$sDist<1,1,df$sDist); df$sDist <- ifelse(df$sDist>1,2,df$sDist); df$sDist <- as.factor(df$sDist)
df$gDist <- ifelse(df$gDist<1,1,df$gDist); df$gDist <- ifelse(df$gDist>1,2,df$gDist); df$gDist <- as.factor(df$gDist)
df$fDist <- ifelse(df$fDist<1,1,df$fDist); df$fDist <- ifelse(df$fDist>1,2,df$fDist); df$fDist <- as.factor(df$fDist)
df$oDist <- ifelse(df$oDist<1,1,df$oDist); df$oDist <- ifelse(df$oDist>1,2,df$oDist); df$oDist <- as.factor(df$oDist)

# make new dataframe to hold results
#NewDF <- as.data.frame(df[,1]); colnames(NewDF) <- c('Accession')

#NewDF$Kingdom <- df$kReal
#if (args[4] == "Kingdom") {
#	NewDF$pKingdom <- round((1 - predict(Kmod, df, type = "response",allow.new.levels=TRUE)),2)
#} else { NewDF$pKingdom <- NA }

#NewDF$Phylum <- df$pReal
#if (args[4] == "Phylum") {
#	NewDF$pPhylum <- round((1 - predict(Pmod, df, type = "response",allow.new.levels=TRUE)),2)
#} else { NewDF$pPhylum <- NA }
#
#NewDF$Class <- df$cReal
#if (args[4] == "Class") {
#	NewDF$pClass <- round((1 - predict(Cmod, df, type = "response",allow.new.levels=TRUE)),2)
#} else { NewDF$pClass <- NA }
#
#NewDF$Order <- df$oReal
#NewDF$pOrder <- round((1 - predict(Omod, df, type = "response",allow.new.levels=TRUE)),2)
#
#NewDF$Family <- df$fReal
#NewDF$pFamily <- round((1 - predict(Fmod, df, type = "response",allow.new.levels=TRUE)),2)

# if -re == Family
#if (args[5] == "Family") {
#	NewDF$Genus <- df$gReal
#	NewDF$pGenus <- round((1 - predict(Gmod, df, type = "response",allow.new.levels=TRUE)),2)
#	NewDF$pGenus <- ifelse(NewDF$pFamily > 0.95, NewDF$pGenus, NA)
#	NewDF$Species <- df$sReal
#	NewDF$pSpecies <- round((1 - predict(Smod, df, type = "response",allow.new.levels=TRUE)),2)
#	NewDF$pSpecies <- ifelse(NewDF$pGenus>0.95, NewDF$pSpecies, NA)}

# if -re == SpeedGenus
if (args[5] == "SpeedGenus") { cat("\n"); print('Loading data'); cat("\n")
	srG <- read.csv(file = paste(args[2], "/srGenera.csv", sep=''), header = TRUE)
	DF <- df[!(df$gReal%in%srG[,2]),]
	srDF <- df[(df$gReal%in%srG[,2]),]
	if (args[4] == "Kingdom") {
		cat("\n"); print('Classifying at kingdom'); cat("\n")
	        DF$pKingdom <- round((1 - predict(Kmod, Df, type = "response",allow.new.levels=TRUE)),2)
		srDF$pKingdom <- round((1 - predict(Kmod, srDF, type = "response",allow.new.levels=TRUE)),2)
	} else { DF$pKingdom <- NA; srDF$pKingdom <- NA }

	if (args[4] == "Phylum") {
		cat("\n"); print('Classifying at phylum'); cat("\n")
	        DF$pPhylum <- round((1 - predict(Pmod, DF, type = "response",allow.new.levels=TRUE)),2)
		srDF$pPhylum <- round((1 - predict(Pmod, srDF, type = "response",allow.new.levels=TRUE)),2)
	} else { DF$pPhylum <- NA; srDF$pPhylum <- NA }
	
	if (args[4] == "Class") {
		cat("\n"); print('Classifying at class'); cat("\n")
	        DF$pClass <- round((1 - predict(Cmod, DF, type = "response",allow.new.levels=TRUE)),2)
		srDF$pClass <- round((1 - predict(Cmod, srDF, type = "response",allow.new.levels=TRUE)),2)
	} else { DF$pClass <- NA; srDF$pClass <- NA }
	cat("\n"); print('Classifying at order'); cat("\n")
	DF$pOrder <- round((1 - predict(Omod, DF, type = "response",allow.new.levels=TRUE)),2)
	srDF$pOrder <- round((1 - predict(Omod, srDF, type = "response",allow.new.levels=TRUE)),2)	
	cat("\n"); print('Classifying at family'); cat("\n")
	DF$pFamily <- round((1 - predict(Fmod, DF, type = "response",allow.new.levels=TRUE)),2)
	DF$pFamily <- ifelse(DF$pOrder > 0.85, DF$pFamily, NA)
	srDF$pFamily <- round((1 - predict(Fmod, srDF, type = "response",allow.new.levels=TRUE)),2)
	srDF$pFamily <- ifelse(srDF$pOrder > 0.85, srDF$pFamily, NA)
	cat("\n"); print('Classifying at genus'); cat("\n")
	DF$pGenus <- round((1 - predict(Gmod, DF, type = "response",allow.new.levels=TRUE)),2)
	DF$pGenus <- ifelse(DF$pFamily > 0.9, DF$pGenus, NA)
	srDF$pGenus <- round((1 - predict(srGmod, srDF, type = "response",allow.new.levels=TRUE)),2)
	srDF$pGenus <- ifelse(srDF$pFamily > 0.9, srDF$pGenus, NA)
	cat("\n"); print('Classifying at species'); cat("\n")
	DF$pSpecies <- round((1 - predict(Smod, DF, type = "response",allow.new.levels=TRUE)),2)
	DF$pSpecies <- ifelse(DF$pGenus > 0.95, DF$pSpecies, NA)
	srDF$pSpecies <- round((1 - predict(srSmod, srDF, type = "response",allow.new.levels=TRUE)),2)
	srDF$pSpecies <- ifelse(srDF$pGenus > 0.95, srDF$pSpecies, NA)
	cat("\n"); print('Formatting results'); cat("\n")
	tmpDF <- rbind(srDF,DF)
	NewDF <- tmpDF[,c(1,2,24,3,25,4,26,5,27,6,28,7,29,8,30)] }

write.table(NewDF, args[3], row.names = FALSE, col.names = FALSE, sep=",", quote=FALSE, )
