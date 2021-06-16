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
load(paste(args[2], "SpeciesGLMM.rda", sep=''))

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
NewDF <- as.data.frame(df[,1]); colnames(NewDF) <- c('Accession')

NewDF$Kingdom <- df$kReal
if (args[4] == "Kingdom") {
	NewDF$pKingdom <- round((1 - predict(Kmod, df, type = "response",allow.new.levels=TRUE)),2)
} else { NewDF$pKingdom <- NA }

NewDF$Phylum <- df$pReal
if (args[4] == "Phylum") {
	NewDF$pPhylum <- round((1 - predict(Pmod, df, type = "response",allow.new.levels=TRUE)),2)
} else { NewDF$pPhylum <- NA }

NewDF$Class <- df$cReal
if (args[4] == "Class") {
	NewDF$pClass <- round((1 - predict(Cmod, df, type = "response",allow.new.levels=TRUE)),2)
} else { NewDF$pClass <- NA }

NewDF$Order <- df$oReal
NewDF$pOrder <- round((1 - predict(Omod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Family <- df$fReal
NewDF$pFamily <- round((1 - predict(Fmod, df, type = "response",allow.new.levels=TRUE)),2)

# if -re == Family
if (args[5] == "Family") {
	NewDF$Genus <- df$gReal
	NewDF$pGenus <- round((1 - predict(Gmod, df, type = "response",allow.new.levels=TRUE)),2)
	NewDF$pGenus <- ifelse(NewDF$pFamily > 0.95, NewDF$pGenus, NA)
	NewDF$Species <- df$sReal
	NewDF$pSpecies <- round((1 - predict(Smod, df, type = "response",allow.new.levels=TRUE)),2)
	NewDF$pSpecies <- ifelse(NewDF$pGenus>0.95, NewDF$pSpecies, NA)}

# if -re == SpeedGenus
if (args[5] == "SpeedGenus") {
	# import genus list
	read.csv()
	# split Q based on occurence in genus list
	# run Qs outside genus list
	# run Qs within genus list
	# recombine
        NewDF$Genus <- df$gReal
        NewDF$pGenus <- round((1 - predict(Gmod, df, type = "response",allow.new.levels=TRUE)),2)
        NewDF$pGenus <- ifelse(NewDF$pFamily > 0.95, NewDF$pGenus, NA)
        NewDF$Species <- df$sReal
        NewDF$pSpecies <- round((1 - predict(Smod, df, type = "response",allow.new.levels=TRUE)),2)
        NewDF$pSpecies <- ifelse(NewDF$pGenus>0.95, NewDF$pSpecies, NA)}

write.table(NewDF, args[3], row.names = FALSE, col.names = FALSE, sep=",", quote=FALSE, )
