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

#df$V18 <- df$V10 #^(1/2)
#df$V17 <- df$V11

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

NewDF$Genus <- df$gReal
NewDF$pGenus <- round((1 - predict(Gmod, df, type = "response",allow.new.levels=TRUE)),2)
NewDF$pGenus <- ifelse(NewDF$pFamily > 0.95, NewDF$pGenus, NA)

NewDF$Species <- df$sReal
NewDF$pSpecies <- round((1 - predict(Smod, df, type = "response",allow.new.levels=TRUE)),2)
NewDF$pSpecies <- ifelse(NewDF$pGenus>0.95, NewDF$pSpecies, NA)

write.table(NewDF, args[3], row.names = FALSE, col.names = FALSE, sep=",", quote=FALSE, )
