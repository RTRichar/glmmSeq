#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1]: Metaxa2 formated tax lineages, args[2] abatase path, args[3] outfile, args[4] highest rank analyzed 

library(MuMIn)
library(glmmTMB)

if (args[4] == "Kingdom") {
	load(paste(args[2], "KingdomGLMM.rda", sep=''))
	load(paste(args[2], "PhylumGLMM.rda", sep=''))
}
if (args[4] == "Phylum") {
        load(paste(args[2], "PhylumGLMM.rda", sep=''))
}
if (args[4] == "Class") {
	load(paste(args[2], "ClassGLMM.rda", sep=''))
}
load(paste(args[2], "OrderGLMM.rda", sep=''))
load(paste(args[2], "FamilyGLMM.rda", sep=''))
load(paste(args[2], "GenusGLMM.rda", sep=''))
load(paste(args[2], "SpeciesGLMM.rda", sep=''))

df <- read.csv(args[1], header=FALSE, na.strings=c("","NA"))
df$V18 <- df$V10
df$V17 <- df$V11

# relationship b/w training LogReg file and tax lineages for undergoing classification
#g,k, p, c, o, f,g,s,i,l
#1,2, 3, 4, 5, 6,7,8,9,10
#     p, c, o, f
#    10,11,12,13

# create dataframe to hold output
NewDF <- as.data.frame(df[,1]); colnames(NewDF) <- c('Accession')

NewDF$Kingdom <- df$V2
if (args[4] == "Kingdom") {
	NewDF$pKingdom <- round((1 - predict(Kmod, df, type = "response",allow.new.levels=TRUE)),2)
} else {
	NewDF$pKingdom <- NA
}

NewDF$Phylum <- df$V3
if (args[4] == "Phylum") {
	NewDF$pPhylum <- round((1 - predict(Pmod, df, type = "response",allow.new.levels=TRUE)),2)
} else {
	NewDF$pPhylum <- NA
}

NewDF$Class <- df$V4
if (args[4] == "Class") {
	NewDF$pClass <- round((1 - predict(Cmod, df, type = "response",allow.new.levels=TRUE)),2)
} else {
	NewDF$pClass <- NA
}

NewDF$Order <- df$V5
NewDF$pOrder <- round((1 - predict(Omod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Family <- df$V6
NewDF$pFamily <- round((1 - predict(Fmod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Genus <- df$V7
NewDF$pGenus <- round((1 - predict(Gmod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Species <- df$V8
NewDF$pSpecies <- round((1 - predict(Smod, df, type = "response",allow.new.levels=TRUE)),2)

write.csv(NewDF, args[3])

# remove NA clauses and just go with assumption that labels are corrent
# make --reLevel specifiable to species genus or family during training
# make --Scope specifiable to Kingdom, Phylum or Class (default)
