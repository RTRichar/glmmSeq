#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1]: Metaxa2 formated tax lineages, dabatase path, outfile

library(MuMIn)
library(glmmTMB)
#load(paste(args[2], "KingdomGLMM.rda", sep=''))
#load(paste(args[2], "PhylumGLMM.rda", sep=''))
load(paste(args[2], "ClassGLMM.rda", sep=''))
load(paste(args[2], "OrderGLMM.rda", sep=''))
load(paste(args[2], "FamilyGLMM.rda", sep=''))
load(paste(args[2], "GenusGLMM.rda", sep=''))
load(paste(args[2], "SpeciesGLMM.rda", sep=''))

#g,k, p, c, o, f,g,s,i,l
#1,2, 3, 4, 5, 6,7,8,9,10
#     p, c, o, f
#    10,11,12,13

df <- read.csv(args[1], header=FALSE, na.strings=c("","NA"))
print(head(df))
df$V18 <- df$V10
df$V17 <- df$V11
df$V10 <- df$V6
df$V11 <- df$V5
df$V12 <- df$V4
df$V13 <- df$V3
print(head(df))

NewDF <- as.data.frame(df[,1]); colnames(NewDF) <- c('Accession')

NewDF$Kingdom <- df$V2
#NewDF$pKingdom <- round((1 - predict(Kmod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Phylum <- df$V3
#NewDF$pPhylum <- round((1 - predict(Pmod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Class <- df$V4
NewDF$pClass <- round((1 - predict(Cmod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Order <- df$V5
NewDF$pOrder <- round((1 - predict(Omod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Family <- df$V6
NewDF$pFamily <- round((1 - predict(Fmod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Genus <- df$V7
NewDF$pGenus <- round((1 - predict(Gmod, df, type = "response",allow.new.levels=TRUE)),2)

NewDF$Species <- df$V8
NewDF$pSpecies <- round((1 - predict(Smod, df, type = "response",allow.new.levels=TRUE)),2)

write.csv(NewDF, args[3])
