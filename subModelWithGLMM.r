#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1]: LogReg formated lineages, DB_directory, lowest reLevel, highest rank 

# load libraries and set some functions
library(MuMIn)
library(glmmTMB)

# could use function to adjust verbosity
ModSummary <- function(x) {
	return(paste('Conditional R-sqr: ', unname(r.squaredGLMM(x)[,2]), sep=''))
	# if 'verbose' == T
		#print(summary(Kmod))
		#print(r.squaredGLMM(Kmod)) 
}

# import data, name columns
CV <- read.csv(args[1], header=FALSE, na.strings=c("","NA"))
colnames(CV) <- c('GI','kReal','pReal','cReal','oReal','fReal','gReal','sReal','kPred','pPred','cPred','oPred','fPred','gPred','sPred','NA','Length','TopID','oScndID','fScndID','gScndID','sScndID','oTpHts','fTpHts','gTpHts','sTpHts')

# Call binomial outcomes
CV$s <- ifelse(as.character(CV$sPred)!=as.character(CV$sReal) & !is.na(CV$sPred) & !is.na(CV$sReal), 1, 0) # change V8, V15, etc to sPred, sReal, etc
CV$g <- ifelse(as.character(CV$gPred)!=as.character(CV$gReal) & !is.na(CV$gPred) & !is.na(CV$gReal), 1, 0)
CV$f <- ifelse(as.character(CV$fPred)!=as.character(CV$fReal) & !is.na(CV$fPred) & !is.na(CV$fReal), 1, 0)
CV$o <- ifelse(as.character(CV$oPred)!=as.character(CV$oReal) & !is.na(CV$oPred) & !is.na(CV$oReal), 1, 0)
if (args[4] == "Class") {
	CV$c <- ifelse(as.character(CV$cPred)!=as.character(CV$cReal) & !is.na(CV$cPred) & !is.na(CV$cReal), 1, 0) }
if (args[4] == "Phylum") {
	CV$p <- ifelse(as.character(CV$pPred)!=as.character(CV$pReal) & !is.na(CV$pPred) & !is.na(CV$pReal), 1, 0) }
if (args[4] == "Kingdom") {
	CV$k <- ifelse(as.character(CV$kPred)!=as.character(CV$kReal) & !is.na(CV$kPred) & !is.na(CV$kReal), 1, 0) }

# calculate distance between 1st and second aligned taxa
CV$gScndID <- sqrt(CV$TopID-CV$gScndID) # could add option for un-transformed distance

# Transform if specified
#CV$TopID <- CV$TopID^(1/2)

# write df used for glmm modeling
write.csv(CV, file=paste(args[2], '/', 'LogReg_rDF.csv', sep=''))

# running kingdom modelling
if (args[4] == "Kingdom") {
	cat("\n"); print('producing kingdom model'); cat("\n")
	CV <- subset(CV, CV$TopID >= 10)
	Kmod <- glmmTMB(k ~ TopID + (1|kReal), data = CV, family = binomial()) 
	save(Kmod, file = paste(args[2], '/', "KingdomGLMM.rda", sep='')); print(ModSummary(Kmod)) }
# running phylum modelling
if (args[4] == "Phylum" | args[4] == "Kingdom") {
	cat("\n"); print('producing phylum model'); cat("\n")
	CV <- subset(CV, CV$TopID >= 10)
	Pmod <- glmmTMB(p ~ TopID + (1|kReal/pReal), data = CV, family = binomial()) 
	save(Pmod, file = paste(args[2], '/', "PhylumGLMM.rda", sep='')); print(ModSummary(Pmod)) }
# running class modelling
if (args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
	cat("\n"); print('producing class model'); cat("\n")
	CV <- subset(CV, CV$TopID >= 20^(1/2))
	Cmod <- glmmTMB(c ~ TopID + (1|pReal/cReal), data = CV, family = binomial()) # nest to class
	save(Cmod, file = paste(args[2], '/', "ClassGLMM.rda", sep='')); print(ModSummary(Cmod)) }
# running order modelling
cat("\n"); print('producing order model'); cat("\n")
CV <- subset(CV, CV$TopID >= 60^(1/2))
Omod <- glmmTMB(o ~ TopID + oScndID + oTpHts + (1|cReal/oReal), data = CV, family = binomial()) # nest to order ### + V17
save(Omod, file = paste(args[2], '/', "OrderGLMM.rda", sep='')); print(ModSummary(Omod))
# running family modelling
cat("\n"); print('producing family model'); cat("\n")
CV <- subset(CV, CV$TopID >= 70^(1/2))
Fmod <- glmmTMB(f ~ TopID + fScndID + fTpHts + (1|cReal/oReal/fReal), data = CV, family = binomial()) # nest glmm down to family
save(Fmod, file = paste(args[2], '/', "FamilyGLMM.rda", sep='')); print(ModSummary(Fmod))
# running Genus modelling
if (args[3] == "Family") {
	cat("\n"); print('producing genus model w/ family as lowest random intercept term'); cat("\n")
	CV <- subset(CV, CV$TopID >= 80^(1/2))
	Gmod <- glmmTMB(g ~ TopID + gScndID + gTpHts + (1|cReal/oReal/fReal), data = CV, family = binomial()) 
	save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod)) }
# running species modelling
if (args[3] == "Family") {
	cat("\n"); print('producing species model w/ family as lowest random intercept term'); cat("\n\n")
	CV <- subset(CV, CV$TopID >= 85^(1/2))
	Smod <- glmmTMB(s ~ TopID + sScndID + sTpHts + (1|cReal/oReal/fReal), data = CV, family = binomial())
	save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod)) }

## default to family but include options for accelerated genus or species level random effect implementations
