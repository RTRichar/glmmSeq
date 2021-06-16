#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1]: LogReg formated lineages, DB_directory, lowest reLevel, highest rank 

# load libraries and set some functions
library(MuMIn)
library(glmmTMB)
library(pROC)

# could use function to adjust verbosity
ModSummary <- function(x) {
	Rsq <- r.squaredGLMM(x)
	Smry <- summary(x)
	MyList <- list(Rsq,Smry)
	return(MyList)}

# import data, name columns
CV <- read.csv(args[1], header=FALSE, na.strings=c("","NA"))
colnames(CV) <- c('GI','kReal','pReal','cReal','oReal','fReal','gReal','sReal','kPred','pPred','cPred','oPred','fPred','gPred','sPred','NA','Length','TopID','oScndID','fScndID','gScndID','sScndID','oTpHts','fTpHts','gTpHts','sTpHts')

# Transform if specified
CV$TopID <- CV$TopID^(1/2)

# mark hits with only one hit to the top taxon
CV$sTpHts <- ifelse(CV$sTpHts>1&CV$sTpHts<4,2,CV$sTpHts); CV$sTpHts <- ifelse(CV$sTpHts>3,3,CV$sTpHts); CV$sTpHts <- as.factor(CV$sTpHts)
CV$gTpHts <- ifelse(CV$gTpHts>1&CV$gTpHts<4,2,CV$gTpHts); CV$gTpHts <- ifelse(CV$gTpHts>3,3,CV$gTpHts); CV$gTpHts <- as.factor(CV$gTpHts)
CV$fTpHts <- ifelse(CV$fTpHts>1&CV$fTpHts<4,2,CV$fTpHts); CV$fTpHts <- ifelse(CV$fTpHts>3,3,CV$fTpHts); CV$fTpHts <- as.factor(CV$fTpHts)
CV$oTpHts <- ifelse(CV$oTpHts>1&CV$oTpHts<4,2,CV$oTpHts); CV$oTpHts <- ifelse(CV$oTpHts>3,3,CV$oTpHts); CV$oTpHts <- as.factor(CV$oTpHts)

# if second hit option is RandomEffect, bin ScndID distances
CV$sDist <- (CV$TopID - CV$sScndID)
CV$gDist <- (CV$TopID - CV$gScndID)
CV$fDist <- (CV$TopID - CV$fScndID)
CV$oDist <- (CV$TopID - CV$oScndID)
CV$sDist <- ifelse(CV$sDist<1,1,CV$sDist); CV$sDist <- ifelse(CV$sDist>1,2,CV$sDist); CV$sDist <- as.factor(CV$sDist)
CV$gDist <- ifelse(CV$gDist<1,1,CV$gDist); CV$gDist <- ifelse(CV$gDist>1,2,CV$gDist); CV$gDist <- as.factor(CV$gDist)
CV$fDist <- ifelse(CV$fDist<1,1,CV$fDist); CV$fDist <- ifelse(CV$fDist>1,2,CV$fDist); CV$fDist <- as.factor(CV$fDist)
CV$oDist <- ifelse(CV$oDist<1,1,CV$oDist); CV$oDist <- ifelse(CV$oDist>1,2,CV$oDist); CV$oDist <- as.factor(CV$oDist)

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
	save(Pmod, file = paste(args[2], "/PhylumGLMM.rda", sep='')); print(ModSummary(Pmod)) }
# running class modelling
if (args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
	cat("\n"); print('producing class model'); cat("\n")
	CV <- subset(CV, CV$TopID >= 20^(1/2))
	Cmod <- glmmTMB(c ~ TopID + (1|pReal/cReal), data = CV, family = binomial()) # nest to class
	save(Cmod, file = paste(args[2], '/', "ClassGLMM.rda", sep='')); print(ModSummary(Cmod)) }
# running order modelling
cat("\n"); print('producing order model'); cat("\n")
CV <- subset(CV, CV$TopID >= 60^(1/2))
Omod <- glmmTMB(o ~ TopID + (1|oTpHts) + (1|oDist) + (1|cReal/oReal), data = CV, family = binomial()) # nest to order ### + V17
save(Omod, file = paste(args[2], '/', "OrderGLMM.rda", sep='')); print(ModSummary(Omod))
# running family modelling
cat("\n"); print('producing family model'); cat("\n")
CV <- subset(CV, CV$TopID >= 70^(1/2))
Fmod <- glmmTMB(f ~ TopID + (1|fTpHts) + (1|fDist) + (1|cReal/oReal/fReal), data = CV, family = binomial()) # nest glmm down to family
save(Fmod, file = paste(args[2], '/', "FamilyGLMM.rda", sep='')); print(ModSummary(Fmod))#; print(summary(Fmod))
# running Genus modelling
if (args[3] == "Family") {
	cat("\n"); print('producing genus model w/ family as lowest random intercept term'); cat("\n")
	CV <- subset(CV, CV$TopID >= 80^(1/2))
	Gmod <- glmmTMB(g ~ TopID + (1|gTpHts) + (1|gDist) + (1|cReal/oReal/fReal), data = CV, family = binomial()) 
	#CV$Pred <- predict(Gmod,type=c("response")); roccurve<-roc(CV$g~CV$Pred); print(auc(roccurve))
	save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod)) }
# running species modelling
if (args[3] == "Family") {
	cat("\n"); print('producing species model w/ family as lowest random intercept term'); cat("\n\n")
	CV <- subset(CV, CV$TopID >= 85^(1/2))
	Smod <- glmmTMB(s ~ TopID + (1|sTpHts) + (1|sDist) + (1|cReal/oReal/fReal), data = CV, family = binomial())
	save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod)) }
# default to family but include options for accelerated genus level random effect implementation
if (args[3] == "SpeedGenus") {
        cat("\n"); print('producing genus models w/ speed genus option'); cat("\n")
        CV <- subset(CV, CV$TopID >= 90^(1/2))
        Gmod <- glmmTMB(g ~ TopID + (1|gTpHts) + (1|gDist) + (1|cReal/oReal/fReal), data = CV, family = binomial())
        save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod))
	Genera <- unique(levels(CV$gReal))
	GenusDF <- as.data.frame(matrix(nrow=0,ncol=ncol(CV)))
	colnames(GenusDF) <- colnames(CV)
	print(head(GenusDF))
	for (i in Genera) {
		tmpDF <- CV[CV$gReal==i,]
		NumSpp <- length(unique(levels(droplevels(tmpDF$sReal))))
		NumCases <- nrow(tmpDF)
		if (NumSpp >= 5 & NumCases >= 15) { GenusDF <- rbind(GenusDF,tmpDF) } } # GenusDF is df used to model species rich genera
	cat("\n"); print(paste('Number of time-consuming, species rich genera: ',length(unique(levels(droplevels(GenusDF$gReal)))),sep='')); cat("\n")
	srGmod <- glmmTMB(g ~ TopID + (1|gTpHts) + (1|gDist) + (1|oReal/fReal/gReal), data = GenusDF, family = binomial())
	save(srGmod, file = paste(args[2], '/', "srGenusGLMM.rda", sep='')); print(ModSummary(srGmod))
	cat("\n"); print('producing species models w/ speed genus option'); cat("\n")
	CV <- subset(CV, CV$TopID >= 95^(1/2))
	GenusDF <- subset(GenusDF, GenusDF$TopID >= 95^(1/2))
	Smod <- glmmTMB(s ~ TopID + (1|sTpHts) + (1|sDist) + (1|cReal/oReal/fReal), data = CV, family = binomial())
	save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod))
	srSmod <- glmmTMB(s ~ TopID + (1|sTpHts) + (1|sDist) + (1|oReal/fReal/gReal), data = GenusDF, family = binomial())
	save(srSmod, file = paste(args[2], '/', "srSpeciesGLMM.rda", sep='')); print(ModSummary(srSmod))}
# write species rich genera to file
srG <- unique(levels(droplevels(GenusDF$gReal)))
write.csv(srG, file = paste(args[2], "/srGenera.csv", sep=''))


