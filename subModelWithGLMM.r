#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#str(CTEMPDIR+'/CV.LogReg.csv'), DBDIR, args.reStructure, args.HighestRank, args.ModelStructure, str(args.idCutoffs), str(args.gRichness), str(args.sqrt)

# 1 = str(CTEMPDIR+'/CV.LogReg.csv')
# 2 = DBDIR 
# 3 = args.reStructure ()
# 4 = highest rank ()
# 5 = ModelStructure (FullModel, MidModel, ReducedModel and MinModel)
# 6 = idCuttofs ()
# 7 = gRichness ()
# 8 = sqrt ("True/False" not T/F or FTRUE/FALSE)

# load libraries and set some functions
library(MuMIn)
library(glmmTMB)
library(pROC)

# could use function to adjust verbosity
ModSummary <- function(x,r) {
	Rsq <- r.squaredGLMM(x)
	#Smry <- summary(x)
	CV[,'Pred'] <- predict(x,newdata=CV,type=c("response"), allow.new.levels = T)
	roccurve<-roc(CV[,r]~CV[,'Pred'])
	MyList <- list(Rsq, auc(roccurve)) #,Smry)
	return(MyList)}

# import data, name columns
CV <- read.csv(args[1], header=FALSE, na.strings=c("","NA"))
colnames(CV) <- c('GI','kReal','pReal','cReal','oReal','fReal','gReal','sReal','kPred','pPred','cPred','oPred','fPred','gPred','sPred','NA','Length','TopID','oScndID','fScndID','gScndID','sScndID','oTpHts','fTpHts','gTpHts','sTpHts')

# Transform if specified
for (i in 1:8) { print(args[i]) }

if (args[8] == "True") {
	x <- c('TopID','oScndID','fScndID','gScndID','sScndID')
	for (i in x) { CV[,i] <- sqrt(CV[,i]) } 
	idCutoffs <- sqrt(as.numeric(as.character(unlist(strsplit(args[6], ",")))))
	print(idCutoffs)
} else { idCutoffs <- as.numeric(as.character(unlist(strsplit(args[6], ",")))) }

# mark hits with only one hit to the top taxon
#### Make function for doing this and split to more categories
CV$sTpHts <- ifelse(CV$sTpHts>1&CV$sTpHts<4,2,CV$sTpHts); CV$sTpHts <- ifelse(CV$sTpHts>3,3,CV$sTpHts); CV$sTpHts <- as.factor(CV$sTpHts)
CV$gTpHts <- ifelse(CV$gTpHts>1&CV$gTpHts<4,2,CV$gTpHts); CV$gTpHts <- ifelse(CV$gTpHts>3,3,CV$gTpHts); CV$gTpHts <- as.factor(CV$gTpHts)
CV$fTpHts <- ifelse(CV$fTpHts>1&CV$fTpHts<4,2,CV$fTpHts); CV$fTpHts <- ifelse(CV$fTpHts>3,3,CV$fTpHts); CV$fTpHts <- as.factor(CV$fTpHts)
CV$oTpHts <- ifelse(CV$oTpHts>1&CV$oTpHts<4,2,CV$oTpHts); CV$oTpHts <- ifelse(CV$oTpHts>3,3,CV$oTpHts); CV$oTpHts <- as.factor(CV$oTpHts)

# if second hit option is RandomEffect, bin ScndID distances
#### Make function for doing this and split to more categories
CV$sDist <- (CV$TopID - CV$sScndID)
CV$gDist <- (CV$TopID - CV$gScndID)
CV$fDist <- (CV$TopID - CV$fScndID)
CV$oDist <- (CV$TopID - CV$oScndID)
CV$sDist <- ifelse(CV$sDist<1,1,CV$sDist); CV$sDist <- ifelse(CV$sDist>1,2,CV$sDist); CV$sDist <- as.factor(CV$sDist)
CV$gDist <- ifelse(CV$gDist<1,1,CV$gDist); CV$gDist <- ifelse(CV$gDist>1,2,CV$gDist); CV$gDist <- as.factor(CV$gDist)
CV$fDist <- ifelse(CV$fDist<1,1,CV$fDist); CV$fDist <- ifelse(CV$fDist>1,2,CV$fDist); CV$fDist <- as.factor(CV$fDist)
CV$oDist <- ifelse(CV$oDist<1,1,CV$oDist); CV$oDist <- ifelse(CV$oDist>1,2,CV$oDist); CV$oDist <- as.factor(CV$oDist)

# Call binomial outcomes
CV$s <- ifelse(as.character(CV$sPred)!=as.character(CV$sReal) & !is.na(CV$sPred) & !is.na(CV$sReal), 1, 0) 
CV$g <- ifelse(as.character(CV$gPred)!=as.character(CV$gReal) & !is.na(CV$gPred) & !is.na(CV$gReal), 1, 0)
CV$f <- ifelse(as.character(CV$fPred)!=as.character(CV$fReal) & !is.na(CV$fPred) & !is.na(CV$fReal), 1, 0)
CV$o <- ifelse(as.character(CV$oPred)!=as.character(CV$oReal) & !is.na(CV$oPred) & !is.na(CV$oReal), 1, 0)
if (args[4] == "Class") {
	CV$c <- ifelse(as.character(CV$cPred)!=as.character(CV$cReal) & !is.na(CV$cPred) & !is.na(CV$cReal), 1, 0) }
if (args[4] == "Phylum") {
	CV$p <- ifelse(as.character(CV$pPred)!=as.character(CV$pReal) & !is.na(CV$pPred) & !is.na(CV$pReal), 1, 0) }
if (args[4] == "Kingdom") {
	CV$k <- ifelse(as.character(CV$kPred)!=as.character(CV$kReal) & !is.na(CV$kPred) & !is.na(CV$kReal), 1, 0) }

#print(head(CV))
#print(class(CV$oTpHts))
#print(range(CV$oTpHts, na.rm = T))
#print(range(CV$fTpHts, na.rm = T))

# write df used for glmm modeling
write.csv(CV, file=paste(args[2], '/', 'LogReg_rDF.csv', sep=''))

# FullModel, MidModel, ReducedModel and MinModel
#####################################################################################################################################################
if (args[5] == "FullModel") {
        if (args[4] == "Kingdom") {
                cat("\n"); print('producing kingdom model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[1])
                Kmod <- glmmTMB(k ~ TopID + (TopID|kPred), data = CV, family = binomial())
                save(Kmod, file = paste(args[2], '/', "KingdomGLMM.rda", sep='')) }
        if (args[4] == "Phylum" | args[4] == "Kingdom") {
                cat("\n"); print('producing phylum model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[2])
                Pmod <- glmmTMB(p ~ TopID + (TopID|kPred/pPred), data = CV, family = binomial())
                save(Pmod, file = paste(args[2], "/PhylumGLMM.rda", sep='')) }
        if (args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
                cat("\n"); print('producing class model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[3])
                Cmod <- glmmTMB(c ~ TopID + (TopID|pPred/cPred), data = CV, family = binomial())
                save(Cmod, file = paste(args[2], '/', "ClassGLMM.rda", sep='')) }
        if (args[4] == "Order" | args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
                cat("\n"); print('producing order model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[4])
                Omod <- glmmTMB(o ~ TopID + (TopID|cPred/oPred), data = CV, family = binomial()) # for some reason, we cant access oTpHts
                save(Omod, file = paste(args[2], '/', "OrderGLMM.rda", sep='')) }
        cat("\n"); print('producing family model'); cat("\n")
        CV <- subset(CV, CV$TopID >= idCutoffs[5])
        Fmod <- glmmTMB(f ~ TopID + (TopID|fTpHts) + (TopID|fDist) + (TopID|cPred/oPred/fPred), data = CV, family = binomial()) # nest glmm down to family
        save(Fmod, file = paste(args[2], '/', "FamilyGLMM.rda", sep='')); print(ModSummary(Fmod,'f'))#; print(summary(Fmod))
        # running Genus and species modelling, -re == Family
        if (args[3] == "Family") {
                cat("\n"); print('producing genus model w/ family as lowest random intercept term'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[6])
                Gmod <- glmmTMB(g ~ TopID + (TopID|gTpHts) + (TopID|gDist) + (TopID|cPred/oPred/fPred), data = CV, family = binomial())
                save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
                cat("\n"); print('producing species model w/ family as lowest random intercept term'); cat("\n\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[7])
                Smod <- glmmTMB(s ~ TopID + (TopID|sTpHts) + (TopID|sDist) + (TopID|cPred/oPred/fPred), data = CV, family = binomial())
                save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s')) }
        # running Genus and species modelling, -re == Family
        if (args[3] == "Genus") {
                cat("\n"); print('producing genus model w/ genus as lowest random intercept term'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[6])
                Gmod <- glmmTMB(g ~ TopID + (TopID|gTpHts) + (TopID|gDist) + (TopID|oPred/fPred/gPred), data = CV, family = binomial())
                save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
                cat("\n"); print('producing species model w/ genus as lowest random intercept term'); cat("\n\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[7])
                Smod <- glmmTMB(s ~ TopID + (TopID|sTpHts) + (TopID|sDist) + (TopID|oPred/fPred/gPred), data = CV, family = binomial())
                save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s')) }
        # default to family but include options for accelerated genus level random effect implementation
        if (args[3] == "SpeedGenus") {
                cat("\n"); print('producing genus models w/ speed genus option'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[6])
                Gmod <- glmmTMB(g ~ TopID + (TopID|gTpHts) + (TopID|gDist) + (TopID|cPred/oPred/fPred), data = CV, family = binomial())
                save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
                Genera <- unique(levels(CV$gReal))
                GenusDF <- as.data.frame(matrix(nrow=0,ncol=ncol(CV)))
                colnames(GenusDF) <- colnames(CV)
                print(head(GenusDF))
                for (i in Genera) {
                        tmpDF <- CV[CV$gReal==i,]
                        NumSpp <- length(unique(levels(droplevels(tmpDF$sReal))))
                        NumCases <- nrow(tmpDF)
                        if (NumSpp >= 15 & NumCases >= 30) { GenusDF <- rbind(GenusDF,tmpDF) } } # GenusDF is df used to model species rich genera
                cat("\n"); print(paste('Number of time-consuming, species rich genera: ',length(unique(levels(droplevels(GenusDF$gReal)))),sep='')); cat("\n")
                srGmod <- glmmTMB(g ~ TopID + (TopID|gTpHts) + (TopID|gDist) + (TopID|oPred/fPred/gPred), data = GenusDF, family = binomial())
                save(srGmod, file = paste(args[2], '/', "srGenusGLMM.rda", sep='')); print(ModSummary(srGmod,'g'))
                cat("\n"); print('producing species models w/ speed genus option'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[7])
                GenusDF <- subset(GenusDF, GenusDF$TopID >= idCutoffs[7])
                Smod <- glmmTMB(s ~ TopID + (TopID|sTpHts) + (TopID|sDist) + (TopID|cPred/oPred/fPred), data = CV, family = binomial())
                save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s'))
                srSmod <- glmmTMB(s ~ TopID + (TopID|sTpHts) + (TopID|sDist) + (TopID|oPred/fPred/gPred), data = GenusDF, family = binomial())
                save(srSmod, file = paste(args[2], '/', "srSpeciesGLMM.rda", sep='')); print(ModSummary(srSmod,'s'))
                srG <- unique(levels(droplevels(GenusDF$gReal)))
                write.csv(srG, file = paste(args[2], "/srGenera.csv", sep=''))} # write species rich genera to file
}
###############################################################################################################################################
if (args[5] == "MidModel") {
        if (args[4] == "Kingdom") {
                cat("\n"); print('producing kingdom model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[1])
                Kmod <- glmmTMB(k ~ TopID + (1|kPred), data = CV, family = binomial())
                save(Kmod, file = paste(args[2], '/', "KingdomGLMM.rda", sep='')) }
        if (args[4] == "Phylum" | args[4] == "Kingdom") {
                cat("\n"); print('producing phylum model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[2])
                Pmod <- glmmTMB(p ~ TopID + (1|kPred/pPred), data = CV, family = binomial())
                save(Pmod, file = paste(args[2], "/PhylumGLMM.rda", sep='')) }
        if (args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
                cat("\n"); print('producing class model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[3])
                Cmod <- glmmTMB(c ~ TopID + (1|pPred/cPred), data = CV, family = binomial())
                save(Cmod, file = paste(args[2], '/', "ClassGLMM.rda", sep='')) }
        if (args[4] == "Order" | args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
                cat("\n"); print('producing order model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[4])
                Omod <- glmmTMB(o ~ TopID + (1|cPred/oPred), data = CV, family = binomial()) # for some reason, we cant access oTpHts
                save(Omod, file = paste(args[2], '/', "OrderGLMM.rda", sep='')) }
        cat("\n"); print('producing family model'); cat("\n")
        CV <- subset(CV, CV$TopID >= idCutoffs[5])
        Fmod <- glmmTMB(f ~ TopID + (TopID|fTpHts) + (TopID|fDist) + (1|cPred/oPred/fPred), data = CV, family = binomial()) # nest glmm down to family
        save(Fmod, file = paste(args[2], '/', "FamilyGLMM.rda", sep='')); print(ModSummary(Fmod,'f'))#; print(summary(Fmod))
        # running Genus and species modelling, -re == Family
        if (args[3] == "Family") {
                cat("\n"); print('producing genus model w/ family as lowest random intercept term'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[6])
                Gmod <- glmmTMB(g ~ TopID + (TopID|gTpHts) + (TopID|gDist) + (1|cPred/oPred/fPred), data = CV, family = binomial())
                save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
                cat("\n"); print('producing species model w/ family as lowest random intercept term'); cat("\n\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[7])
                Smod <- glmmTMB(s ~ TopID + (TopID|sTpHts) + (TopID|sDist) + (1|cPred/oPred/fPred), data = CV, family = binomial())
                save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s')) }
        # running Genus and species modelling, -re == Family
        if (args[3] == "Genus") {
                cat("\n"); print('producing genus model w/ genus as lowest random intercept term'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[6])
                Gmod <- glmmTMB(g ~ TopID + (TopID|gTpHts) + (TopID|gDist) + (1|oPred/fPred/gPred), data = CV, family = binomial())
                save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
                cat("\n"); print('producing species model w/ genus as lowest random intercept term'); cat("\n\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[7])
                Smod <- glmmTMB(s ~ TopID + (TopID|sTpHts) + (TopID|sDist) + (1|oPred/fPred/gPred), data = CV, family = binomial())
                save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s')) }
        # default to family but include options for accelerated genus level random effect implementation
        if (args[3] == "SpeedGenus") {
                cat("\n"); print('producing genus models w/ speed genus option'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[6])
                Gmod <- glmmTMB(g ~ TopID + (TopID|gTpHts) + (TopID|gDist) + (1|cPred/oPred/fPred), data = CV, family = binomial())
                save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
                Genera <- unique(levels(CV$gReal))
                GenusDF <- as.data.frame(matrix(nrow=0,ncol=ncol(CV)))
                colnames(GenusDF) <- colnames(CV)
                print(head(GenusDF))
                for (i in Genera) {
                        tmpDF <- CV[CV$gReal==i,]
                        NumSpp <- length(unique(levels(droplevels(tmpDF$sReal))))
                        NumCases <- nrow(tmpDF)
                        if (NumSpp >= 15 & NumCases >= 30) { GenusDF <- rbind(GenusDF,tmpDF) } } # GenusDF is df used to model species rich genera
                cat("\n"); print(paste('Number of time-consuming, species rich genera: ',length(unique(levels(droplevels(GenusDF$gReal)))),sep='')); cat("\n")
                srGmod <- glmmTMB(g ~ TopID + (TopID|gTpHts) + (TopID|gDist) + (1|oPred/fPred/gPred), data = GenusDF, family = binomial())
                save(srGmod, file = paste(args[2], '/', "srGenusGLMM.rda", sep='')); print(ModSummary(srGmod,'g'))
                cat("\n"); print('producing species models w/ speed genus option'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[7])
                GenusDF <- subset(GenusDF, GenusDF$TopID >= idCutoffs[7])
                Smod <- glmmTMB(s ~ TopID + (TopID|sTpHts) + (TopID|sDist) + (1|cPred/oPred/fPred), data = CV, family = binomial())
                save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s'))
                srSmod <- glmmTMB(s ~ TopID + (TopID|sTpHts) + (TopID|sDist) + (1|oPred/fPred/gPred), data = GenusDF, family = binomial())
                save(srSmod, file = paste(args[2], '/', "srSpeciesGLMM.rda", sep='')); print(ModSummary(srSmod,'s'))
                srG <- unique(levels(droplevels(GenusDF$gReal)))
                write.csv(srG, file = paste(args[2], "/srGenera.csv", sep=''))} # write species rich genera to file
}
#####################################################################################################################################################
if (args[5] == "ReducedModel") {
        if (args[4] == "Kingdom") {
                cat("\n"); print('producing kingdom model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[1])
                Kmod <- glmmTMB(k ~ TopID + (TopID|kPred), data = CV, family = binomial())
                save(Kmod, file = paste(args[2], '/', "KingdomGLMM.rda", sep='')) }
        if (args[4] == "Phylum" | args[4] == "Kingdom") {
                cat("\n"); print('producing phylum model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[2])
                Pmod <- glmmTMB(p ~ TopID + (TopID|kPred/pPred), data = CV, family = binomial())
                save(Pmod, file = paste(args[2], "/PhylumGLMM.rda", sep='')) }
        if (args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
                cat("\n"); print('producing class model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[3])
                Cmod <- glmmTMB(c ~ TopID + (TopID|pPred/cPred), data = CV, family = binomial())
                save(Cmod, file = paste(args[2], '/', "ClassGLMM.rda", sep='')) }
        if (args[4] == "Order" | args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
                cat("\n"); print('producing order model'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[4])
                Omod <- glmmTMB(o ~ TopID + (TopID|cPred/oPred), data = CV, family = binomial()) # for some reason, we cant access oTpHts
                save(Omod, file = paste(args[2], '/', "OrderGLMM.rda", sep='')) }
        cat("\n"); print('producing family model'); cat("\n")
        CV <- subset(CV, CV$TopID >= idCutoffs[5])
        Fmod <- glmmTMB(f ~ TopID + (TopID|cPred/oPred/fPred), data = CV, family = binomial()) # nest glmm down to family
        save(Fmod, file = paste(args[2], '/', "FamilyGLMM.rda", sep='')); print(ModSummary(Fmod,'f'))#; print(summary(Fmod))
        # running Genus and species modelling, -re == Family
        if (args[3] == "Family") {
                cat("\n"); print('producing genus model w/ family as lowest random intercept term'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[6])
                Gmod <- glmmTMB(g ~ TopID + (TopID|cPred/oPred/fPred), data = CV, family = binomial())
                save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
                cat("\n"); print('producing species model w/ family as lowest random intercept term'); cat("\n\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[7])
                Smod <- glmmTMB(s ~ TopID + (TopID|cPred/oPred/fPred), data = CV, family = binomial())
                save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s')) }
        # running Genus and species modelling, -re == Family
        if (args[3] == "Genus") {
                cat("\n"); print('producing genus model w/ genus as lowest random intercept term'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[6])
                Gmod <- glmmTMB(g ~ TopID + (TopID|oPred/fPred/gPred), data = CV, family = binomial())
                save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
                cat("\n"); print('producing species model w/ genus as lowest random intercept term'); cat("\n\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[7])
                Smod <- glmmTMB(s ~ TopID + (TopID|oPred/fPred/gPred), data = CV, family = binomial())
                save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s')) }
        # default to family but include options for accelerated genus level random effect implementation
        if (args[3] == "SpeedGenus") {
                cat("\n"); print('producing genus models w/ speed genus option'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[6])
                Gmod <- glmmTMB(g ~ TopID + (TopID|cPred/oPred/fPred), data = CV, family = binomial())
                save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
                Genera <- unique(levels(CV$gReal))
                GenusDF <- as.data.frame(matrix(nrow=0,ncol=ncol(CV)))
                colnames(GenusDF) <- colnames(CV)
                print(head(GenusDF))
                for (i in Genera) {
                        tmpDF <- CV[CV$gReal==i,]
                        NumSpp <- length(unique(levels(droplevels(tmpDF$sReal))))
                        NumCases <- nrow(tmpDF)
                        if (NumSpp >= 15 & NumCases >= 30) { GenusDF <- rbind(GenusDF,tmpDF) } } # GenusDF is df used to model species rich genera
                cat("\n"); print(paste('Number of time-consuming, species rich genera: ',length(unique(levels(droplevels(GenusDF$gReal)))),sep='')); cat("\n")
                srGmod <- glmmTMB(g ~ TopID + (TopID|oPred/fPred/gPred), data = GenusDF, family = binomial())
                save(srGmod, file = paste(args[2], '/', "srGenusGLMM.rda", sep='')); print(ModSummary(srGmod,'g'))
                cat("\n"); print('producing species models w/ speed genus option'); cat("\n")
                CV <- subset(CV, CV$TopID >= idCutoffs[7])
                GenusDF <- subset(GenusDF, GenusDF$TopID >= idCutoffs[7])
                Smod <- glmmTMB(s ~ TopID + (TopID|cPred/oPred/fPred), data = CV, family = binomial())
                save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s'))
                srSmod <- glmmTMB(s ~ TopID + (TopID|oPred/fPred/gPred), data = GenusDF, family = binomial())
                save(srSmod, file = paste(args[2], '/', "srSpeciesGLMM.rda", sep='')); print(ModSummary(srSmod,'s'))
                srG <- unique(levels(droplevels(GenusDF$gReal)))
                write.csv(srG, file = paste(args[2], "/srGenera.csv", sep=''))} # write species rich genera to file
}
#####################################################################################################################################################
if (args[5] == "MinModel") {
	if (args[4] == "Kingdom") {
		cat("\n"); print('producing kingdom model'); cat("\n")
		CV <- subset(CV, CV$TopID >= idCutoffs[1]) 
		Kmod <- glmmTMB(k ~ TopID + (1|kPred), data = CV, family = binomial()) 
		save(Kmod, file = paste(args[2], '/', "KingdomGLMM.rda", sep='')) }
	if (args[4] == "Phylum" | args[4] == "Kingdom") {
		cat("\n"); print('producing phylum model'); cat("\n")
		CV <- subset(CV, CV$TopID >= idCutoffs[2])
		Pmod <- glmmTMB(p ~ TopID + (1|kPred/pPred), data = CV, family = binomial()) 
		save(Pmod, file = paste(args[2], "/PhylumGLMM.rda", sep='')) }
	if (args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
		cat("\n"); print('producing class model'); cat("\n")
		CV <- subset(CV, CV$TopID >= idCutoffs[3])
		Cmod <- glmmTMB(c ~ TopID + (1|pPred/cPred), data = CV, family = binomial()) 
		save(Cmod, file = paste(args[2], '/', "ClassGLMM.rda", sep='')) }
	if (args[4] == "Order" | args[4] == "Class" | args[4] == "Phylum" | args[4] == "Kingdom"){
		cat("\n"); print('producing order model'); cat("\n")
		CV <- subset(CV, CV$TopID >= idCutoffs[4])
		Omod <- glmmTMB(o ~ TopID + (1|cPred/oPred), data = CV, family = binomial()) # for some reason, we cant access oTpHts
		save(Omod, file = paste(args[2], '/', "OrderGLMM.rda", sep='')) }
	cat("\n"); print('producing family model'); cat("\n")
	CV <- subset(CV, CV$TopID >= idCutoffs[5])
	Fmod <- glmmTMB(f ~ TopID + (1|cPred/oPred/fPred), data = CV, family = binomial()) # nest glmm down to family
	save(Fmod, file = paste(args[2], '/', "FamilyGLMM.rda", sep='')); print(ModSummary(Fmod,'f'))#; print(summary(Fmod))
	# running Genus and species modelling, -re == Family
	if (args[3] == "Family") {
		cat("\n"); print('producing genus model w/ family as lowest random intercept term'); cat("\n")
		CV <- subset(CV, CV$TopID >= idCutoffs[6])
		Gmod <- glmmTMB(g ~ TopID + (1|cPred/oPred/fPred), data = CV, family = binomial()) 
		save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
		cat("\n"); print('producing species model w/ family as lowest random intercept term'); cat("\n\n")
		CV <- subset(CV, CV$TopID >= idCutoffs[7])
		Smod <- glmmTMB(s ~ TopID + (1|cPred/oPred/fPred), data = CV, family = binomial())
		save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s')) }
	# running Genus and species modelling, -re == Family
	if (args[3] == "Genus") {
		cat("\n"); print('producing genus model w/ genus as lowest random intercept term'); cat("\n")
		CV <- subset(CV, CV$TopID >= idCutoffs[6])
		Gmod <- glmmTMB(g ~ TopID + (1|oPred/fPred/gPred), data = CV, family = binomial())
		save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
		cat("\n"); print('producing species model w/ genus as lowest random intercept term'); cat("\n\n")
		CV <- subset(CV, CV$TopID >= idCutoffs[7])
		Smod <- glmmTMB(s ~ TopID + (1|oPred/fPred/gPred), data = CV, family = binomial())
		save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s')) }
	# default to family but include options for accelerated genus level random effect implementation
	if (args[3] == "SpeedGenus") {
	        cat("\n"); print('producing genus models w/ speed genus option'); cat("\n")
	        CV <- subset(CV, CV$TopID >= idCutoffs[6])
	        Gmod <- glmmTMB(g ~ TopID + (1|cPred/oPred/fPred), data = CV, family = binomial())
	        save(Gmod, file = paste(args[2], '/', "GenusGLMM.rda", sep='')); print(ModSummary(Gmod,'g'))
		Genera <- unique(levels(CV$gReal))
		GenusDF <- as.data.frame(matrix(nrow=0,ncol=ncol(CV)))
		colnames(GenusDF) <- colnames(CV)
		print(head(GenusDF))
		for (i in Genera) {
			tmpDF <- CV[CV$gReal==i,]
			NumSpp <- length(unique(levels(droplevels(tmpDF$sReal))))
			NumCases <- nrow(tmpDF)
			if (NumSpp >= 15 & NumCases >= 30) { GenusDF <- rbind(GenusDF,tmpDF) } } # GenusDF is df used to model species rich genera
		cat("\n"); print(paste('Number of time-consuming, species rich genera: ',length(unique(levels(droplevels(GenusDF$gReal)))),sep='')); cat("\n")
		srGmod <- glmmTMB(g ~ TopID + (1|oPred/fPred/gPred), data = GenusDF, family = binomial())
		save(srGmod, file = paste(args[2], '/', "srGenusGLMM.rda", sep='')); print(ModSummary(srGmod,'g'))
		cat("\n"); print('producing species models w/ speed genus option'); cat("\n")
		CV <- subset(CV, CV$TopID >= idCutoffs[7])
		GenusDF <- subset(GenusDF, GenusDF$TopID >= idCutoffs[7])
		Smod <- glmmTMB(s ~ TopID + (1|cPred/oPred/fPred), data = CV, family = binomial())
		save(Smod, file = paste(args[2], '/', "SpeciesGLMM.rda", sep='')); print(ModSummary(Smod,'s'))
		srSmod <- glmmTMB(s ~ TopID + (1|oPred/fPred/gPred), data = GenusDF, family = binomial())
		save(srSmod, file = paste(args[2], '/', "srSpeciesGLMM.rda", sep='')); print(ModSummary(srSmod,'s'))
		srG <- unique(levels(droplevels(GenusDF$gReal)))
		write.csv(srG, file = paste(args[2], "/srGenera.csv", sep=''))} # write species rich genera to file
}
