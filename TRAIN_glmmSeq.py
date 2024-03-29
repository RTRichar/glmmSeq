#!/usr/bin/env python

import subprocess, argparse, sys, time, os

parser = argparse.ArgumentParser(description="")
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# required
required.add_argument('-i', '--InputFasta', required = True, help = "\n Input file of reference sequences used to train the classifier. Sequence headers should include a unique sequence identifier, NCBI accessions are recomended, preceded only by '>'. Sequences should be contiguous\n")
required.add_argument('-it', '--InputTax', required = True, help = "\n Input taxonomic lineages (a tab-delimited file with the first column containing unique sequence identifiers compatible with the files provided under '-i'). Normally, this would be in Metaxa2-compatible format.\n")
required.add_argument('-odb', '--OutputDB', required = True, help = "\n Name of database to be produced by TRAIN_glmmSeq. This is the name of the directory that will hold the R GLMMs, fasta and taxonomy files of the resulting database. It will be written in /DBs/, within the directory that holds the glmmSeq executables.\n")
#optional
optional.add_argument('-t', '--Threads', required = False, default = 1, help = "\n Number of processors (only affects vsearch alignment speed)\n")
optional.add_argument('--SaveTemp', default=False, type=lambda x: (str(x).lower() == 'true'), help = "\n Option for saving intermediate files produced during training (e.g. fasta k-fold partitions and csv formatted cross validation resulst used for glmm fitting). Must be 'True/False,' not 'TRUE/FALSE' or 'T/F.'\n")
optional.add_argument('-id', '--idCutoffs', required = False, default = '30,40,50,60,70,80,85', help = "\n Option for specifying the minimum percent identity of Vsearch alignment matches to be used for glmm fitting. This is adjustable at each taxonomic rank and consists of a comma-separated list of 7 numbers representing the threshold used from kingdom to species (e.g the default of 30,40,50,60,70,80,85 means alignments of < 30 percent ID will not be used during glmm modelling at the kigdom rank) \n")
optional.add_argument('-k', '--kFolds', required = False, default = '3,12', help = "\n Number of k-folds partitions to use during reference cross-alignment and glmm modelling")
#### change -pcv  to -k 
optional.add_argument('-hr', '--HighestRank', required = False, type=str, default = 'Order', help = "\n Specifies the lowest resolution rank to be classified (Kingdom, Phylum, Class, or Order). Default = Class")
optional.add_argument('-re', '--reStructure', required = False, type=str, default = 'Family', help = "\n Specifies the random effect strucutre to used during modelling. Options include 'Family', 'Genus' and 'SpeedGenus'. 'SpeedGenus' splits cases according to whether they are in species-rich genera. For cases within species rich genera, genus and species models include (1|Order/Family/Genus) as random effect. For non-species-rich genera, random effect term is (1|Class/Order/Family). Default = SpeedGenus")
optional.add_argument('-m', '--ModelStructure', required = False, type=str, default = 'FullModel', help = "\n Provides some options of varying model structure with respect to both random effects and fixed effects. Options include FullModel, MidModel, ReducedModel and MinModel. FullModel: Includes categorically binned 'mumber of top hits' and 'distance between 1st and 2nd top taxon' as random intercept and slope effects with subject taxonomy as random intercept and slope. MidModel: Includes categorically binned 'mumber of top hits' and 'distance between 1st and 2nd top taxon' as random intercept and slope effects with subject taxonomy as random intercept. ReducedModel: Includes subject taxonomy as random intercept and slope. MinModel: Includes subject taxonomy as random intercept.")
optional.add_argument('-g', '--gRichness', required = False, default = 20, help = "\n Minimum number of species/genera for a genera to qualify as species-rich")
optional.add_argument('-s', '--sqrt', required = False, type=str, default = 'True', help = "\n Option for square-root transforming alignment score percent identities. Should be 'True/False,' not 'TRUE/FALSE' or 'T/F.' Default = True")
args = parser.parse_args()

# Test if SaveTemp 'True' 'False,' not T/F or TRUE/FALSE
# Test if -re Family/SpeedGenus

# create temp directory
CTEMPDIR = str('TRAIN_tmp_' + '_'.join(time.ctime(time.time()).replace(':','_').split()[1:]))
subprocess.call(['mkdir', CTEMPDIR])

# create directory to store database
DBDIR = str(os.path.abspath(os.path.dirname(sys.argv[0]))+'/DBs/'+args.OutputDB)
subprocess.call(['mkdir', DBDIR])

# set glmmSeq location
glmmSeqDirectory = os.path.abspath(os.path.dirname(sys.argv[0]))

# make file to hold info on highest rank to analyze to and lowest rank for RE specification
InfoFile = str(DBDIR + '/' + 'InfoFile.txt')
with open(InfoFile, 'w') as File:
	File.write(args.reStructure+';'+args.HighestRank+';'+args.ModelStructure+';'+str(args.idCutoffs)+';'+str(args.gRichness)+';'+str(args.sqrt))

kFolds = str(args.kFolds).split(',')
for k in kFolds: # for i in k-fold series
	# Split into test train
	sys.stderr.write('\n'+time.ctime(time.time())+': Partitioning k-fold testing/training sets -- k-folds = '+str(k)+'\n')
	subprocess.call(['GetTestTrain.py', str(args.InputFasta), str(CTEMPDIR+'/'+str(k)+'Fold'), str(k)])
	# Run Vsearch algnmnt
	sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Running vsearch alignments ###\n\n')
	for i in range(0,int(k)): # get files for inferring percent ID of second-to-top hit taxon
		subprocess.call(['vsearch', '--usearch_global', str(CTEMPDIR+'/'+str(k)+'Fold'+str(i)+'_Test.fasta'), '--db', \
		str(CTEMPDIR+'/'+str(k)+'Fold'+str(i)+'_Train.fasta'), '--id', '0.6', '--maxaccepts', '100', '--maxrejects', '50', '--maxhits', '50', \
		'--gapopen', '0TE', '--gapext', '0TE', '--userout', str(CTEMPDIR+'/'+str(k)+'Fold'+str(i)+'_CV_2nd.vsrch.txt'), '--userfields', \
		'query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov', '--query_cov', '0.95', '--threads', str(args.Threads)])
	# combine vsrch outs into single file
	subprocess.call(['CombineCVs.py', str(CTEMPDIR), str(k), '_CV_2nd.vsrch.txt']) # need to change inputs to explicit input output names
	# get vsearch single hit file from 50 hit file 
	subprocess.call(['GetVsrchTopHit.py', str(CTEMPDIR+'/'+str(k)+'Fold_CV_2nd.vsrch.txt'), str(CTEMPDIR+'/'+str(k)+'Fold_CV.vsrch.txt')])
	# Get LogReg file
	sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Formatting vsearch outputs for GLMM modelling ###\n')
	subprocess.call(['VsearchToMetaxa2.py', '-v',str(CTEMPDIR+'/'+str(k)+'Fold_CV.vsrch.txt'),'-t',args.InputTax,'-o', str(CTEMPDIR+'/'+str(k)+'Fold_CV.mtxa.tax')])
	# 
	subprocess.call(['Get2ndHitTaxID_TRAIN.py', args.InputTax, str(CTEMPDIR+'/'+str(k)+'Fold_CV_2nd.vsrch.txt'), str(CTEMPDIR+'/'+str(k)+'Fold_CV_2nd.vsrch.csv')])
	#
	subprocess.call(['CurateForLogReg.py',str(CTEMPDIR+'/'+str(k)+'Fold_CV.mtxa.tax'),args.InputTax,str(CTEMPDIR+'/'+str(k)+'Fold_CV_2nd.vsrch.csv'),str(CTEMPDIR+'/'+str(k)+'Fold_CV.LogReg.csv'),str(k)])

# Combine LogReg files from each kFold used
subprocess.call(['CombineLogRegs.py',str(CTEMPDIR),str(args.kFolds)])

# combine all x LogReg files
# adjust subModel file with rando effect for k-fold size
# if k-fold series greater than 1, combine into one file, No headers

# R file to train on LogReg (must save ModGLMM:
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Modelling data with binomial GLMMs ###\n\n')
subprocess.call(['subModelWithGLMM.r', str(CTEMPDIR+'/FinalLogReg.csv'), DBDIR, args.reStructure, args.HighestRank, \
	args.ModelStructure, str(args.idCutoffs), str(args.gRichness), str(args.sqrt)])

# save consensus filtered fasta and tax to db directory under 'DB.fa' and 'DB.tax'
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Writing database files and cleaning up tmp files ###\n\n')
subprocess.call(['TaxFastaConsensus.py', '-it', str(args.InputTax), '-if', str(args.InputFasta), '-ot', str(DBDIR+'/DB.tax'), '-of', str(DBDIR+'/DB.fa')])

# Clean up tmp
if bool(args.SaveTemp) == False:
	subprocess.call(['rm', '-rf', str(CTEMPDIR)]) 
