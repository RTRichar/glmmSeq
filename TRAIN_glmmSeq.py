#!/usr/bin/env python

import subprocess, argparse, sys, time, os

parser = argparse.ArgumentParser(description="")
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('universal optional arguments')
# required
required.add_argument('-i', '--InputFasta', required = True, help = "\n Input file of reference sequences used to train the classifier. Sequence headers should include a unique sequence identifier, NCBI accessions are recomended, preceded only by '>'. Sequences should be contiguous\n")
required.add_argument('-it', '--InputTax', required = True, help = "\n Input taxonomic lineages (a tab-delimited file with the first column containing unique sequence identifiers compatible with the files provided under '-i'). Normally, this would be in Metaxa2-compatible format.\n")
required.add_argument('-odb', '--OutputDB', required = True, help = "\n Name of database to be produced by TRAIN_glmmSeq. This is the name of the directory that will hold the R GLMMs, fasta and taxonomy files of the resulting database. It will be written in /DBs/, within the directory that holds the glmmSeq executables.\n")
#optional
optional.add_argument('-t', '--Threads', required = False, default = 1, help = "\n Number of processors for Vsearch alignment\n")
optional.add_argument('--SaveTemp', default=False, type=lambda x: (str(x).lower() == 'true'), help = "\n Option for saving intermediate files produced during training (e.g. fasta k-fold partitions and csv formatted cross validation resulst used for glmm fitting). Must be 'True/False,' not 'TRUE/FALSE' or 'T/F.'\n")
optional.add_argument('-id', '--idCutoffs', required = False, default = 1, help = "\n Option for specifying the minimum percent identity of Vsearch alignment matches to be used for glmm fitting. This is adjustable at each taxonomic rank and consists of a comma-separated list of 7 numbers representing the threshold used from kingdom to species (e.g the default of 50,50,60,60,65,75,85 means alignments of <= 50 percent ID will not be used during glmm modelling of the kigdom and phylum ranks) \n")
optional.add_argument('-pcv', '--ProportionForCV', required = False, default = 0.2, help = "\n Size of the reference data partition used for cross-validation")
#### change -pcv  to -k 
optional.add_argument('-hr', '--HighestRank', required = False, type=str, default = 'Class', help = "\n Specifies the lowest resolution rank to be classified (Kingdom, Phylum, or Class). Default = Class")
optional.add_argument('-re', '--reStructure', required = False, type=str, default = 'SpeedGenus', help = "\n Specifies the random effect strucutre to used during modelling. Options include 'Family', 'Genus' and 'SpeedGenus'. 'SpeedGenus' splits cases according to whether they are in species-rich genera. For cases within species rich genera, genus and species models include (1|Order/Family/Genus) as random effect. For non-species-rich genera, random effect term is (1|Class/Order/Family). Default = SpeedGenus")
optional.add_argument('-fe', '--feStructure', required = False, type=str, default = 'CategoricalScndID', help = "\n Specifies fixed effects for the model as either 'CategoricalScndID' or 'ContinuousScndID,' where the percent idenity of the alignmnet to the second highest scoring taxon is either included as a continuous variable or a binned categorical variable. If categorical, the percent identify distance between the first and second highest taxa is calculated. Default = 'CategoricalScndID'")
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
	File.write(args.reStructure+';'+args.HighestRank+';'+args.feStructure+';'+args.idCutoffs+';'+args.gRichness)

# Split into test train
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Partitioning k-fold testing/training sets ###\n')
subprocess.call(['GetTestTrain.py', str(args.InputFasta), str(CTEMPDIR+'/'), str(args.ProportionForCV)])

# Run Vsearch algnmnt
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Running vsearch alignments ###\n\n')
Parts = str(int(1/float(args.ProportionForCV)))
for i in range(0,int(1/float(args.ProportionForCV))): # get files for inferring percent ID of top hit taxon
	subprocess.call(['vsearch', '--usearch_global', str(CTEMPDIR+'/'+str(i)+'_Test.fasta'), '--db', str(CTEMPDIR+'/'+str(i)+'_Train.fasta'), '--id', \
	'0.6', '--maxaccepts', '100', '--maxrejects', '50', '--maxhits', '1', '--gapopen', '0TE', '--gapext', '0TE', '--userout', \
	str(CTEMPDIR+'/'+str(i)+'_CV.vsrch.txt'), '--userfields', 'query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov', '--query_cov', \
	'0.95', '--threads', str(args.Threads)]) # could skip this and pull the first hit from the top 50 VSEARCH run (next line)
for i in range(0,int(1/float(args.ProportionForCV))): # get files for inferring percent ID of second-to-top hit taxon
	subprocess.call(['vsearch', '--usearch_global', str(CTEMPDIR+'/'+str(i)+'_Test.fasta'), '--db', str(CTEMPDIR+'/'+str(i)+'_Train.fasta'), '--id', \
	'0.6', '--maxaccepts', '100', '--maxrejects', '50', '--maxhits', '50', '--gapopen', '0TE', '--gapext', '0TE', '--userout', \
	str(CTEMPDIR+'/'+str(i)+'_CV_2nd.vsrch.txt'), '--userfields', 'query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov', '--query_cov', \
	'0.95', '--threads', str(args.Threads)])

# Get LogReg file
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Formatting vsearch outputs for GLMM modelling ###\n')
subprocess.call(['CombineCVs.py', str(CTEMPDIR), Parts, 'CV.vsrch.txt'])
subprocess.call(['VsearchToMetaxa2.py', '-v', str(CTEMPDIR+'/CV.vsrch.txt'), '-t', args.InputTax, '-o', str(CTEMPDIR+'/CV.mtxa.tax')])

##############
subprocess.call(['CombineCVs.py', str(CTEMPDIR), Parts, 'CV_2nd.vsrch.txt'])
subprocess.call(['Get2ndHitTaxID_TRAIN.py', args.InputTax, str(CTEMPDIR+'/CV_2nd.vsrch.txt'), str(CTEMPDIR+'/CV_2nd.vsrch.csv')])
##############

subprocess.call(['CurateForLogReg.py', str(CTEMPDIR+'/CV.mtxa.tax'), args.InputTax, str(CTEMPDIR+'/CV_2nd.vsrch.csv'), str(CTEMPDIR+'/CV.LogReg.csv')])

# R file to train on LogReg (must save ModGLMM:
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Modelling data with binomial GLMMs ###\n\n')
subprocess.call(['subModelWithGLMM.r', str(CTEMPDIR+'/CV.LogReg.csv'), DBDIR, args.reStructure, args.HighestRank, \
	args.feStructure, args.idCutoffs, args.gRichness, args.sqrt])

# save consensus filtered fasta and tax to db directory under 'DB.fa' and 'DB.tax'
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Writing database files and cleaning up tmp files ###\n\n')
subprocess.call(['TaxFastaConsensus.py', '-it', str(args.InputTax), '-if', str(args.InputFasta), '-ot', str(DBDIR+'/DB.tax'), '-of', str(DBDIR+'/DB.fa')])

# Clean up tmp
if bool(args.SaveTemp) == False:
	subprocess.call(['rm', '-rf', str(CTEMPDIR)]) 
