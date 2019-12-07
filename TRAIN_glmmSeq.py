#!/usr/bin/env python

import subprocess, argparse, sys, time, os

parser = argparse.ArgumentParser(description="")
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('universal optional arguments')
# required
required.add_argument('-i', '--InputFasta', required = True, help = "\n Input file of reference sequences to be extracted (sequence headers should include a unique sequence identifier, NCBI accessions are recomended, preceded only by '>')\n")
required.add_argument('-it', '--InputTax', required = True, help = "\nInput taxonomic lineages (a tab-delimited file with the first column containing unique sequence identifiers compatible with the files provided under '-r' and '-i'). Normally, this would be in Metaxa2-compatible format, however, if lineages come directly from Taxonomizr, you can skip reformatting and declare '-tf True.' See below for more details.\n")
required.add_argument('-odb', '--OutputDB', required = True, help = "\n Name of database to be produced by TRAIN_glmmSeq. This is the name of the directory that will hold the R GLMMs, fasta and taxonomy files of the resulting database. It will be written in the directory that holds the glmmSeq executables.\n")
#optional
optional.add_argument('-t', '--Threads', required = False, default = 1, help = "\n Number of processors for Vsearch alignment\n")
optional.add_argument('--SaveTemp', default=False, type=lambda x: (str(x).lower() == 'true'), help = "\nOption for \n")
optional.add_argument('-id', '--idCutoffs', required = False, default = 1, help = "\n Option for")
optional.add_argument('-pcv', '--ProportionForCV', required = False, default = 0.1, help = "\n Size of the reference data partition used for cross-validation")
optional.add_argument('-hr', '--HighestRank', required = False, type=str, default = 'Class', help = "\n Analyze to King, Phyl, or Class? default Class")
optional.add_argument('-re', '--reStructure', required = False, type=str, default = 'Family', help = "\n reLevels to Fam, Gen or Sp? default Fam")
args = parser.parse_args()

# create temp directory
CTEMPDIR = str('TRAIN_tmp_' + '_'.join(time.ctime(time.time()).replace(':','_').split()[1:]))
subprocess.call(['mkdir', CTEMPDIR])

# create directory to store database
DBDIR = str(os.path.abspath(os.path.dirname(sys.argv[0]))+'/'+args.OutputDB)
subprocess.call(['mkdir', DBDIR])

# set glmmSeq location
glmmSeqDirectory = os.path.abspath(os.path.dirname(sys.argv[0]))

# make file to hold info on highest rank to analyze to and lowest rank for RE specification
InfoFile = str(DBDIR + '/' + 'InfoFile.txt')
with open(InfoFile, 'w') as File:
	File.write(args.reStructure + '\n' + args.HighestRank)

# Split into test train
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Partitioning k-fold testing/training sets ###\n')
subprocess.call(['GetTestTrain.py', str(args.InputFasta), str(CTEMPDIR+'/'), str(args.ProportionForCV)])

# Run Vsearch algnmnt
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Running vsearch alignments ###\n\n')
for i in range(0,10):
	subprocess.call(['vsearch', '--usearch_global', str(CTEMPDIR+'/'+str(i)+'_Test.fasta'), '--db', str(CTEMPDIR+'/'+str(i)+'_Train.fasta'), '--id', \
	'0.75', '--maxaccepts', '100', '--maxrejects', '50', '--maxhits', '5', '--gapopen', '0TE', '--gapext', '0TE', '--userout', \
	str(CTEMPDIR+'/'+str(i)+'_CV.vsrch.txt'), '--userfields', 'query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov', '--query_cov', \
	'0.95', '--threads', str(args.Threads)])

# Get LogReg file
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Formatting vsearch outputs for GLMM modelling ###\n')
subprocess.call(['CombineCVs.py', str(CTEMPDIR)])
subprocess.call(['VsearchToMetaxa2.py', '-v', str(CTEMPDIR+'/CV.vsrch.txt'), '-t', args.InputTax, '-o', str(CTEMPDIR+'/CV.mtxa.tax')])
subprocess.call(['CurateForLogReg.py', str(CTEMPDIR+'/CV.mtxa.tax'), args.InputTax, str(CTEMPDIR+'/CV.LogReg.csv')])

# R file to train on LogReg (must save ModGLMM:
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Modelling data with binomial GLMMs ###\n\n')
subprocess.call(['ModelWithGLMM.r', str(CTEMPDIR+'/CV.LogReg.csv'), DBDIR, args.reStructure, args.HighestRank])

# save consensus filtered fasta and tax to db directory under 'DB.fa' and 'DB.tax'
sys.stderr.write('\n### ' + time.ctime(time.time()) + ': Writing database files and cleaning up tmp files ###\n\n')
subprocess.call(['TaxFastaConsensus.py', '-it', str(args.InputTax), '-if', str(args.InputFasta), '-ot', str(DBDIR+'/DB.tax'), '-of', str(DBDIR+'/DB.fa')])

# Clean up tmp
if bool(args.SaveTemp) == False:
	subprocess.call(['rm', '-rf', str(CTEMPDIR)]) 
